#!/usr/bin/env python3
"""ve-frequency — a frequency provenance gate for the ve-skills pipeline.

This is NOT another annotator. It answers three questions and refuses to emit a
rarity call unless all three have answers:

  1. Is there a usable population-frequency layer at all?
  2. In which reference build?
  3. Is any AF-looking field actually *cohort* frequency in disguise?

Contract A in, Contract A out with `.freq` populated:

    "freq": {"af": ..., "source": ..., "build": ..., "class": ...,
             "provenance_warning": ...}

    freq.class in RARE | COMMON | NO_DATA

The challenge VCF's INFO/AF is GATK cohort allele frequency over 4 samples, not
population frequency. Reading it as population AF gives a confidently wrong
answer on every variant. This skill names that instead of using it.
"""

import argparse
import gzip
import json
import re
import sys
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
VERSION = "0.1.0"
CHAIN_PATH = SKILL_DIR / "data" / "GRCh37_to_GRCh38.chain.gz"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

FREQ_CLASSES = ("RARE", "COMMON", "NO_DATA")

# INFO keys that genuinely carry population allele frequency. `AF` is deliberately
# absent: in a GATK callset that is cohort frequency over the samples in the file.
POPULATION_AF_KEYS = (
    "AF_TGP", "AF_EXAC", "AF_ESP", "AF_grpmax", "AF_popmax",
    "gnomAD_AF", "gnomad_af", "GNOMAD_AF", "MAX_AF",
)
COHORT_AF_KEYS = ("AF", "MLEAF")

# Rarity boundary applied to a population AF. Deliberately conservative: this
# gate exists to avoid over-claiming, and 1% is the widest defensible "rare".
RARE_MAX_AF = 0.01

# A site counts as looked-at when this fraction of reference samples reached 20x.
MIN_COVERAGE_OVER_20 = 0.80
MIN_ALLELE_NUMBER = 10_000

CONTIG_LENGTHS = {
    "GRCh37": {"1": 249250621, "2": 243199373, "20": 63025520, "X": 155270560},
    "GRCh38": {"1": 248956422, "2": 242193529, "20": 64444167, "X": 156040895},
}
REFERENCE_PATTERNS = {
    "GRCh37": (r"hg19", r"grch37", r"\bb37\b", r"human_g1k_v37"),
    "GRCh38": (r"hg38", r"grch38", r"\bb38\b"),
}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------------------
# Question 3: is this AF actually population frequency?
# ---------------------------------------------------------------------------


def assess_af_provenance(info: dict, sample_count: int | None = None) -> dict:
    """Decide whether any AF-looking INFO field is usable as population frequency.

    Returns the field actually usable, or refuses and says why. A cohort AF is
    named explicitly rather than silently ignored, because "no frequency" and
    "a frequency that means something else" are different problems for the
    downstream reader.
    """
    for key in POPULATION_AF_KEYS:
        if key in info and info[key] not in (None, "", "."):
            try:
                value = float(str(info[key]).split(",")[0])
            except ValueError:
                continue
            return {
                "usable_population_af": True,
                "af": value,
                "source": key,
                "provenance_warning": None,
            }

    present_cohort = [k for k in COHORT_AF_KEYS if k in info]
    if present_cohort:
        samples = f"{sample_count} samples" if sample_count else "the cohort"
        return {
            "usable_population_af": False,
            "af": None,
            "source": None,
            "provenance_warning": (
                f"INFO/{present_cohort[0]} present but is GATK cohort AF over "
                f"{samples}, not population AF"
            ),
        }

    return {
        "usable_population_af": False,
        "af": None,
        "source": None,
        "provenance_warning": "no population allele-frequency field in INFO",
    }


# ---------------------------------------------------------------------------
# Question 2: which build?
# ---------------------------------------------------------------------------


def detect_build(header_lines) -> str | None:
    """Identify the reference build from VCF header lines, or return None.

    Returns None rather than a default. Downstream, an unknown build forces
    NO_DATA: 'match the build explicitly or declare that you could not'.
    """
    contig_votes: set[str] = set()
    reference_votes: set[str] = set()

    for line in header_lines:
        if line.startswith("##contig"):
            match = re.search(r"ID=([^,>]+).*?length=(\d+)", line)
            if match:
                contig = normalise_contig(match.group(1))
                length = int(match.group(2))
                for build, lengths in CONTIG_LENGTHS.items():
                    if lengths.get(contig) == length:
                        contig_votes.add(build)
        elif line.startswith("##reference"):
            lowered = line.lower()
            for build, patterns in REFERENCE_PATTERNS.items():
                if any(re.search(p, lowered) for p in patterns):
                    reference_votes.add(build)

    if len(contig_votes) > 1:
        return None
    if contig_votes and reference_votes and contig_votes != reference_votes:
        return None
    votes = contig_votes or reference_votes
    return votes.pop() if len(votes) == 1 else None


def normalise_contig(name: str) -> str:
    contig = str(name).strip()
    if contig.lower().startswith("chr"):
        contig = contig[3:]
    if contig in {"M", "m"}:
        return "MT"
    return contig.upper() if contig.lower() in {"x", "y", "mt"} else contig


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Liftover: our data is b37, gnomAD v4 is GRCh38-only
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MappedCoord:
    chrom: str
    pos: int
    strand: str


@dataclass
class _ChainBlock:
    t_start: int
    t_end: int
    q_name: str
    q_start: int
    q_strand: str
    q_size: int


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


class ChainMapper:
    """Minimal chain reader — stdlib only, no pyliftover/CrossMap dependency."""

    def __init__(self, blocks_by_chrom):
        self._blocks = {}
        self._starts = {}
        for chrom, blocks in blocks_by_chrom.items():
            ordered = sorted(blocks, key=lambda b: b.t_start)
            self._blocks[chrom] = ordered
            self._starts[chrom] = [b.t_start for b in ordered]

    @classmethod
    def from_file(cls, path: Path) -> "ChainMapper":
        blocks: dict = {}
        t_name = q_name = q_strand = None
        t_pos = q_pos = q_size = 0
        with _open_text(Path(path)) as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("chain"):
                    parts = line.split()
                    t_name = normalise_contig(parts[2])
                    t_pos = int(parts[5])
                    q_name = normalise_contig(parts[7])
                    q_size = int(parts[8])
                    q_strand = parts[9]
                    q_pos = int(parts[10])
                    continue
                if t_name is None:
                    continue
                fields = line.split()
                size = int(fields[0])
                blocks.setdefault(t_name, []).append(
                    _ChainBlock(t_pos, t_pos + size, q_name, q_pos, q_strand, q_size)
                )
                if len(fields) >= 3:
                    t_pos += size + int(fields[1])
                    q_pos += size + int(fields[2])
                else:
                    t_name = None
        return cls(blocks)

    def map(self, chrom: str, pos: int) -> MappedCoord | None:
        contig = normalise_contig(chrom)
        blocks = self._blocks.get(contig)
        if not blocks:
            return None
        zero_based = int(pos) - 1
        index = bisect_right(self._starts[contig], zero_based) - 1
        if index < 0:
            return None
        block = blocks[index]
        if not (block.t_start <= zero_based < block.t_end):
            return None
        offset = zero_based - block.t_start
        q_offset = block.q_start + offset
        mapped = block.q_size - 1 - q_offset if block.q_strand == "-" else q_offset
        return MappedCoord(chrom=block.q_name, pos=mapped + 1, strand=block.q_strand)


@dataclass
class LiftResult:
    status: str
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    strand: str | None = None
    reason: str | None = None


def lift_variant(chrom: str, pos: int, ref: str, alt: str, mapper: ChainMapper) -> LiftResult:
    """Lift one variant; 32.5% of chain blocks are minus-strand and need revcomp."""
    mapped = mapper.map(chrom, pos)
    if mapped is None:
        return LiftResult(status="unresolvable",
                          reason=f"{normalise_contig(chrom)}:{pos} not covered by the chain")
    if mapped.strand == "-":
        return LiftResult(status="lifted", chrom=mapped.chrom, pos=mapped.pos,
                          ref=reverse_complement(ref), alt=reverse_complement(alt), strand="-")
    return LiftResult(status="lifted", chrom=mapped.chrom, pos=mapped.pos,
                      ref=ref, alt=alt, strand="+")


# ---------------------------------------------------------------------------
# Question 1: is there a usable frequency layer, and what does it say?
# ---------------------------------------------------------------------------


def classify_lookup(record: dict | None) -> dict:
    """Turn a reference-database record into a freq class.

    Absence only counts as evidence of rarity when the site was actually
    covered. Absence at an uncovered site is NO_DATA, never RARE.
    """
    if not record:
        return {"class": "NO_DATA", "upper_bound_af": None,
                "reason": "variant not present in the reference dataset"}

    allele_count = record.get("ac") or 0
    allele_number = record.get("an") or 0
    coverage = record.get("coverage_over_20")

    if allele_count > 0 and allele_number:
        af = allele_count / allele_number
        return {
            "class": "RARE" if af < RARE_MAX_AF else "COMMON",
            "upper_bound_af": None,
            "reason": f"observed at AF {af:.3g}",
        }

    well_covered = (
        coverage is not None
        and coverage >= MIN_COVERAGE_OVER_20
        and allele_number >= MIN_ALLELE_NUMBER
    )
    if well_covered:
        return {
            "class": "RARE",
            "upper_bound_af": 3.0 / allele_number,
            "reason": (f"not observed in {allele_number:,} called alleles at a "
                       f"well-covered site; bounded above by ~3/AN"),
        }
    return {
        "class": "NO_DATA",
        "upper_bound_af": None,
        "reason": "not observed, but the site was not adequately covered",
    }


def resolve_freq(variant_key, info, build, sample_count=None, mapper=None, source=None) -> dict:
    """Produce the Contract A `.freq` block for one variant.

    Order matters: provenance is checked before any number is emitted. A missing
    build short-circuits to NO_DATA even when an AF-looking field exists.
    """
    provenance = assess_af_provenance(info or {}, sample_count=sample_count)

    if build is None:
        warning = "reference build could not be determined from the VCF header"
        if provenance["provenance_warning"]:
            warning = f"{warning}; {provenance['provenance_warning']}"
        return {"af": None, "source": None, "build": None,
                "class": "NO_DATA", "provenance_warning": warning}

    if provenance["usable_population_af"]:
        af = provenance["af"]
        return {
            "af": af,
            "source": provenance["source"],
            "build": build,
            "class": "RARE" if af < RARE_MAX_AF else "COMMON",
            "provenance_warning": None,
        }

    # No usable in-file frequency. A live reference lookup may still answer,
    # but only if the coordinates can be placed in that dataset's build.
    if source is not None:
        lifted_key = _lifted_key(variant_key, build, mapper)
        if lifted_key is None:
            return {"af": None, "source": None, "build": build, "class": "NO_DATA",
                    "provenance_warning": (
                        f"{provenance['provenance_warning']}; and coordinates could not be "
                        "lifted to the reference dataset's build")}
        verdict = classify_lookup(source.fetch(lifted_key))
        return {
            "af": None,
            "source": getattr(source, "source_label", "reference lookup"),
            "build": build,
            "class": verdict["class"],
            "provenance_warning": (
                provenance["provenance_warning"]
                if verdict["class"] == "NO_DATA" else None
            ),
        }

    return {"af": None, "source": None, "build": build, "class": "NO_DATA",
            "provenance_warning": provenance["provenance_warning"]}


def _lifted_key(variant_key: str, build: str, mapper) -> str | None:
    """Translate a b37 variant key into the reference dataset's b38 space."""
    try:
        chrom, pos, ref, alt = str(variant_key).split(":")
    except ValueError:
        return None
    if build == "GRCh38":
        return f"{normalise_contig(chrom)}-{pos}-{ref}-{alt}"
    if mapper is None:
        return None
    lift = lift_variant(chrom, int(pos), ref, alt, mapper)
    if lift.status != "lifted":
        return None
    return f"{lift.chrom}-{lift.pos}-{lift.ref}-{lift.alt}"


# ---------------------------------------------------------------------------
# Contract A record handling
# ---------------------------------------------------------------------------


def parse_info(raw_info: str) -> dict:
    info = {}
    for item in (raw_info or "").split(";"):
        if "=" in item:
            key, _, value = item.partition("=")
            info[key] = value
        elif item:
            info[item] = True
    return info


def annotate_records(records, build, sample_count=None, mapper=None, source=None):
    """Populate `.freq` on Contract A records, leaving every other field alone."""
    for record in records:
        info = record.get("info")
        if not isinstance(info, dict):
            info = parse_info(record.get("raw_info", ""))
        record["freq"] = resolve_freq(
            variant_key=record.get("variant_key"), info=info, build=build,
            sample_count=sample_count, mapper=mapper, source=source,
        )
    return records


def read_vcf_as_contract_a(path: Path):
    """Bootstrap Contract A records straight from a VCF (when run standalone)."""
    header, records, samples = [], [], []
    with _open_text(Path(path)) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header.append(line)
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alts = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            genotypes = {}
            if len(parts) > 9 and samples:
                for name, field in zip(samples, parts[9:]):
                    genotypes[name] = field.split(":")[0]
            for alt in alts.split(","):
                if alt in {".", "<NON_REF>"}:
                    continue
                records.append({
                    "variant_key": f"{normalise_contig(chrom)}:{pos}:{ref}:{alt}",
                    "chrom": normalise_contig(chrom), "pos": pos,
                    "ref": ref, "alt": alt, "id": vid,
                    "genotypes": genotypes,
                    "raw_info": parts[7],
                })
    return header, records, len(samples)


def load_input(path: Path):
    """Accept either a Contract A JSON array or a VCF."""
    text_head = ""
    with _open_text(Path(path)) as handle:
        text_head = handle.read(2048)
    if text_head.lstrip().startswith(("[", "{")):
        with _open_text(Path(path)) as handle:
            payload = json.load(handle)
        records = payload["records"] if isinstance(payload, dict) else payload
        header = payload.get("header", []) if isinstance(payload, dict) else []
        samples = payload.get("sample_count") if isinstance(payload, dict) else None
        if samples is None:
            samples = len(records[0].get("genotypes", {})) if records else None
        return header, records, samples
    return read_vcf_as_contract_a(path)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def write_outputs(records, build, output_dir: Path, input_label: str):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    counts = {cls: 0 for cls in FREQ_CLASSES}
    for record in records:
        counts[record["freq"]["class"]] += 1

    warnings: dict = {}
    for record in records:
        warning = record["freq"].get("provenance_warning")
        if warning:
            warnings[warning] = warnings.get(warning, 0) + 1

    payload = {
        "skill": "ve-frequency",
        "version": VERSION,
        "input": input_label,
        "build": build,
        "record_count": len(records),
        "class_counts": counts,
        "provenance_warnings": warnings,
        "records": records,
        "scope_note": (
            "ve-frequency is a provenance gate, not an annotator. RARE and COMMON are "
            "emitted only when a population-frequency layer and a reference build are "
            "both established. NO_DATA means the gate refused, not that the variant is "
            "absent from the population."
        ),
        "disclaimer": DISCLAIMER,
    }
    (output_dir / "result.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# ve-frequency — frequency provenance gate",
        "",
        f"- **Input**: {input_label}",
        f"- **Build**: {build or 'UNDETERMINED'}",
        f"- **Records**: {len(records)}",
        "",
        "## Frequency classes",
        "",
        "| Class | Count | Meaning |",
        "|---|---|---|",
        f"| RARE | {counts['RARE']} | Population frequency established and below {RARE_MAX_AF:g} |",
        f"| COMMON | {counts['COMMON']} | Population frequency established and at or above {RARE_MAX_AF:g} |",
        f"| NO_DATA | {counts['NO_DATA']} | **Provenance could not be established — no rarity claim made** |",
        "",
    ]
    if warnings:
        lines += ["## Why the gate refused", "", "| Count | Reason |", "|---|---|"]
        for reason, count in sorted(warnings.items(), key=lambda kv: -kv[1]):
            lines.append(f"| {count} | {reason} |")
        lines.append("")

    lines += ["## Records", "", "| Variant | Class | AF | Source |", "|---|---|---|---|"]
    for record in records[:50]:
        freq = record["freq"]
        af = "—" if freq["af"] is None else f"{freq['af']:.3g}"
        lines.append(f"| `{record.get('variant_key')}` | {freq['class']} | {af} "
                     f"| {freq['source'] or '—'} |")
    if len(records) > 50:
        lines.append(f"| … | *{len(records) - 50} further records in result.json* | | |")

    lines += ["", "## Scope", "", payload["scope_note"], "", "---", "", f"*{DISCLAIMER}*", ""]
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    return payload


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="ve-frequency: frequency provenance gate (Contract A in, .freq out)")
    parser.add_argument("--input", type=Path, help="Contract A JSON or VCF")
    parser.add_argument("--output", type=Path, default=Path("ve_frequency_out"))
    parser.add_argument("--demo", action="store_true", help="Run offline on bundled demo data")
    parser.add_argument("--assembly", type=str, default=None,
                        help="Override build detection (GRCh37/GRCh38)")
    args = parser.parse_args(argv)

    input_path = DEMO_INPUT if args.demo else args.input
    if input_path is None:
        print("No --input given. Try --demo to see the gate refuse on cohort-AF-only data:\n"
              "  python skills/ve-frequency/ve_frequency.py --demo --output /tmp/ve_freq_demo",
              file=sys.stderr)
        return 2

    header, records, sample_count = load_input(input_path)
    build = args.assembly or detect_build(header)
    mapper = ChainMapper.from_file(CHAIN_PATH) if (build == "GRCh37" and CHAIN_PATH.exists()) else None

    annotate_records(records, build=build, sample_count=sample_count,
                     mapper=mapper, source=None)
    payload = write_outputs(records, build, args.output, str(input_path))

    print(f"ve-frequency: {payload['record_count']} records, build={build or 'UNDETERMINED'}")
    for cls in FREQ_CLASSES:
        print(f"  {cls:<8} {payload['class_counts'][cls]}")
    for reason, count in payload["provenance_warnings"].items():
        print(f"  refused ({count}): {reason}")
    print(f"Report: {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
