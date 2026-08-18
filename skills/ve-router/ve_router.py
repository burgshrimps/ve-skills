#!/usr/bin/env python3
"""ve-router: assign exactly one routing class per variant x transcript.

Parses legacy SnpEff EFF and standardised ANN annotations, falling back to a
local snpEff binary or the build-aware Ensembl VEP REST API when neither is
present in the input. Emits the routing decision as data: which annotation
was chosen, the rule that fired, and every alternative that was discarded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

# Effect-name -> routing class. Routing is decided from the Effect/Annotation
# name only. FunctionalClass is NEVER used: real challenge data contains a
# STOP_LOST record whose FunctionalClass field reads "MISSENSE".
TRUNCATING_EFF = {"FRAME_SHIFT", "STOP_GAINED", "START_LOST", "STOP_LOST"}
SPLICE_EFF = {"SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR"}
MISSENSE_EFF = {"NON_SYNONYMOUS_CODING"}
NONCODING_EFF = {
    "INTRON", "INTERGENIC", "INTRAGENIC", "UPSTREAM", "DOWNSTREAM", "EXON",
    "UTR_5_PRIME", "UTR_3_PRIME", "UTR_5_DELETED", "UTR_3_DELETED",
    "SYNONYMOUS_CODING", "SYNONYMOUS_STOP", "SYNONYMOUS_START",
    "SPLICE_SITE_REGION", "NON_CODING_EXON", "GENE", "TRANSCRIPT",
    "REGULATION", "MICRO_RNA", "CODON_CHANGE",
}
# Same buckets for the standardised SO-term vocabulary (modern ANN, VEP REST).
TRUNCATING_SO = {"frameshift_variant", "stop_gained", "start_lost", "stop_lost"}
SPLICE_SO = {"splice_donor_variant", "splice_acceptor_variant"}
MISSENSE_SO = {"missense_variant"}
NONCODING_SO = {
    "intron_variant", "intergenic_variant", "upstream_gene_variant",
    "downstream_gene_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "synonymous_variant", "stop_retained_variant", "start_retained_variant",
    "splice_region_variant", "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant", "coding_sequence_variant",
}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
SELECTION_RULE = (
    "Highest annotation impact tier wins (HIGH > MODERATE > LOW > MODIFIER); "
    "ties are broken by the order the annotation appears in the EFF/ANN field "
    "(first transcript listed). The Effect/Annotation name decides the class - "
    "FunctionalClass is read but never used for routing."
)
_BUILD_HINTS = [
    ("grch37", "GRCh37"), ("g1k_v37", "GRCh37"), ("hg19", "GRCh37"), ("b37", "GRCh37"),
    ("grch38", "GRCh38"), ("hg38", "GRCh38"), ("b38", "GRCh38"),
]
_VEP_HOSTS = {"GRCh37": "https://grch37.rest.ensembl.org", "GRCh38": "https://rest.ensembl.org"}
_SNPEFF_DB = {"GRCh37": "GRCh37.75", "GRCh38": "GRCh38.86"}
_EFF_RE = re.compile(r"^([A-Za-z0-9_]+)\((.*)\)$")
_ANN_FIELDS = ["allele", "effect", "impact", "gene", "gene_id", "feature_type",
               "transcript", "biotype", "rank", "hgvs_c", "hgvs_p", "cdna_pos",
               "cds_pos", "aa_pos", "distance", "errors"]

def _parse_eff(value: str) -> list[dict]:
    entries = []
    for idx, raw in enumerate(value.split(",")):
        raw = raw.strip()
        match = _EFF_RE.match(raw)
        if not match:
            entries.append({"effect": raw or "UNPARSEABLE", "impact": None,
                             "functional_class": None, "gene": None, "transcript": None,
                             "raw": raw, "order": idx})
            continue
        name = match.group(1)
        fields = (match.group(2).split("|") + [""] * 11)[:11]
        impact, func_class, _codon, _aa, _aa_len, gene, _biotype, _coding, transcript, _rank, _gt = fields
        entries.append({"effect": name, "impact": impact or None,
                         "functional_class": func_class or None, "gene": gene or None,
                         "transcript": transcript or None, "raw": raw, "order": idx})
    return entries

def _parse_ann(value: str) -> list[dict]:
    entries = []
    for idx, raw in enumerate(value.split(",")):
        raw = raw.strip()
        parts = (raw.split("|") + [""] * len(_ANN_FIELDS))[: len(_ANN_FIELDS)]
        record = dict(zip(_ANN_FIELDS, parts))
        primary_effect = record["effect"].split("&")[0] if record["effect"] else "UNPARSEABLE"
        entries.append({"effect": primary_effect, "impact": record["impact"] or None,
                         "functional_class": None, "gene": record["gene"] or None,
                         "transcript": record["transcript"] or None, "raw": raw, "order": idx})
    return entries

def _classify(effect: str, fmt: str) -> str | None:
    truncating, splice, missense, noncoding = (
        (TRUNCATING_EFF, SPLICE_EFF, MISSENSE_EFF, NONCODING_EFF) if fmt == "legacy"
        else (TRUNCATING_SO, SPLICE_SO, MISSENSE_SO, NONCODING_SO)
    )
    if effect in truncating:
        return "protein_truncating"
    if effect in splice:
        return "splice"
    if effect in missense:
        return "missense"
    if effect in noncoding:
        return "non_coding"
    return None

def _select(entries: list[dict]) -> tuple[dict, list[dict]]:
    selected = min(entries, key=lambda e: (IMPACT_RANK.get(e["impact"], 99), e["order"]))
    discarded = [e for e in entries if e is not selected]
    return selected, discarded

def _slim(entry: dict) -> dict:
    return {k: entry[k] for k in ("effect", "gene", "transcript", "impact", "functional_class")}

def _finalize(variant: dict, entries: list[dict], fmt: str, source: str, build: str | None,
              reason: str | None = None) -> dict:
    if not entries:
        routing = {"annotation_source": source, "build": build, "selection_rule": None,
                   "selected": None, "discarded": [],
                   "unroutable_reason": reason or "no consequence annotation present and no annotator reachable"}
        return {**variant, "gene": None, "transcript": None, "consequence": None,
                "impact": None, "class": "unroutable", "routing": routing}
    selected, discarded = _select(entries)
    cls = _classify(selected["effect"], fmt)
    unroutable_reason = None if cls else f"unrecognized {fmt} effect name: {selected['effect']}"
    routing = {"annotation_source": source, "build": build, "selection_rule": SELECTION_RULE,
               "selected": _slim(selected), "discarded": [_slim(e) for e in discarded],
               "unroutable_reason": unroutable_reason}
    return {**variant, "gene": selected["gene"], "transcript": selected["transcript"],
            "consequence": selected["effect"], "impact": selected["impact"],
            "class": cls or "unroutable", "routing": routing}

def _detect_build(header_lines: list[str], assembly_flag: str | None) -> str | None:
    def normalize(text: str) -> str | None:
        lowered = text.strip().lower()
        for hint, build in _BUILD_HINTS:
            if hint in lowered:
                return build
        return None
    if assembly_flag:
        return normalize(assembly_flag) or assembly_flag.strip()
    for line in header_lines:
        if line.startswith("##reference"):
            found = normalize(line)
            if found:
                return found
    return None

def _parse_info(info_str: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for token in info_str.split(";"):
        if not token:
            continue
        key, sep, value = token.partition("=")
        info[key] = value if sep else "true"
    return info

def _try_snpeff(chrom: str, pos: str, ref: str, alt: str, build: str | None):
    binary = shutil.which("snpEff")
    db = _SNPEFF_DB.get(build or "")
    if not binary or not db:
        return None, "snpEff"
    mini_vcf = ("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    try:
        proc = subprocess.run([binary, "-noStats", "-classic", db], input=mini_vcf,
                               capture_output=True, text=True, timeout=20, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None, "snpEff"
    if proc.returncode != 0:
        return None, "snpEff"
    for line in proc.stdout.splitlines():
        if not line.startswith("#"):
            cols = line.split("\t")
            if len(cols) >= 8 and "EFF=" in cols[7]:
                return _parse_eff(_parse_info(cols[7])["EFF"]), "snpEff"
    return None, "snpEff"

def _try_vep_rest(chrom: str, pos: str, ref: str, alt: str, build: str | None):
    host = _VEP_HOSTS.get(build or "")
    if not host:
        return None, "vep_rest"
    end = int(pos) + max(len(ref) - 1, 0)
    url = f"{host}/vep/human/region/{chrom}:{pos}-{end}/{alt}?content-type=application/json"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, OSError, ValueError, TimeoutError):
        return None, "vep_rest"
    entries = []
    for hit in payload if isinstance(payload, list) else []:
        for tc in hit.get("transcript_consequences", []):
            terms = tc.get("consequence_terms") or ["UNPARSEABLE"]
            entries.append({"effect": terms[0], "impact": tc.get("impact"), "functional_class": None,
                             "gene": tc.get("gene_symbol"), "transcript": tc.get("transcript_id"),
                             "raw": json.dumps(tc), "order": len(entries)})
    return (entries, "vep_rest") if entries else (None, "vep_rest")

def route_one(variant: dict, info: dict[str, str], build: str | None, allow_network: bool) -> dict:
    if "EFF" in info:
        return _finalize(variant, _parse_eff(info["EFF"]), "legacy", "EFF", build)
    if "ANN" in info:
        return _finalize(variant, _parse_ann(info["ANN"]), "standard", "ANN", build)
    entries, source = _try_snpeff(variant["chrom"], variant["pos"], variant["ref"], variant["alt"], build)
    if not entries and allow_network:
        entries, source = _try_vep_rest(variant["chrom"], variant["pos"], variant["ref"], variant["alt"], build)
    fmt = "legacy" if source == "snpEff" else "standard"
    return _finalize(variant, entries or [], fmt, source, build)

def _open_text(path: Path):
    with path.open("rb") as fh:
        magic = fh.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open
    return opener(path, "rt", encoding="utf-8", errors="replace")

def _load_vcf(path: Path, assembly_flag: str | None) -> list[tuple[dict, dict, str | None]]:
    with _open_text(path) as fh:
        lines = fh.readlines()
    header_lines = [line for line in lines if line.startswith("##")]
    column_lines = [line for line in lines if line.startswith("#CHROM")]
    if not column_lines:
        raise ValueError("Input does not look like a VCF file (no #CHROM header found)")
    samples = column_lines[0].rstrip("\n").split("\t")[9:]
    build = _detect_build(header_lines, assembly_flag)
    records = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError(f"Malformed VCF data line (expected >=8 columns): {line[:80]!r}")
        chrom, pos, vid, ref, alt, _qual, _filt, info_str = fields[:8]
        genotypes = {}
        if len(fields) > 9 and samples:
            fmt_keys = fields[8].split(":")
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0
            for sample, sample_field in zip(samples, fields[9:]):
                genotypes[sample] = sample_field.split(":")[gt_idx]
        variant = {"variant_key": f"{chrom}:{pos}:{ref}:{alt}", "chrom": chrom, "pos": int(pos),
                   "ref": ref, "alt": alt, "id": None if vid == "." else vid, "genotypes": genotypes}
        records.append((variant, _parse_info(info_str), build))
    if not records:
        raise ValueError("VCF contains no data records")
    return records

def _load_json(path: Path) -> list[tuple[dict, dict, str | None]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, dict):
        data = data.get("records", [data])
    if not isinstance(data, list) or not data:
        raise ValueError("JSON input must be a non-empty array of Contract A records")
    records = []
    for entry in data:
        if not isinstance(entry, dict):
            raise ValueError("Each JSON record must be an object")
        missing = {"chrom", "pos", "ref", "alt"} - entry.keys()
        if missing:
            raise ValueError(f"JSON record missing required fields: {sorted(missing)}")
        info = entry.get("info", {}) or {}
        info = _parse_info(info) if isinstance(info, str) else info
        variant = {k: v for k, v in entry.items() if k != "info"}
        variant.setdefault("variant_key", f"{entry['chrom']}:{entry['pos']}:{entry['ref']}:{entry['alt']}")
        records.append((variant, info, entry.get("build")))
    return records

def load_variants(path: Path, assembly_flag: str | None) -> list[tuple[dict, dict, str | None]]:
    if not path.exists():
        raise ValueError(f"Input file not found: {path}")
    with _open_text(path) as fh:
        head = fh.read(4096)
    stripped = head.lstrip()
    if stripped.startswith("##fileformat=VCF") or "\n#CHROM\t" in head or head.startswith("#CHROM\t"):
        return _load_vcf(path, assembly_flag)
    if stripped.startswith("[") or stripped.startswith("{"):
        return _load_json(path)
    raise ValueError("Input is not a recognised VCF (##fileformat=VCF...) or Contract A JSON file")

def summarize(records: list[dict]) -> dict:
    by_class: dict[str, int] = {}
    discarded_records = discarded_total = 0
    for record in records:
        by_class[record["class"]] = by_class.get(record["class"], 0) + 1
        count = len(record["routing"]["discarded"])
        if count:
            discarded_records += 1
            discarded_total += count
    return {"total": len(records), "by_class": by_class,
            "records_with_discarded_alternatives": discarded_records,
            "total_discarded_annotations": discarded_total}

_CSV_FIELDS = ["variant_key", "chrom", "pos", "ref", "alt", "id", "gene", "transcript",
               "consequence", "impact", "class", "annotation_source", "build",
               "discarded_count", "unroutable_reason"]
_RECORD_KEYS = _CSV_FIELDS[:11]

def _csv_row(r: dict) -> dict:
    row = {k: r.get(k) for k in _RECORD_KEYS}
    row.update(annotation_source=r["routing"]["annotation_source"], build=r["routing"]["build"],
               discarded_count=len(r["routing"]["discarded"]), unroutable_reason=r["routing"]["unroutable_reason"])
    return row

def _write_csv(path: Path, records: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        writer.writeheader()
        writer.writerows(_csv_row(r) for r in records)

def _write_report(path: Path, result: dict, input_path: Path, demo: bool) -> None:
    s = result["summary"]
    lines = ["# ve-router Report", "", f"**Input**: `{input_path}`",
              f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
              f"**Variants routed**: {s['total']}", "", "## Routing class counts", "",
              "| Class | Count |", "|---|---:|"]
    for cls in ("protein_truncating", "splice", "missense", "non_coding", "unroutable"):
        lines.append(f"| {cls} | {s['by_class'].get(cls, 0)} |")
    lines += ["", "## Consequence-selection transparency", "",
              f"{s['records_with_discarded_alternatives']} of {s['total']} variants carried more than one "
              f"variant x transcript annotation; {s['total_discarded_annotations']} alternative annotations "
              "were discarded in favour of the highest-impact-tier pick (see Domain Decisions in SKILL.md). "
              "Every discarded alternative is kept in `result.json` under each record's `routing.discarded` "
              "- nothing is silently collapsed.", "",
              "| Variant | Chosen consequence | Gene | Transcript | Impact | Discarded alternatives |",
              "|---|---|---|---|---|---:|"]
    for r in result["records"]:
        lines.append(f"| {r['variant_key']} | {r.get('consequence') or '-'} | {r.get('gene') or '-'} | "
                      f"{r.get('transcript') or '-'} | {r.get('impact') or '-'} | "
                      f"{len(r['routing']['discarded'])} |")
    unroutable = [r for r in result["records"] if r["class"] == "unroutable"]
    if unroutable:
        lines += ["", "## Unroutable variants", ""]
        for r in unroutable:
            lines.append(f"- `{r['variant_key']}`: {r['routing']['unroutable_reason']}")
    lines += ["", DISCLAIMER, ""]
    path.write_text("\n".join(lines), encoding="utf-8")

def write_outputs(result: dict, input_path: Path, output_dir: Path, command: list[str], demo: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _write_csv(output_dir / "tables" / "results.csv", result["records"])
    _write_report(output_dir / "report.md", result, input_path, demo)
    (output_dir / "reproducibility" / "commands.sh").write_text(
        "#!/usr/bin/env bash\n" + " ".join(command) + "\n", encoding="utf-8")

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="ve-router: one routing class per variant x transcript")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("ve_router_out"))
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--assembly", type=str, default=None,
                         help="Override build detection, e.g. GRCh37 or GRCh38")
    args = parser.parse_args(argv)
    input_path = SKILL_DIR / "demo_input.txt" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")
    try:
        raw_records = load_variants(input_path, args.assembly)
        routed = [route_one(v, i, b, allow_network=not args.demo) for v, i, b in raw_records]
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    result = {"skill": "ve-router", "summary": summarize(routed), "records": routed, "disclaimer": DISCLAIMER}
    write_outputs(result, input_path, args.output, [sys.executable, __file__, *sys.argv[1:]], args.demo)
    print(f"ve-router wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
