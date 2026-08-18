#!/usr/bin/env python3
"""ve-splice - ref-vs-alt splice-site strength delta for Contract A splice records.

Scores the CHANGE in splice-site strength between the reference and the alternate
allele using the MaxEntScan maximum-entropy model (Yeo & Burge 2004). Unlike
site-detection tools, which ask "is there a splice site in this sequence?", this
skill always compares two alleles and reports the delta.

Standard library only. No network, no GPU, no model inference.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path
from typing import Iterable

SKILL = "ve-splice"
VERSION = "0.1.0"

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and does "
    "not provide clinical diagnoses. Consult a healthcare professional before making "
    "any medical decisions."
)

# ---------------------------------------------------------------------------
# MaxEntScan model
#
# Ported from maxentpy (https://github.com/kepbod/maxentpy, MIT), which is itself
# a port of the MaxEntScan perl scripts by Gene Yeo and Christopher Burge.
#   Yeo G, Burge C. Maximum entropy modeling of short sequence motifs with
#   applications to RNA splicing signals. J Comput Biol 2004; 11:377-94.
# The bundled score matrices are the maxentpy data files; see data/MAXENTPY-LICENSE.txt.
# ---------------------------------------------------------------------------

BGD = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
CONS1_5 = {"A": 0.004, "C": 0.0032, "G": 0.9896, "T": 0.0032}
CONS2_5 = {"A": 0.0034, "C": 0.0039, "G": 0.0042, "T": 0.9884}
CONS1_3 = {"A": 0.9903, "C": 0.0032, "G": 0.0034, "T": 0.0030}
CONS2_3 = {"A": 0.0027, "C": 0.0037, "G": 0.9905, "T": 0.0030}

DONOR_LEN = 9  # 3 exonic + GT + 4 intronic; the GT sits at index 3..4
DONOR_GT_OFFSET = 3
ACCEPTOR_LEN = 23  # 20 intronic + AG + 3 exonic; the AG sits at index 18..19
ACCEPTOR_AG_OFFSET = 18

_MATRIX5: dict[str, float] | None = None
_MATRIX3: dict[int, dict[int, float]] | None = None


def load_matrix5() -> dict[str, float]:
    global _MATRIX5
    if _MATRIX5 is None:
        m: dict[str, float] = {}
        with gzip.open(DATA / "score5_matrix.txt.gz", "rt") as fh:
            for line in fh:
                key, val = line.split()
                m[key] = float(val)
        _MATRIX5 = m
    return _MATRIX5


def load_matrix3() -> dict[int, dict[int, float]]:
    global _MATRIX3
    if _MATRIX3 is None:
        m: dict[int, dict[int, float]] = {}
        with gzip.open(DATA / "score3_matrix.txt.gz", "rt") as fh:
            for line in fh:
                n, k, s = line.split()
                m.setdefault(int(n), {})[int(k)] = float(s)
        _MATRIX3 = m
    return _MATRIX3


def _hashseq(seq: str) -> int:
    table = str.maketrans("ACGT", "0123")
    digits = seq.translate(table)
    return sum(int(d) * 4 ** (len(digits) - i - 1) for i, d in enumerate(digits))


def score5(seq: str) -> float:
    """MaxEnt score of a 9-mer donor (5' splice site) window, in bits."""
    if len(seq) != DONOR_LEN:
        raise ValueError(f"donor window must be {DONOR_LEN} nt, got {len(seq)}")
    seq = seq.upper()
    key = seq[3:5]
    score = CONS1_5[key[0]] * CONS2_5[key[1]] / (BGD[key[0]] * BGD[key[1]])
    rest = seq[:3] + seq[5:]
    return math.log(score * load_matrix5()[rest], 2)


def score3(seq: str) -> float:
    """MaxEnt score of a 23-mer acceptor (3' splice site) window, in bits."""
    if len(seq) != ACCEPTOR_LEN:
        raise ValueError(f"acceptor window must be {ACCEPTOR_LEN} nt, got {len(seq)}")
    seq = seq.upper()
    key = seq[18:20]
    score = CONS1_3[key[0]] * CONS2_3[key[1]] / (BGD[key[0]] * BGD[key[1]])
    rest = seq[:18] + seq[20:]
    m = load_matrix3()
    s = 1.0
    s *= m[0][_hashseq(rest[:7])]
    s *= m[1][_hashseq(rest[7:14])]
    s *= m[2][_hashseq(rest[14:])]
    s *= m[3][_hashseq(rest[4:11])]
    s *= m[4][_hashseq(rest[11:18])]
    s /= m[5][_hashseq(rest[4:7])]
    s /= m[6][_hashseq(rest[7:11])]
    s /= m[7][_hashseq(rest[11:14])]
    s /= m[8][_hashseq(rest[14:18])]
    return math.log(score * s, 2)


# ---------------------------------------------------------------------------
# Domain decisions - every threshold here is documented in SKILL.md
# ---------------------------------------------------------------------------

SEARCH_RADIUS = 20          # nt either side of the variant to look for a canonical site
REGION_PAD = 80             # nt of reference sequence fetched either side of the variant
MIN_REF_SITE_BITS = 0.0     # a reference window must score above this to be a credible site
DELTA_DAMAGING_BITS = -2.0  # delta at or below this counts as a weakened site
DELTA_GAIN_BITS = 2.0       # delta at or above this counts as a strengthened site
PCT_DAMAGING = -20.0        # and the relative drop must be at least this large
SCORE_SCALE_BITS = 12.0     # |delta| mapped onto 0..1 by dividing by this
CONF_HIGH_REF_BITS = 3.0    # reference site this strong -> high confidence in the delta

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Sequence providers
# ---------------------------------------------------------------------------


class FastaProvider:
    """Random access to an indexed FASTA using its .fai, stdlib only."""

    def __init__(self, path: Path):
        self.path = Path(path)
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(
                f"missing FASTA index {fai}; create it with 'samtools faidx {self.path}'"
            )
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with fai.open() as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split()[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fh = self.path.open("rb")

    def contigs(self) -> set[str]:
        return set(self.index)

    def fetch(self, chrom: str, start: int, end: int) -> str | None:
        """1-based inclusive fetch. Returns None if the contig or range is unavailable."""
        key = chrom if chrom in self.index else _alt_contig(chrom, self.index)
        if key is None:
            return None
        length, offset, linebases, linewidth = self.index[key]
        start = max(1, start)
        end = min(length, end)
        if end < start:
            return None
        s0 = start - 1
        byte_start = offset + (s0 // linebases) * linewidth + (s0 % linebases)
        e0 = end - 1
        byte_end = offset + (e0 // linebases) * linewidth + (e0 % linebases)
        self._fh.seek(byte_start)
        raw = self._fh.read(byte_end - byte_start + 1)
        return raw.decode("ascii", "replace").replace("\n", "").replace("\r", "").upper()

    def close(self) -> None:
        self._fh.close()


def _alt_contig(chrom: str, index: dict) -> str | None:
    """Tolerate chr-prefix mismatches between the VCF and the reference."""
    for cand in (chrom, chrom.removeprefix("chr"), "chr" + chrom):
        if cand in index:
            return cand
    return None


class BundledProvider:
    """Pre-extracted reference windows, so --demo runs with no FASTA and no network."""

    def __init__(self, path: Path):
        payload = json.loads(Path(path).read_text())
        self.build = payload.get("build")
        self.regions = payload["regions"]  # variant_key -> {chrom,start,end,seq}

    def _slice(self, r: dict, start: int, end: int) -> str | None:
        if start < r["start"] or end > r["end"]:
            return None
        off = start - r["start"]
        return r["seq"][off : off + (end - start + 1)].upper()

    def fetch_for(self, variant_key: str, chrom: str, start: int, end: int) -> str | None:
        r = self.regions.get(variant_key)
        if r is not None:
            return self._slice(r, start, end)
        # Fall back to any bundled window that covers the requested span, so a record
        # whose key differs from the one used at extraction time still resolves.
        for cand in self.regions.values():
            if str(cand["chrom"]) == str(chrom):
                got = self._slice(cand, start, end)
                if got:
                    return got
        return None


# ---------------------------------------------------------------------------
# Core scoring
# ---------------------------------------------------------------------------


def _kind_from_consequence(consequence: str | None) -> str | None:
    c = (consequence or "").upper()
    if "SPLICE_SITE_DONOR" in c or "SPLICE_DONOR" in c:
        return "donor"
    if "SPLICE_SITE_ACCEPTOR" in c or "SPLICE_ACCEPTOR" in c:
        return "acceptor"
    return None


def _candidate_sites(oriented: str, var_lo: int, var_hi: int, kind: str) -> list[dict]:
    """Find canonical GT/AG anchors in a transcript-oriented sequence.

    var_lo/var_hi bracket the reference allele's footprint in `oriented` coordinates.
    Only windows that actually contain the variant are returned - a window that does
    not overlap the variant cannot change between ref and alt.
    """
    dinuc, win_len, anchor_off = (
        ("GT", DONOR_LEN, DONOR_GT_OFFSET) if kind == "donor" else ("AG", ACCEPTOR_LEN, ACCEPTOR_AG_OFFSET)
    )
    out = []
    for p in range(len(oriented) - 1):
        if oriented[p : p + 2] != dinuc:
            continue
        w0 = p - anchor_off
        w1 = w0 + win_len
        if w0 < 0 or w1 > len(oriented):
            continue
        if var_hi < w0 or var_lo >= w1:
            continue  # window does not overlap the variant
        window = oriented[w0:w1]
        if any(b not in "ACGT" for b in window):
            continue
        out.append({"anchor": p, "w0": w0, "w1": w1, "window": window})
    return out


def _score_window(window: str, kind: str) -> float:
    return score5(window) if kind == "donor" else score3(window)


def _apply(plus: str, i: int, ref: str, alt: str) -> str:
    return plus[:i] + alt + plus[i + len(ref) :]


def equivalent_indel_span(plus: str, vi: int, ref: str, alt: str,
                          reach: int = 12) -> tuple[int, int, int]:
    """Return (span_start, span_end, n_representations) for an indel in a repeat.

    An indel inside a repeat can be written at several positions that all produce the
    same alternate sequence. VCF left-alignment is done on the plus strand, so for a
    minus-strand transcript the chosen representation is arbitrary with respect to the
    splice site. Anything inside this span cannot be placed exactly.
    """
    if len(ref) == len(alt):
        return (vi, vi + len(ref) - 1, 1)
    target = _apply(plus, vi, ref, alt)
    hits = [vi]
    for shift in range(-reach, reach + 1):
        j = vi + shift
        if shift == 0 or j < 0 or j + len(ref) > len(plus):
            continue
        # keep the anchor-base convention: same first base, same indel payload
        if len(alt) > len(ref):
            cand = plus[j] + alt[1:]
            if _apply(plus, j, plus[j], cand) == target:
                hits.append(j)
        else:
            cand_ref = plus[j : j + len(ref)]
            if _apply(plus, j, cand_ref, plus[j]) == target:
                hits.append(j)
    return (min(hits), max(hits) + len(ref) - 1, len(hits))


def evaluate_record(rec: dict, fetch_region, matrices_ready: bool = True) -> dict:
    """Score one Contract A record, returning a Contract B result."""
    vkey = rec.get("variant_key") or f"{rec.get('chrom')}:{rec.get('pos')}:{rec.get('ref')}:{rec.get('alt')}"
    out = {
        "skill": SKILL,
        "variant_key": vkey,
        "score": None,
        "direction": None,
        "confidence": None,
        "in_domain": False,
        "abstain_reason": None,
        "evidence": {},
    }

    if rec.get("class") != "splice":
        out["abstain_reason"] = (
            f"ve-splice scores changes in splice-site strength at annotated splice sites; "
            f"this record was routed as class={rec.get('class')!r}, which is outside that domain."
        )
        return out

    kind = _kind_from_consequence(rec.get("consequence"))
    if kind is None:
        out["abstain_reason"] = (
            f"consequence {rec.get('consequence')!r} does not identify the site as a donor "
            f"or an acceptor, so the correct MaxEntScan model cannot be selected."
        )
        return out

    chrom, pos = str(rec.get("chrom")), int(rec.get("pos"))
    ref, alt = str(rec.get("ref", "")).upper(), str(rec.get("alt", "")).upper()

    region_start = pos - REGION_PAD
    region_end = pos + REGION_PAD
    plus = fetch_region(vkey, chrom, region_start, region_end)
    if not plus:
        out["abstain_reason"] = (
            f"no reference sequence available for {chrom}:{region_start}-{region_end} on the "
            f"declared build, so a ref-vs-alt comparison could not be made."
        )
        return out

    vi = pos - region_start  # index of the variant's first base on the plus strand
    if plus[vi : vi + len(ref)] != ref:
        out["abstain_reason"] = (
            f"reference allele {ref!r} does not match the reference sequence "
            f"{plus[vi:vi + len(ref)]!r} at {chrom}:{pos}; the VCF and the FASTA are not the "
            f"same build, so no comparison was attempted."
        )
        return out

    alt_plus = plus[:vi] + alt + plus[vi + len(ref) :]
    delta_len = len(alt) - len(ref)
    amb_start, amb_end, n_repr = equivalent_indel_span(plus, vi, ref, alt)

    # Search both orientations. The strand is inferred from the sequence rather than
    # taken on trust, and the winning orientation is reported in the evidence.
    win_len = DONOR_LEN if kind == "donor" else ACCEPTOR_LEN
    anchor_off = DONOR_GT_OFFSET if kind == "donor" else ACCEPTOR_AG_OFFSET

    candidates = []
    n = len(plus)
    for strand in ("+", "-"):
        oriented = plus if strand == "+" else revcomp(plus)
        if strand == "+":
            lo, hi = vi, vi + max(len(ref), 1) - 1
        else:
            lo, hi = n - 1 - (vi + max(len(ref), 1) - 1), n - 1 - vi
        for c in _candidate_sites(oriented, lo, hi, kind):
            # Leftmost plus-strand index of the canonical dinucleotide.
            a_plus = c["anchor"] if strand == "+" else n - 2 - c["anchor"]
            if abs(a_plus - vi) > SEARCH_RADIUS:
                continue
            # The alt window is re-anchored on the SAME genomic dinucleotide. An indel
            # shifts everything downstream of the variant, so the anchor moves with it;
            # without this the GT/AG slides out of its fixed position in the window and
            # the model scores a motif that is not the site being tested.
            if a_plus > amb_end:
                shifts = [delta_len]
            elif a_plus + 1 < amb_start:
                shifts = [0]
            else:
                # The site lies inside the span where the indel cannot be placed exactly.
                shifts = sorted({0, delta_len})

            scored_variants = []
            for sh in shifts:
                a_alt = a_plus + sh
                ws = a_alt - anchor_off if strand == "+" else a_alt + anchor_off - win_len + 2
                if ws < 0 or ws + win_len > len(alt_plus):
                    continue
                sl = alt_plus[ws : ws + win_len]
                if any(b not in "ACGT" for b in sl):
                    continue
                w = sl if strand == "+" else revcomp(sl)
                try:
                    scored_variants.append((_score_window(w, kind), w))
                except (KeyError, ValueError):
                    continue
            if not scored_variants:
                continue
            try:
                mes_ref = _score_window(c["window"], kind)
            except (KeyError, ValueError):
                continue
            # Where the representation is ambiguous, report the least dramatic outcome
            # rather than the most, and carry the range so the ambiguity stays visible.
            mes_alt, alt_window = min(scored_variants, key=lambda t: abs(t[0] - mes_ref))
            alt_range = [round(min(s for s, _ in scored_variants), 4),
                         round(max(s for s, _ in scored_variants), 4)]
            candidates.append(
                {
                    "strand": strand,
                    "site_pos": region_start + a_plus,
                    "site_offset_from_variant": a_plus - vi,
                    "ref_window": c["window"],
                    "alt_window": alt_window,
                    "anchor_intact": alt_window[anchor_off : anchor_off + 2]
                    == ("GT" if kind == "donor" else "AG"),
                    "mes_ref": round(mes_ref, 4),
                    "mes_alt": round(mes_alt, 4),
                    "mes_alt_range": alt_range,
                    "representation_ambiguous": len(scored_variants) > 1,
                    "delta_bits": round(mes_alt - mes_ref, 4),
                }
            )

    credible = [c for c in candidates if c["mes_ref"] > MIN_REF_SITE_BITS]
    if not credible:
        out["abstain_reason"] = (
            f"no canonical {'GT' if kind == 'donor' else 'AG'} {kind} motif scoring above "
            f"{MIN_REF_SITE_BITS} bits was found within {SEARCH_RADIUS} nt of {chrom}:{pos} in "
            f"either orientation, so there is no reference site whose strength a delta could be "
            f"measured against."
        )
        out["evidence"] = {
            "site_kind": kind,
            "candidates_examined": len(candidates),
            "note": "the HIGH impact on this record rests on an annotation call, not on a "
                    "motif this model can locate.",
        }
        return out

    # The annotated site is the strongest credible reference motif nearby.
    best = max(credible, key=lambda c: c["mes_ref"])
    others = [c for c in credible if c is not best]

    delta = best["delta_bits"]
    pct = (delta / abs(best["mes_ref"]) * 100.0) if best["mes_ref"] else 0.0

    if delta <= DELTA_DAMAGING_BITS and pct <= PCT_DAMAGING:
        direction = "damaging"
    elif delta >= DELTA_GAIN_BITS:
        direction = "strengthening"
    else:
        direction = "neutral"

    score = min(1.0, abs(delta) / SCORE_SCALE_BITS)

    if best["mes_ref"] >= CONF_HIGH_REF_BITS:
        confidence = "high"
    elif best["mes_ref"] >= MIN_REF_SITE_BITS:
        confidence = "medium"
    else:
        confidence = "low"
    if len(ref) != len(alt):
        # An indel changes the window length; the window is re-read from the mutated
        # sequence at the same anchor, which is a modelling choice, not a measurement.
        confidence = "medium" if confidence == "high" else confidence
    if best["representation_ambiguous"]:
        # The delta depends on which of several equivalent VCF spellings was used.
        confidence = "low"

    flags = []
    if len(ref) != len(alt):
        flags.append("indel: alt window re-anchored on the same genomic GT/AG")
    if not best["anchor_intact"]:
        flags.append(
            f"the canonical {'GT' if kind == 'donor' else 'AG'} dinucleotide itself is "
            f"altered by this variant"
        )
    if best["site_offset_from_variant"] != 0:
        flags.append(
            f"variant sits {abs(best['site_offset_from_variant'])} nt from the "
            f"{'GT' if kind == 'donor' else 'AG'}, inside the scored window but not on it"
        )
    if n_repr > 1:
        flags.append(
            f"indel has {n_repr} equivalent representations in this repeat "
            f"({chrom}:{region_start + amb_start}-{region_start + amb_end}); VCF "
            f"left-alignment is done on the plus strand and is arbitrary with respect "
            f"to this site"
        )
    if best["representation_ambiguous"]:
        lo, hi = best["mes_alt_range"]
        flags.append(
            f"the site falls inside that ambiguous span, so MaxEnt(alt) is only bounded "
            f"to [{lo}, {hi}] bits; the least dramatic outcome is reported"
        )
    if others:
        flags.append(f"{len(others)} other canonical {kind} motif(s) overlapped the variant")

    out.update(
        {
            "score": round(score, 4),
            "direction": direction,
            "confidence": confidence,
            "in_domain": True,
            "evidence": {
                "model": "MaxEntScan (Yeo & Burge 2004), maximum-entropy splice-site model",
                "site_kind": kind,
                "inferred_strand": best["strand"],
                "site_pos": best["site_pos"],
                "site_offset_from_variant": best["site_offset_from_variant"],
                "anchor_intact": best["anchor_intact"],
                "gene": rec.get("gene"),
                "transcript": rec.get("transcript"),
                "consequence": rec.get("consequence"),
                "mes_ref_bits": best["mes_ref"],
                "mes_alt_bits": best["mes_alt"],
                "delta_bits": best["delta_bits"],
                "pct_change": round(pct, 2),
                "ref_window": best["ref_window"],
                "alt_window": best["alt_window"],
                "alternatives_considered": [
                    {k: c[k] for k in ("strand", "mes_ref", "mes_alt", "delta_bits")} for c in others
                ],
                "flags": flags,
                "interpretation_boundary": (
                    "This is a change in predicted motif strength. It is not evidence that "
                    "splicing is altered in any tissue, and this data set contains no RNA "
                    "evidence with which to test it."
                ),
            },
        }
    )
    return out


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def load_records(path: Path) -> list[dict]:
    payload = json.loads(Path(path).read_text())
    if isinstance(payload, dict):
        for key in ("records", "variants", "results"):
            if key in payload:
                return list(payload[key])
        raise ValueError("input JSON object has no 'records' key")
    return list(payload)


def write_outputs(outdir: Path, results: list[dict], summary: dict, commands: list[str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    (outdir / "reproducibility").mkdir(exist_ok=True)

    (outdir / "result.json").write_text(
        json.dumps(
            {"skill": SKILL, "version": VERSION, "summary": summary,
             "results": results, "disclaimer": DISCLAIMER},
            indent=2,
        )
        + "\n"
    )

    with (outdir / "tables" / "results.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["variant_key", "gene", "site_kind", "strand", "mes_ref_bits", "mes_alt_bits",
             "delta_bits", "pct_change", "score", "direction", "confidence", "in_domain",
             "abstain_reason"]
        )
        for r in results:
            e = r.get("evidence") or {}
            w.writerow(
                [r["variant_key"], e.get("gene", ""), e.get("site_kind", ""),
                 e.get("inferred_strand", ""), e.get("mes_ref_bits", ""), e.get("mes_alt_bits", ""),
                 e.get("delta_bits", ""), e.get("pct_change", ""), r.get("score", ""),
                 r.get("direction", ""), r.get("confidence", ""), r["in_domain"],
                 r.get("abstain_reason") or ""]
            )

    (outdir / "reproducibility" / "commands.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands) + "\n"
    )
    (outdir / "report.md").write_text(render_report(results, summary))


def render_report(results: list[dict], summary: dict) -> str:
    scored = [r for r in results if r["in_domain"]]
    abstained = [r for r in results if not r["in_domain"]]
    scored.sort(key=lambda r: r["score"], reverse=True)

    lines = [
        "# ve-splice report",
        "",
        f"Records received: **{summary['n_input']}**  |  scored: **{len(scored)}**  |  "
        f"abstained: **{len(abstained)}**",
        "",
        f"Model: MaxEntScan (Yeo & Burge 2004). Reference: `{summary.get('reference')}` "
        f"(build `{summary.get('build')}`).",
        "",
        "`delta` is MaxEnt(alt) - MaxEnt(ref) in bits. Negative means the motif is weaker "
        "with the alternate allele.",
        "",
        "## Scored",
        "",
    ]
    if scored:
        lines += [
            "| variant | gene | site | strand | MES ref | MES alt | delta | % | direction | conf |",
            "|---|---|---|---|---:|---:|---:|---:|---|---|",
        ]
        for r in scored:
            e = r["evidence"]
            lines.append(
                f"| `{r['variant_key']}` | {e.get('gene') or '-'} | {e['site_kind']} | "
                f"{e['inferred_strand']} | {e['mes_ref_bits']:.2f} | {e['mes_alt_bits']:.2f} | "
                f"{e['delta_bits']:+.2f} | {e['pct_change']:+.0f}% | {r['direction']} | "
                f"{r['confidence']} |"
            )
    else:
        lines.append("_No record could be scored._")

    lines += ["", "## Abstentions", ""]
    if abstained:
        for r in abstained:
            lines.append(f"- `{r['variant_key']}` — {r['abstain_reason']}")
    else:
        lines.append("_None._")

    lines += [
        "",
        "## What this does not establish",
        "",
        "- A MaxEnt delta is a change in **motif strength**, not an observation of splicing.",
        "- This data set carries no RNA evidence, so no prediction here is tested against a",
        "  measured transcript.",
        "- The HIGH impact tier on these records is SnpEff's annotation call. Where no canonical",
        "  motif could be located, that call is not supported by anything this model can measure,",
        "  and the record is abstained on rather than scored.",
        "- No statement about frequency, phase, or clinical significance is made or implied.",
        "",
        "## Attribution",
        "",
        "MaxEntScan model and score matrices: Yeo G, Burge C. *Maximum entropy modeling of short",
        "sequence motifs with applications to RNA splicing signals.* J Comput Biol 2004;11:377-94.",
        "Matrices vendored from maxentpy (MIT); see `data/MAXENTPY-LICENSE.txt`.",
        "",
        DISCLAIMER,
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def run(records: Iterable[dict], fasta: Path | None = None,
        bundled: Path | None = None) -> list[dict]:
    """Importable entrypoint: Contract A records in, Contract B results out."""
    records = list(records)
    if fasta:
        provider = FastaProvider(fasta)

        def fetch(vkey, chrom, start, end):
            return provider.fetch(chrom, start, end)
    elif bundled:
        b = BundledProvider(bundled)

        def fetch(vkey, chrom, start, end):
            return b.fetch_for(vkey, chrom, start, end)
    else:
        raise ValueError("one of fasta= or bundled= is required")

    return [evaluate_record(r, fetch) for r in records]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="ve_splice.py",
        description="Score ref-vs-alt splice-site strength deltas for Contract A splice records.",
    )
    p.add_argument("--input", type=Path, help="Contract A records (JSON array or {records:[...]})")
    p.add_argument("--output", type=Path, default=Path("ve_splice_out"))
    p.add_argument("--reference", type=Path, help="indexed reference FASTA matching the VCF build")
    p.add_argument("--demo", action="store_true",
                   help="run on the bundled demo records with bundled sequence, no network")
    p.add_argument("--build", default="GRCh37", help="declared build of the input coordinates")
    args = p.parse_args(argv)

    if args.demo:
        records = load_records(HERE / "demo_input.txt")
        results = run(records, bundled=HERE / "demo_sequences.json")
        reference = "bundled demo windows (pre-extracted from human_g1k_v37)"
        cmd = f"python {Path(__file__).name} --demo --output {args.output}"
    else:
        if not args.input:
            p.error("--input is required unless --demo is given")
        if not args.reference:
            p.error("--reference is required unless --demo is given; ve-splice cannot "
                    "score a delta without the reference allele's sequence context")
        records = load_records(args.input)
        results = run(records, fasta=args.reference)
        reference = str(args.reference)
        cmd = (f"python {Path(__file__).name} --input {args.input} "
               f"--reference {args.reference} --output {args.output}")

    scored = [r for r in results if r["in_domain"]]
    summary = {
        "n_input": len(records),
        "n_scored": len(scored),
        "n_abstained": len(results) - len(scored),
        "n_damaging": sum(1 for r in scored if r["direction"] == "damaging"),
        "build": args.build,
        "reference": reference,
        "model": "MaxEntScan (Yeo & Burge 2004)",
    }
    write_outputs(args.output, results, summary, [cmd])
    print(f"{SKILL}: {summary['n_scored']} scored, {summary['n_abstained']} abstained "
          f"-> {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
