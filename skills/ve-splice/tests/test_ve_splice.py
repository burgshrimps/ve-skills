"""Tests for ve-splice.

Run:  python -m pytest skills/ve-splice/tests -q
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import ve_splice as v  # noqa: E402
from api import run_demo  # noqa: E402


# --------------------------------------------------------------------------
# The MaxEntScan port must reproduce the published reference values exactly.
# Source: maxentpy doctests, themselves from Yeo & Burge's perl scripts.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("cagGTAAGT", 10.86),
    ("gagGTAAGT", 11.08),
    ("taaATAAGT", -0.12),
])
def test_score5_reference_values(seq, expected):
    assert round(v.score5(seq), 2) == expected


@pytest.mark.parametrize("seq,expected", [
    ("ttccaaacgaacttttgtAGgga", 2.89),
    ("tgtctttttctgtgtggcAGtgg", 8.19),
    ("ttctctcttcagacttatAGcaa", -0.08),
])
def test_score3_reference_values(seq, expected):
    assert round(v.score3(seq), 2) == expected


def test_window_lengths_are_enforced():
    with pytest.raises(ValueError):
        v.score5("ACGT")
    with pytest.raises(ValueError):
        v.score3("ACGT")


# --------------------------------------------------------------------------
# Verified numbers from README.md: the challenge pack has 12 splice records,
# 7 donor / 4 acceptor / 1 both, and 6 of the 12 are indels.
# --------------------------------------------------------------------------

def test_demo_covers_the_twelve_splice_records():
    demo = json.loads((SKILL_DIR / "demo_input.txt").read_text())
    splice = [r for r in demo if r["class"] == "splice"]
    assert len(splice) == 12
    indels = [r for r in splice if len(r["ref"]) != len(r["alt"])]
    assert len(indels) == 6, "README verifies 6 of the 12 splice records are indels"


def test_demo_runs_offline_and_returns_contract_b():
    results = run_demo()
    assert len(results) == 13  # 12 splice + 1 out-of-domain record
    required = {"skill", "variant_key", "score", "direction", "confidence",
                "in_domain", "abstain_reason", "evidence"}
    for r in results:
        assert set(r) == required, f"Contract B shape drift on {r['variant_key']}"
        assert r["skill"] == "ve-splice"


def test_out_of_domain_record_abstains_with_a_reason():
    results = {r["variant_key"]: r for r in run_demo()}
    off = [r for r in results.values()
           if not r["in_domain"] and "outside that domain" in (r["abstain_reason"] or "")]
    assert off, "a protein_truncating record must abstain, not score"
    assert all(r["score"] is None for r in off)


def test_every_abstention_carries_a_reason():
    for r in run_demo():
        if not r["in_domain"]:
            assert r["abstain_reason"], f"{r['variant_key']} abstained without a reason"
            assert len(r["abstain_reason"]) > 40, "reasons must be specific, not a code"


# --------------------------------------------------------------------------
# Strand is inferred from sequence, never taken on trust. Ground truth here is
# the Ensembl GRCh37 REST gene records for these loci.
# --------------------------------------------------------------------------

ENSEMBL_STRAND = {
    "1:145606274:C:T": "-",   # POLR3C
    "1:156354347:TC:T": "+",  # RHBG
    "2:44528267:GT:G": "+",   # SLC3A1
    "4:88231392:T:TA": "-",   # HSD17B13
    "6:132203615:G:A": "+",   # ENPP1
    "11:61165731:C:CA": "+",  # TMEM216
    "11:61165741:G:C": "+",   # TMEM216
    "13:31531009:G:A": "+",   # TEX26
    "14:51378590:CT:C": "-",  # PYGL
    "17:42979026:T:C": "+",   # CCDC103
    "21:11029596:AC:A": "-",  # BAGE2
}


def test_inferred_strand_matches_ensembl():
    results = {r["variant_key"]: r for r in run_demo()}
    checked = 0
    for key, expected in ENSEMBL_STRAND.items():
        r = results[key]
        if not r["in_domain"]:
            continue
        assert r["evidence"]["inferred_strand"] == expected, (
            f"{key}: inferred {r['evidence']['inferred_strand']}, Ensembl says {expected}"
        )
        checked += 1
    assert checked >= 10


def test_canonical_dinucleotide_sits_where_the_model_expects_it():
    for r in run_demo():
        if not r["in_domain"]:
            continue
        e = r["evidence"]
        w = e["ref_window"]
        if e["site_kind"] == "donor":
            assert len(w) == v.DONOR_LEN and w[3:5] == "GT"
        else:
            assert len(w) == v.ACCEPTOR_LEN and w[18:20] == "AG"


# --------------------------------------------------------------------------
# Indel handling: the alt window must stay anchored on the same GT/AG, and
# indels in repeats must be flagged as representation-dependent.
# --------------------------------------------------------------------------

def test_indel_alt_window_stays_anchored():
    """Regression: without re-anchoring, a deletion upstream of the AG slid the
    dinucleotide out of position and produced a spurious -20 bit drop."""
    results = {r["variant_key"]: r for r in run_demo()}
    rhbg = results["1:156354347:TC:T"]
    assert rhbg["in_domain"]
    assert rhbg["evidence"]["alt_window"][18:20] == "AG"
    assert abs(rhbg["evidence"]["delta_bits"]) < 2.0, (
        "a deletion 9 nt from an intact AG must not read as a destroyed site"
    )
    assert rhbg["direction"] == "neutral"


def test_repeat_indel_is_flagged_and_downgraded():
    results = {r["variant_key"]: r for r in run_demo()}
    hsd = results["4:88231392:T:TA"]
    assert hsd["in_domain"]
    assert hsd["confidence"] == "low"
    assert any("equivalent representations" in f for f in hsd["evidence"]["flags"])


def test_equivalent_indel_span_detects_a_homopolymer():
    plus = "AAAA" + "CCCCCCCC" + "TTTT"
    # deleting one C from the run of eight can be spelled with the anchor base at any
    # of several positions; all produce the same alternate sequence
    _, _, n = v.equivalent_indel_span(plus, 3, "AC", "A")
    assert n > 1


def test_equivalent_indel_span_is_one_for_snv():
    assert v.equivalent_indel_span("ACGTACGT", 3, "T", "C")[2] == 1


# --------------------------------------------------------------------------
# Outputs and the interpretation boundary.
# --------------------------------------------------------------------------

def test_cli_demo_writes_the_full_output_tree(tmp_path):
    assert v.main(["--demo", "--output", str(tmp_path)]) == 0
    for rel in ("report.md", "result.json", "tables/results.csv",
                "reproducibility/commands.sh"):
        assert (tmp_path / rel).exists(), f"missing {rel}"
    assert v.DISCLAIMER in (tmp_path / "report.md").read_text()


BANNED = ("pathogenic", "diagnostic", "de novo", "compound heterozygous")


def test_report_makes_no_banned_claim(tmp_path):
    v.main(["--demo", "--output", str(tmp_path)])
    text = (tmp_path / "report.md").read_text().lower()
    for word in BANNED:
        assert word not in text, f"report claims {word!r}"
    blob = json.dumps(run_demo()).lower()
    for word in BANNED:
        assert word not in blob, f"result.json claims {word!r}"


def test_mismatched_reference_allele_abstains_rather_than_scoring():
    """A build mismatch must be named, not silently scored - README trap 3."""
    # real reference base at 13:31531009 is G; claiming A simulates a build mismatch
    rec = {"variant_key": "13:31531009:G:A", "chrom": "13", "pos": 31531009,
           "ref": "A", "alt": "G", "class": "splice",
           "consequence": "SPLICE_SITE_ACCEPTOR"}
    out = v.run([rec], bundled=SKILL_DIR / "demo_sequences.json")[0]
    assert out["in_domain"] is False
    assert "does not match the reference sequence" in out["abstain_reason"]
