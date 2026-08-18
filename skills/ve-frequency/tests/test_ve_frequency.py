"""Tests for ve-frequency — the frequency provenance gate.

The gate's job is to refuse. Most of these tests assert that it declines to emit
a rarity call rather than that it emits one.

Liftover truth pairs were verified against both Ensembl REST /map and gnomAD's
own liftover query before being written down here.
"""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "ve_frequency.py"
CHAIN_PATH = SKILL_DIR / "data" / "GRCh37_to_GRCh38.chain.gz"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("ve_frequency", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Trap 2: INFO/AF is GATK cohort frequency, not population frequency
# ---------------------------------------------------------------------------


def test_cohort_af_is_named_and_refused():
    """The challenge VCF's INFO/AF is cohort AF over 4 samples. Never use it."""
    module = load_module()
    info = {"AC": "2", "AF": "0.250", "AN": "8", "MLEAC": "2", "MLEAF": "0.250", "DP": "120"}
    verdict = module.assess_af_provenance(info, sample_count=4)
    assert verdict["usable_population_af"] is False
    assert verdict["af"] is None
    assert "cohort" in verdict["provenance_warning"].lower()
    assert "not population" in verdict["provenance_warning"].lower()


def test_cohort_af_warning_names_the_sample_count():
    module = load_module()
    verdict = module.assess_af_provenance({"AF": "0.25", "AN": "8"}, sample_count=4)
    assert "4 samples" in verdict["provenance_warning"]


def test_population_af_key_is_accepted():
    """A real population-frequency key is usable; AF alone is not."""
    module = load_module()
    verdict = module.assess_af_provenance(
        {"AF": "0.25", "gnomAD_AF": "0.0031"}, sample_count=4
    )
    assert verdict["usable_population_af"] is True
    assert verdict["af"] == pytest.approx(0.0031)
    assert verdict["source"] == "gnomAD_AF"


@pytest.mark.parametrize("key", ["AF_TGP", "AF_EXAC", "AF_ESP", "AF_grpmax", "gnomAD_AF"])
def test_all_known_population_af_keys_recognised(key):
    module = load_module()
    verdict = module.assess_af_provenance({key: "0.002"}, sample_count=4)
    assert verdict["usable_population_af"] is True
    assert verdict["source"] == key


def test_no_af_at_all_is_no_data_without_a_false_warning():
    module = load_module()
    verdict = module.assess_af_provenance({"DP": "120"}, sample_count=4)
    assert verdict["usable_population_af"] is False
    assert verdict["af"] is None


# ---------------------------------------------------------------------------
# Trap 3: build must be matched explicitly or declared unmatched
# ---------------------------------------------------------------------------


def test_b37_reference_detected():
    module = load_module()
    header = [
        "##fileformat=VCFv4.1",
        "##reference=file:///data/human_g1k_v37.fasta",
        "##contig=<ID=1,length=249250621>",
    ]
    assert module.detect_build(header) == "GRCh37"


def test_b38_reference_detected():
    module = load_module()
    header = ["##reference=file:///refs/GRCh38_full_analysis_set.fa",
              "##contig=<ID=1,length=248956422>"]
    assert module.detect_build(header) == "GRCh38"


def test_undeclarable_build_returns_none_rather_than_guessing():
    module = load_module()
    assert module.detect_build(["##fileformat=VCFv4.1"]) is None


def test_unknown_build_forces_no_data():
    """'Match the build explicitly or declare that you could not.'"""
    module = load_module()
    freq = module.resolve_freq(
        variant_key="1:11906068:A:G", info={"AF": "0.25"}, build=None,
        sample_count=4, mapper=None, source=None,
    )
    assert freq["class"] == "NO_DATA"
    assert freq["build"] is None
    assert "build" in freq["provenance_warning"].lower()


# ---------------------------------------------------------------------------
# Contract A conformance
# ---------------------------------------------------------------------------


REQUIRED_FREQ_KEYS = {"af", "source", "build", "class", "provenance_warning"}


def test_freq_block_matches_contract_a_shape():
    module = load_module()
    freq = module.resolve_freq(
        variant_key="1:11906068:A:G", info={"AF": "0.25"}, build="GRCh37",
        sample_count=4, mapper=None, source=None,
    )
    assert set(freq.keys()) == REQUIRED_FREQ_KEYS


@pytest.mark.parametrize("cls", ["RARE", "COMMON", "NO_DATA"])
def test_freq_class_vocabulary_is_closed(cls):
    module = load_module()
    assert cls in module.FREQ_CLASSES
    assert len(module.FREQ_CLASSES) == 3


def test_records_pass_through_untouched_apart_from_freq():
    """ve-frequency populates .freq and must not disturb upstream fields."""
    module = load_module()
    record = {
        "variant_key": "1:11906068:A:G", "chrom": "1", "pos": 11906068,
        "ref": "A", "alt": "G", "id": "rs5065",
        "genotypes": {"ISDBM322015": "0/1"},
        "segregation": {"pattern": "paternal", "carriers": ["ISDBM322015"],
                        "phased": False, "rule": "proband carries AND exactly one parent carries"},
        "raw_info": "AC=2;AF=0.250;AN=8",
    }
    out = module.annotate_records([dict(record)], build="GRCh37", sample_count=4,
                                  mapper=None, source=None)
    assert out[0]["segregation"] == record["segregation"]
    assert out[0]["variant_key"] == record["variant_key"]
    assert out[0]["genotypes"] == record["genotypes"]
    assert set(out[0]["freq"].keys()) == REQUIRED_FREQ_KEYS


# ---------------------------------------------------------------------------
# Rarity classification, when provenance permits it
# ---------------------------------------------------------------------------


def test_population_af_below_threshold_is_rare():
    module = load_module()
    freq = module.resolve_freq(
        variant_key="1:1:A:G", info={"gnomAD_AF": "0.0004"}, build="GRCh37",
        sample_count=4, mapper=None, source=None,
    )
    assert freq["class"] == "RARE"
    assert freq["af"] == pytest.approx(0.0004)
    assert freq["source"] == "gnomAD_AF"


def test_population_af_above_threshold_is_common():
    module = load_module()
    freq = module.resolve_freq(
        variant_key="1:1:A:G", info={"gnomAD_AF": "0.12"}, build="GRCh37",
        sample_count=4, mapper=None, source=None,
    )
    assert freq["class"] == "COMMON"


def test_absent_at_well_covered_site_is_rare_with_a_bound():
    """Absence only counts when the site was actually looked at."""
    module = load_module()
    record = {"ac": 0, "an": 303936, "coverage_over_20": 0.96, "filters": ["AC0"]}
    assert module.classify_lookup(record)["class"] == "RARE"
    assert module.classify_lookup(record)["upper_bound_af"] == pytest.approx(3 / 303936)


def test_absent_at_uncovered_site_is_no_data_not_rare():
    module = load_module()
    record = {"ac": 0, "an": 1200, "coverage_over_20": 0.32, "filters": []}
    assert module.classify_lookup(record)["class"] == "NO_DATA"


def test_missing_lookup_record_is_no_data():
    module = load_module()
    assert module.classify_lookup(None)["class"] == "NO_DATA"


# ---------------------------------------------------------------------------
# Liftover: b37 data, gnomAD v4 is GRCh38-only
# ---------------------------------------------------------------------------

LIFT_TRUTH = [("17", 41246481, 43094464), ("1", 55505647, 55039974),
              ("13", 32906729, 32332592), ("7", 117199644, 117559590)]


@pytest.mark.skipif(not CHAIN_PATH.exists(), reason="chain not vendored")
@pytest.mark.parametrize("chrom,pos37,pos38", LIFT_TRUTH)
def test_liftover_matches_ground_truth(chrom, pos37, pos38):
    module = load_module()
    mapper = module.ChainMapper.from_file(CHAIN_PATH)
    mapped = mapper.map(chrom, pos37)
    assert mapped is not None and mapped.pos == pos38 and mapped.strand == "+"


def test_minus_strand_block_reverse_complements(tmp_path):
    module = load_module()
    chain = tmp_path / "s.chain"
    chain.write_text("chain 1000 A 300 + 0 150 A 350 + 50 200 1\n150\n\n"
                     "chain 1000 A 300 + 150 300 A 350 - 150 300 2\n150\n\n", encoding="utf-8")
    mapper = module.ChainMapper.from_file(chain)
    flipped = module.lift_variant("A", 200, "G", "A", mapper)
    assert flipped.status == "lifted"
    assert (flipped.ref, flipped.alt) == ("C", "T")


def test_failed_lift_is_no_data_never_rare(tmp_path):
    module = load_module()
    chain = tmp_path / "s.chain"
    chain.write_text("chain 1000 A 300 + 0 150 A 350 + 50 200 1\n150\n\n", encoding="utf-8")
    mapper = module.ChainMapper.from_file(chain)
    lift = module.lift_variant("A", 275, "G", "A", mapper)
    assert lift.status == "unresolvable"


# ---------------------------------------------------------------------------
# CLI / demo
# ---------------------------------------------------------------------------


def test_demo_runs_offline(tmp_path):
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    assert (out / "report.md").exists()
    assert (out / "result.json").exists()


def test_demo_output_is_contract_a(tmp_path):
    out = tmp_path / "out2"
    subprocess.run([sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
                   capture_output=True, text=True, timeout=300, check=True)
    payload = json.loads((out / "result.json").read_text(encoding="utf-8"))
    records = payload["records"]
    assert records, "demo must emit records"
    for record in records:
        assert set(record["freq"].keys()) == REQUIRED_FREQ_KEYS
        assert record["freq"]["class"] in load_module().FREQ_CLASSES
    assert payload["disclaimer"] == DISCLAIMER


def test_demo_refuses_rarity_on_cohort_af_only_input(tmp_path):
    """The headline behaviour: b37 VCF with only GATK AF yields no rarity calls."""
    out = tmp_path / "out3"
    subprocess.run([sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
                   capture_output=True, text=True, timeout=300, check=True)
    payload = json.loads((out / "result.json").read_text(encoding="utf-8"))
    cohort_only = [r for r in payload["records"]
                   if r["freq"]["provenance_warning"]
                   and "cohort" in r["freq"]["provenance_warning"].lower()]
    assert cohort_only, "demo must demonstrate the cohort-AF trap"
    assert all(r["freq"]["class"] == "NO_DATA" for r in cohort_only)
    assert all(r["freq"]["af"] is None for r in cohort_only)


def test_report_states_the_build_and_the_refusal_count(tmp_path):
    out = tmp_path / "out4"
    subprocess.run([sys.executable, str(MODULE_PATH), "--demo", "--output", str(out)],
                   capture_output=True, text=True, timeout=300, check=True)
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "GRCh37" in text
    assert "NO_DATA" in text
    assert DISCLAIMER in text
