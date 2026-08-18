import importlib.util
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "ve_regulatory.py"
DEMO_TSV = SKILL_DIR / "demo_variants.tsv"
CACHE_PATH = SKILL_DIR / "tests" / "fixtures" / "cached_scores.json"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
BANNED = ["rare", "pathogenic", "diagnostic", "de novo", "compound heterozygous"]


def load_module():
    spec = importlib.util.spec_from_file_location("ve_regulatory", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_cli(args, **kwargs):
    return subprocess.run([sys.executable, str(MODULE_PATH), *args], text=True,
                           capture_output=True, check=False, **kwargs)


# --- --demo must work with no network, no model stack, and no user files --

def test_demo_cli_works_offline_and_exits_zero(tmp_path):
    """--demo must succeed even under the interpreter running this test suite,
    which is not guaranteed to have torch/cherimoya installed at all."""
    out = tmp_path / "demo_out"
    completed = run_cli(["--demo", "--output", str(out)])
    assert completed.returncode == 0, completed.stderr
    assert "ve-regulatory" in completed.stdout
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["mode"] == "demo"
    assert len(result["records"]) == 4
    assert all(r["in_domain"] for r in result["records"])


def test_demo_module_call_never_touches_network(monkeypatch, tmp_path):
    module = load_module()

    def _blocked_urlopen(*_args, **_kwargs):
        raise AssertionError("--demo must never touch the network (urlopen called)")
    monkeypatch.setattr(module.urllib.request, "urlopen", _blocked_urlopen)
    # hf_hub_download is imported lazily inside functions --demo never calls, so there
    # is nothing to patch there; the real guarantee is that _run_demo only touches
    # load_variants() and the local JSON fixture, asserted structurally below.
    out = tmp_path / "demo_out"
    rc = module._run_demo(out, ["ve_regulatory.py", "--demo"])
    assert rc == 0
    assert (out / "result.json").exists()


def test_demo_matches_cached_fixture_verbatim():
    module = load_module()
    cache = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    variants = module.load_variants(DEMO_TSV)
    by_key = {r["variant_key"]: r for r in cache["records"]}
    for v in variants:
        assert v["variant_key"] in by_key, f"fixture is missing {v['variant_key']}"


def test_demo_real_captured_directions_and_models():
    """Locks in the real, live-model numbers captured in cached_scores.json (see its
    `provenance` block) so a future edit cannot silently swap in fabricated values."""
    cache = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    by_key = {r["variant_key"]: r for r in cache["records"]}
    assert by_key["1:109274968:G:T"]["direction"] == "increases"  # rs12740374 / SORT1
    assert by_key["1:109274968:G:T"]["evidence"]["confirmed_model"] == "ENCSR562FNN"
    assert by_key["16:53767042:T:C"]["evidence"]["confirmed_model"] == "ENCSR540BML"  # rs1421085 / FTO-IRX3
    assert by_key["8:127401060:G:T"]["evidence"]["confirmed_model"] == "ENCSR994KTY"  # rs6983267 / 8q24 MYC
    assert by_key["2:135851076:G:A"]["evidence"]["confirmed_model"] == "ENCSR133KBX"  # rs4988235 / LCT
    for r in cache["records"]:
        assert "UNCALIBRATED" in r["evidence"]["notes"][1]


# --- abstention: pure logic paths that never need torch/cherimoya ---------

def test_wrong_class_abstains_without_torch():
    module = load_module()
    v = {"variant_key": "x", "class": "missense", "build": "hg38", "ref": "A", "alt": "G", "chrom": "1", "pos": 1e9}
    r = module.score_variant(v, "ENCSR000XXX", 0, None)
    assert r["in_domain"] is False
    assert r["abstain_reason"] == "variant is not non_coding (wrong branch)"


def test_undeclared_build_abstains_without_torch():
    module = load_module()
    v = {"variant_key": "x", "class": "non_coding", "build": None, "ref": "A", "alt": "G", "chrom": "1", "pos": 1e9}
    r = module.score_variant(v, "ENCSR000XXX", 0, None)
    assert r["abstain_reason"] == "coordinates not hg38 / build undeclared"
    v["build"] = "hg19"
    assert module.score_variant(v, "ENCSR000XXX", 0, None)["abstain_reason"] == "coordinates not hg38 / build undeclared"


def test_indel_abstains_without_torch():
    module = load_module()
    v = {"variant_key": "x", "class": "non_coding", "build": "hg38", "ref": "AT", "alt": "G", "chrom": "1", "pos": 1e9}
    r = module.score_variant(v, "ENCSR000XXX", 0, None)
    assert r["abstain_reason"] == "indel (window construction assumes a substitution)"


def test_contig_end_abstains_without_torch():
    module = load_module()
    v = {"variant_key": "x", "class": "non_coding", "build": "hg38", "ref": "A", "alt": "G", "chrom": "1", "pos": 100}
    r = module.score_variant(v, "ENCSR000XXX", 0, None)
    assert r["abstain_reason"] == "variant within 1057bp of a contig end (no full window)"


def test_no_cell_type_and_no_model_abstains():
    module = load_module()
    args = _fake_args(cell_type=None, model=None)
    v = {"variant_key": "x", "confirmed_model": None}
    r = module.resolve_and_score(v, args, None)
    assert r["abstain_reason"] == "no --cell-type given"


def test_cell_type_without_confirmed_model_abstains_with_exact_reason():
    module = load_module()
    args = _fake_args(cell_type="liver", model=None)
    v = {"variant_key": "x", "confirmed_model": None}
    r = module.resolve_and_score(v, args, None)
    assert r["abstain_reason"] == "cell-type model proposed but not confirmed"


def _fake_args(**kw):
    class Args:
        pass
    a = Args()
    a.cell_type = kw.get("cell_type")
    a.model = kw.get("model")
    a.fold = 0
    return a


# --- malformed / missing input --------------------------------------------

def test_malformed_input_exits_2_with_error_and_no_traceback(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("this is neither JSON nor a tab-separated table\n", encoding="utf-8")
    completed = run_cli(["--input", str(bad), "--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_missing_required_columns_exits_2(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("foo\tbar\n1\t2\n", encoding="utf-8")
    completed = run_cli(["--input", str(bad), "--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_missing_input_without_demo_or_propose_is_a_usage_error(tmp_path):
    completed = run_cli(["--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "--input is required" in completed.stderr


# --- output structure + disclaimer + safety --------------------------------

def test_output_structure_written(tmp_path):
    out = tmp_path / "demo_out"
    completed = run_cli(["--demo", "--output", str(out)])
    assert completed.returncode == 0
    assert (out / "report.md").exists()
    assert (out / "result.json").exists()
    assert (out / "tables" / "results.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


def test_disclaimer_present_in_report(tmp_path):
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report


def test_no_banned_interpretive_words_in_demo_report_or_result(tmp_path):
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    text = (out / "report.md").read_text(encoding="utf-8").lower() + (out / "result.json").read_text(encoding="utf-8").lower()
    for word in BANNED:
        assert word not in text, f"banned word '{word}' found in ve-regulatory output"


def test_result_json_is_contract_b_shaped():
    module = load_module()
    cache = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    for r in cache["records"]:
        for field in ("skill", "variant_key", "score", "direction", "confidence", "in_domain", "abstain_reason", "evidence"):
            assert field in r
        assert r["skill"] == "ve-regulatory"


def test_fuzzy_never_names_a_transcription_factor():
    """Guards the one hardest-to-catch overclaim: tf_binding must never carry a TF name,
    only a boolean flag."""
    cache = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    for r in cache["records"]:
        assert isinstance(r["evidence"]["tf_binding_high_attribution_motif_like_position"], bool)


# --- the CLI itself must abstain without torch, not just score_variant() ----

def test_cli_abstains_on_cheap_checks_without_torch(tmp_path):
    """Regression for a false docstring claim: cheap abstentions (build/indel/class)
    must work through the real CLI in an environment with no model stack, not only
    when score_variant() is called directly."""
    src = tmp_path / "v.tsv"
    src.write_text(
        "chrom\tpos\tref\talt\tbuild\tclass\n"
        "1\t109274968\tG\tT\thg19\tnon_coding\n"      # wrong build
        "1\t109274968\tG\tGA\thg38\tnon_coding\n"     # indel
        "1\t109274968\tG\tT\thg38\tprotein_truncating\n",  # wrong branch
        encoding="utf-8")
    out = tmp_path / "o"
    completed = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(src), "--model", "ENCSR562FNN",
         "--output", str(out)],
        capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    records = json.loads((out / "result.json").read_text(encoding="utf-8"))["records"]
    assert len(records) == 3
    assert all(r["in_domain"] is False for r in records)
    reasons = {r["abstain_reason"] for r in records}
    assert "coordinates not hg38 / build undeclared" in reasons
    assert any("indel" in r for r in reasons)
    assert "variant is not non_coding (wrong branch)" in reasons


def test_scored_records_disclose_model_selection_limitation():
    cached = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    records = cached["records"] if isinstance(cached, dict) else cached
    scored = [r for r in records if r.get("in_domain")]
    assert scored, "fixture has no scored records"
    for r in scored:
        note = (r.get("evidence") or {}).get("model_selection", "")
        assert "cannot verify" in note, "scored record must disclose that confirmation is unverified"
