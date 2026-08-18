import importlib.util
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "ve_router.py"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"
CHALLENGE_VCF = Path("/Users/burgshrimps/project/personal/clawbio-hack/challenge1-b37-segregation.vcf.gz")
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("ve_router", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_cli(args, **kwargs):
    return subprocess.run([sys.executable, str(MODULE_PATH), *args], text=True,
                           capture_output=True, check=False, **kwargs)


# --- Acceptance test: this is the definition of done ---------------------

def test_acceptance_challenge_vcf_routing_counts():
    assert CHALLENGE_VCF.exists(), f"challenge VCF not found at {CHALLENGE_VCF}"
    out = Path("/tmp/ve_router_test_acceptance")
    completed = run_cli(["--input", str(CHALLENGE_VCF), "--output", str(out)])
    assert completed.returncode == 0, completed.stderr
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    summary = result["summary"]
    assert summary["total"] == 68
    by_class = summary["by_class"]
    assert by_class.get("protein_truncating", 0) == 56
    assert by_class.get("splice", 0) == 12
    assert by_class.get("missense", 0) == 0
    assert by_class.get("non_coding", 0) == 0
    assert by_class.get("unroutable", 0) == 0


def test_functional_class_trap_never_drives_routing():
    """A real record: STOP_LOST(HIGH|MISSENSE|...). Must route on the Effect
    name (STOP_LOST -> protein_truncating), never on FunctionalClass=MISSENSE."""
    module = load_module()
    info = {"EFF": "STOP_LOST(HIGH|MISSENSE|Tga/Cga|*152R|151|NPPA||CODING|NM_006172.3|3|1)"}
    variant = {"variant_key": "1:11906068:A:G", "chrom": "1", "pos": 11906068,
               "ref": "A", "alt": "G", "id": "rs5065", "genotypes": {}}
    record = module.route_one(variant, info, "GRCh37", allow_network=False)
    assert record["class"] == "protein_truncating"
    assert record["consequence"] == "STOP_LOST"
    assert record["routing"]["selected"]["functional_class"] == "MISSENSE"


# --- --demo must work with no network and no user files -------------------

def test_demo_cli_works_offline_and_exits_zero(monkeypatch):
    def _blocked(*_args, **_kwargs):
        raise AssertionError("--demo must never touch the network")
    monkeypatch.setattr(urllib.request, "urlopen", _blocked)
    out = Path("/tmp/ve_router_test_demo")
    completed = run_cli(["--demo", "--output", str(out)])
    assert completed.returncode == 0, completed.stderr
    assert "ve-router" in completed.stdout
    result = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert result["summary"]["total"] == 10
    assert set(result["summary"]["by_class"]) == {
        "protein_truncating", "splice", "missense", "non_coding", "unroutable",
    }


def test_demo_module_call_never_touches_network(monkeypatch):
    module = load_module()

    def _blocked(*_args, **_kwargs):
        raise AssertionError("network touched")
    monkeypatch.setattr(module.urllib.request, "urlopen", _blocked)
    records = module.load_variants(DEMO_INPUT, None)
    routed = [module.route_one(v, i, b, allow_network=False) for v, i, b in records]
    unroutable = [r for r in routed if r["class"] == "unroutable"]
    assert len(unroutable) == 1
    assert unroutable[0]["routing"]["unroutable_reason"] == (
        "no consequence annotation present and no annotator reachable"
    )


# --- malformed input --------------------------------------------------

def test_malformed_input_exits_2_with_error_and_no_traceback(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("this is not a vcf and not json either\n", encoding="utf-8")
    completed = run_cli(["--input", str(bad), "--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_malformed_vcf_data_line_exits_2(tmp_path):
    bad = tmp_path / "bad.vcf"
    bad.write_text("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "1\t100\n", encoding="utf-8")
    completed = run_cli(["--input", str(bad), "--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_missing_input_file_exits_2(tmp_path):
    completed = run_cli(["--input", str(tmp_path / "nope.vcf"), "--output", str(tmp_path / "out")])
    assert completed.returncode == 2
    assert "ERROR:" in completed.stderr
    assert "Traceback" not in completed.stderr


# --- disclaimer -------------------------------------------------------

def test_disclaimer_present_in_report():
    out = Path("/tmp/ve_router_test_disclaimer")
    completed = run_cli(["--demo", "--output", str(out)])
    assert completed.returncode == 0
    report = (out / "report.md").read_text(encoding="utf-8")
    assert DISCLAIMER in report


# --- discarded alternatives are emitted, never silently collapsed ------

def test_discarded_alternatives_emitted_for_multi_annotation_record():
    module = load_module()
    records = module.load_variants(DEMO_INPUT, None)
    routed = {v["variant_key"]: module.route_one(v, i, b, allow_network=False) for v, i, b in records}
    multi_gene_record = routed["1:11906068:A:G"]
    assert multi_gene_record["class"] == "protein_truncating"
    discarded = multi_gene_record["routing"]["discarded"]
    assert len(discarded) == 4
    discarded_genes = {d["gene"] for d in discarded}
    assert discarded_genes == {"CLCN6", "NPPA-AS1"}
    assert multi_gene_record["routing"]["selection_rule"] is not None


def test_output_structure_written():
    out = Path("/tmp/ve_router_test_structure")
    completed = run_cli(["--demo", "--output", str(out)])
    assert completed.returncode == 0
    assert (out / "report.md").exists()
    assert (out / "result.json").exists()
    assert (out / "tables" / "results.csv").exists()
    assert (out / "reproducibility" / "commands.sh").exists()


# --- routing classes covered by the demo file --------------------------

def test_demo_covers_all_five_classes_and_annotation_formats():
    module = load_module()
    records = module.load_variants(DEMO_INPUT, None)
    routed = [module.route_one(v, i, b, allow_network=False) for v, i, b in records]
    by_key = {r["variant_key"]: r for r in routed}
    assert by_key["2:1000000:C:T"]["class"] == "missense"
    assert by_key["2:2000000:G:A"]["class"] == "non_coding"
    assert by_key["2:3000000:A:T"]["class"] == "unroutable"
    ann_record = by_key["3:500000:G:A"]
    assert ann_record["class"] == "protein_truncating"
    assert ann_record["routing"]["annotation_source"] == "ANN"
    both_splice = by_key["1:156354347:TC:T"]
    assert both_splice["class"] == "splice"
    assert both_splice["consequence"] == "SPLICE_SITE_ACCEPTOR"
