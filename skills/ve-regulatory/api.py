"""Importable API entry point for ve-regulatory.

Mirrors the CLI: scoring requires an explicitly supplied CATv1 accession. This
module will not choose a cell-type model for you -- call propose() first, review
the shortlist, then pass the accession you picked.
"""
from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_MOD_PATH = Path(__file__).resolve().parent / "ve_regulatory.py"
_spec = spec_from_file_location("ve_regulatory", _MOD_PATH)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)

VERSION = "0.1.0"


def propose(cell_type: str, top_n: int = 10) -> list[dict]:
    """Rank CATv1 experiments against a biosample term. Chooses nothing.

    Returns a ranked shortlist of dicts (accession, biosample, assay, name_match,
    count_pearson_fold0) for a human to review.
    """
    return _mod.propose(cell_type, top_n)


def list_cell_types() -> list[dict]:
    """Every GRCh38 human biosample with a CATv1 model, with experiment counts."""
    return _mod.list_cell_types()


def run(input_path: str, model: str, output_dir: str = "/tmp/ve-regulatory") -> dict:
    """Score variants against one caller-supplied CATv1 accession.

    Args:
        input_path: Contract A JSON array, or a TSV of variants.
        model: CATv1 experiment_accession, chosen from a propose() shortlist.
            Required -- this function will not select one for you.
        output_dir: where report.md / result.json / tables / reproducibility go.

    Returns the result dict, whose `records` are Contract B objects. Records that
    fall outside the skill's domain carry in_domain=False and an abstain_reason
    rather than a score.
    """
    if not model:
        raise ValueError("model is required: call propose() and pass a reviewed accession")
    args = _mod.argparse.Namespace(
        input=Path(input_path), output=Path(output_dir), model=model,
        cell_type=None, fold=0, demo=False, propose=False, list_cell_types=False, top=10)
    variants = _mod.load_variants(Path(input_path))
    try:
        performance = _mod.load_performance()
    except Exception:
        performance = None
    records = [_mod.resolve_and_score(v, args, performance) for v in variants]
    result = {"skill": "ve-regulatory", "version": VERSION, "mode": "score",
              "records": records, "disclaimer": _mod.DISCLAIMER}
    _mod.write_outputs(
        Path(output_dir), result,
        lambda p: _mod._write_score_report(p, records, input_path, "Called via api.run()"),
        ["python", str(_MOD_PATH), "--input", input_path, "--model", model])
    return result
