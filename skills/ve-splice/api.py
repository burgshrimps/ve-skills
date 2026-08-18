"""Importable entrypoint for ve-splice.

    from api import run
    results = run(contract_a_records, reference="human_g1k_v37.fasta")

Returns a list of Contract B records, one per input record.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from ve_splice import HERE, evaluate_record, run as _run

__all__ = ["run", "run_demo", "evaluate_record"]


def run(records: Iterable[dict], reference: str | Path | None = None,
        bundled: str | Path | None = None) -> list[dict]:
    """Score Contract A records, returning Contract B results.

    Exactly one of `reference` (an indexed FASTA matching the input build) or
    `bundled` (pre-extracted windows) must be given.
    """
    if bool(reference) == bool(bundled):
        raise ValueError("pass exactly one of reference= or bundled=")
    return _run(records, fasta=Path(reference) if reference else None,
                bundled=Path(bundled) if bundled else None)


def run_demo() -> list[dict]:
    """Score the bundled demo records with no network and no reference FASTA."""
    import json

    records = json.loads((HERE / "demo_input.txt").read_text())
    return _run(records, bundled=HERE / "demo_sequences.json")
