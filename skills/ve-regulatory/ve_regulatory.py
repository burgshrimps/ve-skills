#!/usr/bin/env python3
"""ve-regulatory: CATv1 (Cherimoya) ref-vs-alt accessibility delta for non-coding variants.

Human-in-the-loop cell-type selection. --propose ranks ENCODE DNase/ATAC experiments
against a biosample term and CHOOSES NOTHING. Scoring requires an accession supplied
explicitly, via --model (applies to every row) or a per-row `confirmed_model` column.

The review step is a workflow convention, not a technical guarantee: this tool cannot
tell whether a caller-supplied accession came from a reviewed --propose shortlist or
was invented. It refuses to pick one itself, and records that limitation in
evidence.model_selection on every scored record. Everything else is Contract A in,
Contract B out.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import math
import sys
import urllib.error
import urllib.request
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DEMO_TSV = SKILL_DIR / "demo_variants.tsv"
DEMO_CACHE = SKILL_DIR / "tests" / "fixtures" / "cached_scores.json"
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
HF_REPO = "programmable-genomics/CATv1"
WINDOW = 2114
FLANK_L, FLANK_R = 1057, 1056  # pos-FLANK_L .. pos+FLANK_R inclusive == 2114bp, variant at offset FLANK_L
ALPHABET = "ACGT"
ENSEMBL_URL = "https://rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}?content-type=text/plain"
ATTR_HALF_WINDOW = 25       # +/-25bp local attribution window: CPU cost tradeoff, see SKILL.md
ATTR_BATCH_SIZE = 8         # fastest measured per-sequence CPU throughput for this checkpoint size
DELTA_NOISE_FLOOR = 0.1     # natural-log units; direction threshold, NOT a significance test
FUZZY_MIN = 0.35
CONFIDENCE_BINS = [(0.85, "high"), (0.6, "medium"), (0.0, "low")]
INSTALL_HELP = (
    "torch/cherimoya not importable. Verified recipe, Python 3.13 CPU:\n"
    "  python3 -m venv venv-cat && venv-cat/bin/pip install --upgrade pip\n"
    "  venv-cat/bin/pip install torch\n"
    "  venv-cat/bin/pip install --no-deps cherimoya bpnet-lite\n"
    "  venv-cat/bin/pip install tangermeme numpy scipy pandas h5py tqdm huggingface_hub\n"
    "bpnet-lite is an unguarded import inside cherimoya/losses.py; without it, "
    "`import cherimoya` fails even though pip never lists it as a dependency."
)

# --- metadata / performance / propose (network; never touched by --demo) ---

def _load_tsv(filename: str) -> list[dict]:
    from huggingface_hub import hf_hub_download
    path = hf_hub_download(HF_REPO, filename)
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_metadata() -> list[dict]:
    return _load_tsv("CATv1-metadata.tsv")


def load_performance() -> dict[str, float]:
    rows = _load_tsv("performance.tsv")
    return {r["experiment_accession"]: float(r["count_pearson"])
            for r in rows if r["fold"] == "0" and r["count_pearson"]}


def _fuzzy(term: str, name: str) -> float:
    """Token-level difflib similarity. Deliberately NOT a synonym map -- see SKILL.md
    for verified cases where the top-ranked row by this score is the wrong tissue."""
    term_l, name_l = term.lower().strip(), name.lower().strip()
    if term_l == name_l:
        return 1.0
    if term_l in name_l:
        return 0.85 + 0.15 * (len(term_l) / len(name_l))
    best = 0.0
    for tt in term_l.split():
        for tok in name_l.replace("-", " ").split():
            best = max(best, difflib.SequenceMatcher(None, tt, tok).ratio())
    return best * 0.8


def propose(cell_type: str, top_n: int = 10) -> list[dict]:
    meta, perf = load_metadata(), load_performance()
    rows = []
    for r in meta:
        if r["assembly"] != "GRCh38" or r["organism"] != "Homo sapiens":
            continue
        acc = r["experiment_accession"]
        score = _fuzzy(cell_type, r["experiment_biosample_term_name"])
        if score < FUZZY_MIN:
            continue
        rows.append({"accession": acc, "biosample": r["experiment_biosample_term_name"],
                     "assay": r["assay_term_name"], "name_match": round(score, 3),
                     "count_pearson_fold0": perf.get(acc)})
    rows.sort(key=lambda r: (-r["name_match"], -(r["count_pearson_fold0"] or 0)))
    return rows[:top_n]


def list_cell_types() -> list[dict]:
    meta = load_metadata()
    seen: dict[str, dict] = {}
    for r in meta:
        if r["assembly"] != "GRCh38" or r["organism"] != "Homo sapiens":
            continue
        entry = seen.setdefault(r["experiment_biosample_term_name"],
                                 {"biosample": r["experiment_biosample_term_name"], "n_experiments": 0, "assays": set()})
        entry["n_experiments"] += 1
        entry["assays"].add(r["assay_term_name"])
    out = [{"biosample": v["biosample"], "n_experiments": v["n_experiments"], "assays": sorted(v["assays"])} for v in seen.values()]
    out.sort(key=lambda r: (-r["n_experiments"], r["biosample"]))
    return out

# --- variant loading ---------------------------------------------------

def load_variants(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8")
    stripped = text.lstrip()
    if stripped.startswith("[") or stripped.startswith("{"):
        data = json.loads(text)
        data = data.get("records", [data]) if isinstance(data, dict) else data
        if not isinstance(data, list) or not data:
            raise ValueError("JSON input must be a non-empty array of Contract A records")
        rows = data
    elif "\t" in stripped.splitlines()[0]:
        rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        if not rows:
            raise ValueError("TSV input contains no data rows")
    else:
        raise ValueError("Input is not a recognised Contract A JSON array or tab-separated variant table")
    variants = []
    for row in rows:
        missing = {"chrom", "pos", "ref", "alt"} - row.keys()
        if missing:
            raise ValueError(f"Variant record missing required fields: {sorted(missing)}")
        routing_build = row.get("routing", {}).get("build") if isinstance(row.get("routing"), dict) else None
        variants.append({
            "chrom": str(row["chrom"]), "pos": int(row["pos"]), "ref": str(row["ref"]), "alt": str(row["alt"]),
            "id": row.get("id"), "gene": row.get("gene"), "class": row.get("class") or "non_coding",
            "build": row.get("build") or routing_build, "confirmed_model": row.get("confirmed_model") or None,
            "variant_key": row.get("variant_key") or f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}",
        })
    return variants

# --- sequence + model ---------------------------------------------------

def fetch_window(chrom: str, pos: int) -> str:
    """One retry on timeout: the Ensembl REST API has been observed to hang on a first
    request and succeed immediately on a second (see SKILL.md Gotchas)."""
    start, end = pos - FLANK_L, pos + FLANK_R
    if start < 1:
        raise ValueError("variant within 1057bp of a contig end (no full window)")
    url = ENSEMBL_URL.format(chrom=chrom, start=start, end=end)
    req = urllib.request.Request(url, headers={"Accept": "text/plain"})
    last_exc: Exception | None = None
    for _attempt in range(2):
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                seq = resp.read().decode("utf-8").strip()
            if len(seq) != WINDOW:
                raise ValueError(f"Ensembl fetch failed or returned length != 2114: got {len(seq)} bp")
            return seq
        except (urllib.error.URLError, OSError, TimeoutError, ValueError) as exc:
            last_exc = exc
    raise ValueError(f"Ensembl fetch failed or returned length != 2114: {last_exc}")


def require_model_stack() -> None:
    try:
        import torch  # noqa: F401
        import cherimoya  # noqa: F401
        import tangermeme.saturation_mutagenesis  # noqa: F401
        import huggingface_hub  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(f"{INSTALL_HELP}\n(missing: {exc.name})") from exc


_MODEL_CACHE: dict[str, object] = {}


def load_cherimoya(accession: str, fold: int = 0):
    if accession in _MODEL_CACHE:
        return _MODEL_CACHE[accession]
    from huggingface_hub import hf_hub_download
    from cherimoya import Cherimoya
    try:
        path = hf_hub_download(HF_REPO, f"models/{accession}/cherimoya.fold_{fold}.torch")
    except Exception as exc:
        raise LookupError(f"no CATv1 checkpoint for accession {accession} fold {fold} ({exc})") from exc
    _MODEL_CACHE[accession] = Cherimoya.load(path, device="cpu", compile=False).eval()
    return _MODEL_CACHE[accession]


def one_hot(seq: str):
    import torch
    x = torch.zeros(4, len(seq), dtype=torch.float32)
    for i, c in enumerate(seq.upper()):
        if c in ALPHABET:
            x[ALPHABET.index(c), i] = 1.0
    return x


def _abstain(variant_key: str, reason: str) -> dict:
    return {"skill": "ve-regulatory", "variant_key": variant_key, "score": None, "direction": None,
            "confidence": None, "in_domain": False, "abstain_reason": reason, "evidence": {"flags": []}}


def score_variant(variant: dict, accession: str, fold: int, performance: dict[str, float] | None) -> dict:
    """Score one variant, or abstain with a reason.

    The cheap checks (class/build/indel/window/ref-match) run before the model
    stack is required, so a variant that abstains on those grounds does so even
    when torch/cherimoya are absent. require_model_stack() is called lazily, only
    once a variant has survived the cheap checks and genuinely needs inference.
    """
    key = variant["variant_key"]
    if variant["class"] != "non_coding":
        return _abstain(key, "variant is not non_coding (wrong branch)")
    if (variant["build"] or "").strip().lower() not in {"hg38", "grch38"}:
        return _abstain(key, "coordinates not hg38 / build undeclared")
    ref, alt = variant["ref"], variant["alt"]
    if len(ref) != 1 or len(alt) != 1 or ref.upper() not in ALPHABET or alt.upper() not in ALPHABET:
        return _abstain(key, "indel (window construction assumes a substitution)")
    try:
        seq = fetch_window(variant["chrom"], variant["pos"])
    except ValueError as exc:
        return _abstain(key, str(exc))
    var_idx = FLANK_L
    if seq[var_idx].upper() != ref.upper():
        return _abstain(key, f"reference allele mismatch: declared ref {ref} but hg38 has {seq[var_idx]} at this position")
    require_model_stack()
    try:
        model = load_cherimoya(accession, fold)
    except LookupError as exc:
        return _abstain(key, f"requested cell type has no CATv1 model: {exc}")

    import torch
    from cherimoya import ControlWrapper, LogCountWrapper
    from tangermeme.saturation_mutagenesis import saturation_mutagenesis

    alt_seq = seq[:var_idx] + alt.upper() + seq[var_idx + 1:]
    X = torch.stack([one_hot(seq), one_hot(alt_seq)])
    with torch.no_grad():
        _profile, y_counts = model(X)
    ref_lc, alt_lc = y_counts[0, 0].item(), y_counts[1, 0].item()
    delta = alt_lc - ref_lc
    direction = "increases" if delta > DELTA_NOISE_FLOOR else "decreases" if delta < -DELTA_NOISE_FLOOR else "no meaningful change"
    ref_counts, alt_counts = math.expm1(ref_lc), math.expm1(alt_lc)
    fold_change = alt_counts / ref_counts if ref_counts > 1e-9 else None

    wrapped = LogCountWrapper(ControlWrapper(model))
    lo, hi = var_idx - ATTR_HALF_WINDOW, var_idx + ATTR_HALF_WINDOW + 1
    attr = saturation_mutagenesis(wrapped, X, start=lo, end=hi, batch_size=ATTR_BATCH_SIZE, verbose=False)
    importance = (attr * X[:, :, lo:hi]).sum(dim=1)
    local_idx = var_idx - lo
    ref_attr, alt_attr = importance[0, local_idx].item(), importance[1, local_idx].item()
    window_max_abs = importance[0].abs().max().item()
    motif_like = window_max_abs > 0 and abs(ref_attr) >= 0.5 * window_max_abs

    pearson = (performance or {}).get(accession)
    confidence = next(label for threshold, label in CONFIDENCE_BINS if (pearson or 0) >= threshold)
    evidence = {
        "confirmed_model": accession, "fold": fold, "count_pearson_fold0": pearson,
        "model_selection": ("accession supplied by the caller; ve-regulatory cannot verify that a "
                            "--propose shortlist was produced or reviewed by a human"),
        "window": f"{variant['chrom']}:{variant['pos'] - FLANK_L}-{variant['pos'] + FLANK_R} (hg38)",
        "predicted_accessibility": {"ref_log_count": round(ref_lc, 4), "alt_log_count": round(alt_lc, 4)},
        "magnitude_uncalibrated": {"delta_log_count": round(delta, 4),
                                    "fold_change_approx": round(fold_change, 3) if fold_change is not None else None},
        "contribution_delta": {"ref_attribution_at_variant": round(ref_attr, 5), "alt_attribution_at_variant": round(alt_attr, 5),
                                "delta": round(alt_attr - ref_attr, 5), "window_bp": 2 * ATTR_HALF_WINDOW + 1},
        "tf_binding_high_attribution_motif_like_position": motif_like,
        "notes": [
            "predicted_accessibility is CATv1's own model output: a continuous score, not a called ATAC/DNase peak; ENCODE peak BEDs were not consulted.",
            "magnitude_uncalibrated is UNCALIBRATED: a published evaluation of this model class found predicted fold-changes fall roughly an order of magnitude below measured ones (rs12740374/SORT1: ~12-fold measured). Treat only direction as informative.",
            "contribution_delta is native tangermeme.saturation_mutagenesis attribution on the log-count head -- the one number here with no interpretive layer added.",
            "tf_binding flags a high-attribution, motif-like position only. CATv1 gives no TF identity; no transcription factor is ever named.",
        ],
        "flags": [],
    }
    return {"skill": "ve-regulatory", "variant_key": key, "score": round(abs(delta), 4), "direction": direction,
            "confidence": confidence, "in_domain": True, "abstain_reason": None, "evidence": evidence}


def resolve_and_score(variant: dict, args: argparse.Namespace, performance: dict[str, float] | None) -> dict:
    accession = variant.get("confirmed_model") or args.model
    if not accession:
        reason = "cell-type model proposed but not confirmed" if args.cell_type else "no --cell-type given"
        return _abstain(variant["variant_key"], reason)
    return score_variant(variant, accession, args.fold, performance)

# --- output writers ---------------------------------------------------

_CSV_FIELDS = ["variant_key", "chrom", "pos", "ref", "alt", "id", "confirmed_model", "in_domain",
               "direction", "score", "confidence", "abstain_reason"]


def _write_score_csv(path: Path, variants: list[dict], records: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_FIELDS)
        writer.writeheader()
        for v, r in zip(variants, records):
            confirmed = v.get("confirmed_model") or (r.get("evidence") or {}).get("confirmed_model")
            writer.writerow({"variant_key": r["variant_key"], "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"],
                              "alt": v["alt"], "id": v.get("id"), "confirmed_model": confirmed,
                              "in_domain": r["in_domain"], "direction": r["direction"], "score": r["score"],
                              "confidence": r["confidence"], "abstain_reason": r["abstain_reason"]})


def _write_score_report(path: Path, records: list[dict], input_desc: str, mode: str) -> None:
    scored = [r for r in records if r["in_domain"]]
    row = lambda r: (f"| {r['variant_key']} | {r['in_domain']} | {r['direction'] or '-'} | "
                      f"{r['score'] if r['score'] is not None else '-'} | {r['confidence'] or '-'} | {r['abstain_reason'] or '-'} |")
    lines = ["# ve-regulatory Report", "", f"**Input**: `{input_desc}`", f"**Mode**: {mode}",
              f"**Variants scored**: {len(scored)}", f"**Variants abstained**: {len(records) - len(scored)}", "",
              "## Results", "", "| Variant | In domain | Direction | Score (abs delta log-count, uncalibrated) | Confidence | Abstain reason |",
              "|---|---|---|---:|---|---|", *[row(r) for r in records],
              "", "## What this report does NOT claim", "",
              "- **Not a called peak**: predicted accessibility is CATv1's continuous model output, not an ENCODE ATAC/DNase peak call; no peak BED file was consulted.",
              "- **Not a calibrated effect size**: the magnitude is shown but is explicitly uncalibrated (`evidence.magnitude_uncalibrated`); direction is the supportable claim.",
              "- **No transcription factor is ever named**: a high-attribution, motif-like position is flagged, never identified -- that needs attribution -> seqlets -> a motif database this skill does not run.",
              "", DISCLAIMER, ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_propose_report(path: Path, cell_type: str, rows: list[dict]) -> None:
    lines = ["# ve-regulatory Report", "", "**Mode**: --propose (no scoring; no model was chosen)", f"**Cell-type query**: `{cell_type}`", "",
              f"## Shortlist ({len(rows)} candidates, GRCh38 human DNase-seq/ATAC-seq experiments)", "",
              "| Accession | Biosample | Assay | Name match | count_pearson (fold 0) |", "|---|---|---|---:|---:|",
              *[f"| {r['accession']} | {r['biosample']} | {r['assay']} | {r['name_match']} | {r['count_pearson_fold0']} |" for r in rows],
              "", "No accession was selected automatically. Fuzzy name matching alone is not a safe way to pick a "
              "biological cell type -- review the biosample names above, then re-run with "
              "`--model <accession> --input <variants>` once you have confirmed one.", "", DISCLAIMER, ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_list_report(path: Path, rows: list[dict]) -> None:
    lines = ["# ve-regulatory Report", "", "**Mode**: --list-cell-types", f"**Distinct GRCh38 human biosamples**: {len(rows)}", "",
              "| Biosample | # experiments | Assays |", "|---|---:|---|",
              *[f"| {r['biosample']} | {r['n_experiments']} | {', '.join(r['assays'])} |" for r in rows],
              "", DISCLAIMER, ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(output_dir: Path, result: dict, report_fn, command: list[str]) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report_fn(output_dir / "report.md")
    (output_dir / "reproducibility" / "commands.sh").write_text("#!/usr/bin/env bash\n" + " ".join(command) + "\n", encoding="utf-8")

# --- CLI -----------------------------------------------------------------

def _run_demo(output: Path, command: list[str]) -> int:
    variants = load_variants(DEMO_TSV)
    cache = json.loads(DEMO_CACHE.read_text(encoding="utf-8"))
    by_key = {r["variant_key"]: r for r in cache["records"]}
    records = [by_key[v["variant_key"]] for v in variants]
    result = {"skill": "ve-regulatory", "mode": "demo",
              "provenance": "Cached output from a verified live CATv1/Cherimoya CPU run over these 4 variants "
                             "(see cache_provenance). --demo makes no network or model calls.",
              "cache_provenance": cache.get("provenance"), "records": records, "disclaimer": DISCLAIMER}
    write_outputs(output, result,
                  lambda p: _write_score_report(p, records, str(DEMO_TSV), "Cached demo output (no network, no model call)"), command)
    _write_score_csv(output / "tables" / "results.csv", variants, records)
    print(f"ve-regulatory wrote {output / 'report.md'} (cached demo output, no network)")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="ve-regulatory: CATv1 ref-vs-alt accessibility delta")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("ve_regulatory_out"))
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--propose", action="store_true", help="Rank CATv1 experiments for --cell-type; chooses nothing")
    parser.add_argument("--list-cell-types", action="store_true")
    parser.add_argument("--cell-type", type=str, default=None)
    parser.add_argument("--model", type=str, default=None, help="CATv1 experiment_accession chosen from a --propose shortlist")
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args(argv)
    command = [sys.executable, __file__, *sys.argv[1:]]

    try:
        if args.demo:
            return _run_demo(args.output, command)
        if args.list_cell_types:
            rows = list_cell_types()
            result = {"skill": "ve-regulatory", "mode": "list-cell-types", "biosamples": rows, "disclaimer": DISCLAIMER}
            write_outputs(args.output, result, lambda p: _write_list_report(p, rows), command)
            print(f"ve-regulatory listed {len(rows)} biosamples -> {args.output / 'report.md'}")
            return 0
        if args.propose:
            if not args.cell_type:
                parser.error("--propose requires --cell-type")
            rows = propose(args.cell_type, args.top)
            result = {"skill": "ve-regulatory", "mode": "propose", "cell_type_query": args.cell_type,
                       "shortlist": rows, "model_chosen": None, "disclaimer": DISCLAIMER}
            write_outputs(args.output, result, lambda p: _write_propose_report(p, args.cell_type, rows), command)
            print(f"ve-regulatory proposed {len(rows)} candidates for '{args.cell_type}', chose none -> {args.output / 'report.md'}")
            return 0
        if not args.input:
            parser.error("--input is required unless --demo, --propose, or --list-cell-types is used")
        variants = load_variants(args.input)
        try:
            performance = load_performance()
        except Exception as exc:  # network flake degrades confidence gracefully, never blocks scoring
            print(f"WARNING: could not load CATv1 performance.tsv, confidence will read 'low': {exc}", file=sys.stderr)
            performance = None
        records = [resolve_and_score(v, args, performance) for v in variants]
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    result = {"skill": "ve-regulatory", "mode": "score", "records": records, "disclaimer": DISCLAIMER}
    write_outputs(args.output, result, lambda p: _write_score_report(p, records, str(args.input), "User-provided local data"), command)
    _write_score_csv(args.output / "tables" / "results.csv", variants, records)
    print(f"ve-regulatory wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
