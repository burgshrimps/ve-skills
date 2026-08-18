---
name: ve-regulatory
description: >-
  Scores whether a non-coding hg38 variant plausibly disrupts a regulatory element by
  running live Cherimoya/CATv1 ref-vs-alt accessibility inference against a caller-specified
  ENCODE cell-type model, and reports honestly which of its three outputs are native model
  results versus interpretive add-ons.
license: MIT
metadata:
  version: 0.1.0
  author: ve-skills team
  domain: genomics
  tags:
    - regulatory-variant
    - chromatin-accessibility
    - cherimoya
    - catv1
    - non-coding
  inputs:
    - name: input_file
      type: file
      format:
        - tsv
        - json
      description: Non-coding hg38 variants as a tab-separated table or a Contract A JSON array, each with an explicit build
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Per-variant accessibility-delta report with an explicit "what this does NOT claim" section
    - name: result
      type: file
      format:
        - json
      description: Contract B records (score/direction/confidence/in_domain/abstain_reason/evidence)
  dependencies:
    python: ">=3.10"
    packages:
      - torch
      - cherimoya
      - bpnet-lite
      - tangermeme
      - huggingface_hub
  demo_data:
    - path: demo_variants.tsv
      description: 4 variants with known biology (rs12740374/SORT1, rs1421085/FTO-IRX3, rs6983267/8q24 MYC enhancer, rs4988235/LCT), hg38 coordinates verified against Ensembl, each with a caller-specified CATv1 model
  endpoints:
    cli: python skills/ve-regulatory/ve_regulatory.py --model {accession} --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/burgshrimps/ve-skills
    os:
      - macos
      - linux
    install:
      - python3 -m venv venv-cat
      - venv-cat/bin/pip install torch
      - venv-cat/bin/pip install --no-deps cherimoya bpnet-lite
      - venv-cat/bin/pip install tangermeme numpy scipy pandas h5py tqdm huggingface_hub
    trigger_keywords:
      - regulatory variant scoring
      - chromatin accessibility delta
      - CATv1 cherimoya
      - non-coding variant accessibility
      - ATAC DNase ref vs alt
---

# ve-regulatory

You are **ve-regulatory**, a specialised variant-effect agent that scores whether a
non-coding hg38 variant plausibly disrupts a regulatory element, using live Cherimoya/CATv1
ref-vs-alt accessibility inference against a cell-type model a human has confirmed.

## Trigger

**Fire this skill when the user says any of:**
- "does this variant affect a regulatory element / enhancer / promoter"
- "score chromatin accessibility change for this variant"
- "run CATv1 / Cherimoya on this non-coding variant"
- "propose a cell type for CATv1" / "what ENCODE experiments match &lt;tissue&gt;"
- "ref vs alt accessibility delta"

**Do NOT fire when:**
- The variant is coding (protein-truncating, missense, splice) — that is `ve-lof`,
  `ve-missense`, `ve-splice`, and this skill will abstain with `"variant is not
  non_coding (wrong branch)"` if handed one.
- The user wants a called ATAC/DNase peak, TF ChIP-seq binding, or a named transcription
  factor — CATv1 gives none of these; see "Why This Exists" and the report's own
  "What this report does NOT claim" section.
- Coordinates are GRCh37/hg19, or the build is not stated — CATv1 is GRCh38-only and this
  skill refuses to guess a build.
- The user wants this skill to pick the cell type itself — it never does; see Domain
  Decisions.

## Why This Exists

- **Without it**: a non-coding variant in this pipeline has nowhere to go — `ve-lof`,
  `ve-missense`, `ve-splice` all require a coding or splice-site consequence.
  ClawBio's `gi-chromatin` and `gi-enhancer` take a bare sequence and ask "what is in
  it?"; neither compares a ref sequence to an alt sequence, so neither can answer "did
  *this variant* change *this element*."
- **With it**: one variant, one caller-specified cell type, one live forward pass through
  the model that ENCODE consortium data actually trained, with the three most tempting
  overclaims (a peak call, a calibrated effect size, a named transcription factor)
  explicitly refused.
- **Why ve-skills, not an existing ClawBio skill**: `gi-chromatin` is a DeepSEA-style
  919-track scan of a single sequence; it never does ref-vs-alt. No other skill in the
  library runs Cherimoya/CATv1 at all.

## Core Capabilities

1. **Human-in-the-loop cell-type discovery** (`--propose`): fuzzy-ranks ENCODE
   DNase-seq/ATAC-seq experiments against a biosample term, shows evidence per row, and
   selects nothing.
2. **Live ref-vs-alt accessibility scoring**: builds a 2114bp hg38 window from Ensembl,
   swaps the alt allele in, and runs both sequences through a caller-specified CATv1
   checkpoint.
3. **Native attribution delta**: `tangermeme.saturation_mutagenesis` around the variant
   position, ref vs alt, on the model's own log-count head.
4. **Honest three-output handling**: predicted accessibility (labelled "not a peak"),
   an uncalibrated magnitude (direction only is asserted), and a motif-like-position flag
   that never names a transcription factor.
5. **Abstention by construction**: eight distinct reasons, each returned as
   `in_domain=false` with a specific `abstain_reason` — never a guess.

## Scope

One skill, one task: score a confirmed-cell-type CATv1 accessibility delta for a single
non-coding SNV. It does not call ATAC/DNase peaks, does not run TF motif discovery, does
not perform liftover, does not handle indels (window construction assumes a
substitution), and never selects a cell-type model on its own.

## Input Formats

| Format | Extension | Required Fields | Example |
|---|---|---|---|
| Tab-separated variant table | `.tsv` | `chrom`, `pos`, `ref`, `alt`; strongly recommended: `build`, `class`, `id`, `confirmed_model` | `demo_variants.tsv` |
| Contract A JSON | `.json` | Array of objects with `chrom`, `pos`, `ref`, `alt`; `build` (top-level or `routing.build`); `class` | ve-router output |

A per-row `confirmed_model` column (or JSON field) lets one invocation score variants
against *different* caller-specified cell types in a single run — the demo file uses this
to score liver, adipose, colon, and small-intestine variants together. Rows without one
fall back to the invocation's `--model` flag, and rows with neither abstain.

## Workflow

1. **Discover** (optional): `--propose --cell-type "<term>"` fetches
   `CATv1-metadata.tsv` and `performance.tsv` from the `programmable-genomics/CATv1`
   Hugging Face repo, fuzzy-ranks GRCh38 human DNase-seq/ATAC-seq experiments, and prints
   a shortlist with accession, biosample name, assay, name-match score, and
   `count_pearson` (fold 0). **Chooses nothing. Exits 0.**
2. **Confirm**: a human reads the shortlist and picks an `experiment_accession`.
3. **Load variants**: parse the TSV or Contract A JSON; resolve each row's model as
   `confirmed_model` (per-row) or `--model` (global).
4. **Validate per variant** (see Domain Decisions for the exact gates and their order):
   class is `non_coding`, build is declared hg38/GRCh38, the ref/alt pair is a single-base
   substitution, the 2114bp window fits the contig, the fetched hg38 base matches the
   declared ref.
5. **Fetch** the 2114bp hg38 window from the Ensembl REST sequence API, retrying once on
   timeout.
6. **Score**: one-hot encode ref and alt windows, forward-pass both through the confirmed
   Cherimoya checkpoint (`device="cpu", compile=False`) to get `y_profile, y_counts`.
7. **Attribute**: run `tangermeme.saturation_mutagenesis` (wrapped in
   `cherimoya.ControlWrapper` + `LogCountWrapper`) over a ±25bp window around the variant;
   per-base importance = `(attr * X_onehot).sum(dim=1)`.
8. **Emit** Contract B: `score` = `|Δlog-count|`, `direction` from a fixed noise-floor
   threshold, `confidence` from the model's own `count_pearson`, and `evidence` carrying
   the three honestly-labelled outputs plus their disclaiming notes.
9. **Write** `report.md`, `result.json`, `tables/results.csv`,
   `reproducibility/commands.sh`.

## Domain Decisions

| Decision | Rule | Source |
|---|---|---|
| Window construction | `pos - 1057 .. pos + 1056` inclusive = 2114bp, variant at 0-based offset 1057 | Task spec: CATv1 input shape is `(N,4,2114)`; symmetric split of `2114` around one variant base is `1057` bases on the shorter (left) side and `1056` on the longer (right) side |
| Contig-end abstention | `pos - 1057 < 1` aborts before any network call | Task spec's own abstention list, "variant within 1057bp of a contig end" |
| Cell-type selection | Never automatic. `--propose` only ranks and prints; scoring requires `--model` or a per-row `confirmed_model` | Task spec: "The skill must NEVER silently pick a model," verified concretely below |
| Fuzzy-match scoring | Exact match = 1.0; substring match ≈0.85–1.0; else the best `difflib.SequenceMatcher` ratio between any query word and any name word, ×0.8 | Own implementation, stdlib-only (no synonym dictionary) — **deliberately not tuned to "fix" biology**, see below |
| **Verified fuzzy-match failure #1** | Query `"colorectal"` against the live CATv1 metadata ranks `ENCSR000EOK` ("renal cortical epithelial cell", 0.533) **above** `ENCSR994KTY` ("colonic mucosa", 0.471) — a kidney line outranks the correct colon tissue purely from the substring `cortical`/`colorectal` overlap | Live `--propose --cell-type colorectal` run against `programmable-genomics/CATv1`, 2026-08-18 |
| **Verified fuzzy-match failure #2** | Query `"intestine"` scores `"large intestine"` (0.940) and `"small intestine"` (0.940) within 0.06 of each other; lactase (`LCT`, the demo's own variant) is a small-intestinal brush-border enzyme, so the tissues are not interchangeable despite the near-identical score | Live `--propose --cell-type intestine` run, same date; standard physiology (lactase expression is enterocyte/small-intestine specific) |
| Model checkpoint | `hf_hub_download("programmable-genomics/CATv1", f"models/{accession}/cherimoya.fold_{fold}.torch")`, `fold=0` default | CATv1 repo layout (`manifest.csv`), verified live |
| Confidence bins | `count_pearson_fold0 >= 0.85` → high, `>= 0.60` → medium, else low | CATv1's own `performance.tsv`; thresholds are ours (round, memorable cut points on a Pearson scale), not a published clinical cutoff |
| Direction threshold | `\|Δlog-count\| > 0.1` natural-log units → directional; else "no meaningful change" | Our own noise floor, **not a significance test** — chosen only to avoid reporting a direction on numerical noise near zero |
| Magnitude labelling | Always shown, always tagged `UNCALIBRATED`, `fold_change_approx` computed via `math.expm1` on each side (matches `cherimoya.ExpectedCountsWrapper`'s own count-recovery method) | Task spec: published evaluations of this model class report predicted effects roughly an order of magnitude below measured ones for rs12740374/SORT1 (~12-fold measured). Our own live run on that exact variant predicted 2.13-fold — same direction, ~6x understated |
| Attribution window | ±25bp around the variant (51bp, 153 single-base edits per sequence) | CPU cost tradeoff: a full 2114bp saturation-mutagenesis run measured ~80s for *two* sequences even at a 101bp sub-window; ±25bp keeps a 4-variant `--demo`-equivalent run under ~2 minutes on this reference machine (Apple Silicon, CPU) while still covering local motif context |
| Attribution batch size | 8 | Measured on this machine: CPU forward-pass throughput is *not* monotonic in batch size for this checkpoint (batch 4–8 ≈ 47ms/sequence; batch 256 ≈ 187ms/sequence) — small batches were fastest, the opposite of typical GPU intuition |
| Motif-like flag | `\|attribution at variant\| >= 0.5 × max(\|attribution\|)` within the ±25bp window | Our own heuristic, no published significance threshold exists; intentionally relative to the variant's own local window, never to a genome-wide baseline |
| TF identity | Never computed, never reported | Task spec: naming a TF needs attribution → seqlets → a motif database (e.g. tomtom vs JASPAR); this skill does not run that pipeline |
| Build check | Only literal `hg38`/`grch38` (case-insensitive) passes; anything else, or absent, abstains | README trap 3: a GRCh38-only model silently returns wrong-build answers for mismatched coordinates with no error, so the build must be explicit, never assumed |
| Ensembl host | `rest.ensembl.org` only (GRCh38 host) | CATv1 verified facts: "GRCh38/hg38 ONLY, no liftover in the package" |
| Ensembl retry | One retry on timeout/URLError before abstaining | Empirically reproduced twice on this machine: a first request can hang indefinitely with zero bytes, a retry with the same 20s timeout succeeds immediately |

## CLI Reference

```bash
# 1. Discover — never scores, never chooses
python skills/ve-regulatory/ve_regulatory.py --propose --cell-type "liver" --output /tmp/ve_reg_propose

# 2. List every available biosample
python skills/ve-regulatory/ve_regulatory.py --list-cell-types --output /tmp/ve_reg_list

# 3. Score, one confirmed model for every row
python skills/ve-regulatory/ve_regulatory.py --model ENCSR562FNN --input variants.tsv --output /tmp/ve_reg_out

# 4. Score, per-row confirmed_model column (multi-tissue in one run)
python skills/ve-regulatory/ve_regulatory.py --input variants_with_confirmed_model.tsv --output /tmp/ve_reg_out

# 5. Cached demo, no network, no model stack required
python skills/ve-regulatory/ve_regulatory.py --demo --output /tmp/ve_reg_demo
```

| Flag | Required | Description |
|---|---|---|
| `--input` | Yes, unless `--demo`/`--propose`/`--list-cell-types` | TSV or Contract A JSON |
| `--output` | No (default `ve_regulatory_out`) | Output directory |
| `--demo` | No | Cached 4-variant demo; no network, no torch/cherimoya required |
| `--propose` | No | Requires `--cell-type`; prints a shortlist and scores nothing |
| `--list-cell-types` | No | Dumps every distinct GRCh38 biosample in CATv1 |
| `--cell-type` | No | Free-text biosample query for `--propose`, or to tag an unconfirmed run |
| `--model` | No | Caller-specified `experiment_accession`, applies to rows without their own `confirmed_model` |
| `--fold` | No (default `0`) | Which of the 5 trained folds to load |
| `--top` | No (default `10`) | Shortlist size for `--propose` |

## Demo

```bash
python skills/ve-regulatory/ve_regulatory.py --demo --output /tmp/ve_reg_demo
```

Expected: 4 variants, all `in_domain=true`, reproduced verbatim from a real live CPU run
(`tests/fixtures/cached_scores.json`, see its `provenance` block) — rs12740374/SORT1
(liver, `ENCSR562FNN`) increases predicted accessibility (Δlog-count 0.75, ~2.1-fold
uncalibrated); rs4988235/LCT (small intestine, `ENCSR133KBX`) increases modestly (0.13,
~1.15-fold); rs1421085/FTO-IRX3 (adipose, `ENCSR540BML`) and rs6983267/8q24 MYC enhancer
(colon, `ENCSR994KTY`) both show no meaningful change at their own variant position. No
network call and no torch/cherimoya import happens in this mode.

## Example Output

```markdown
# ve-regulatory Report

**Mode**: Cached demo output (no network, no model call)
**Variants scored**: 4

## Results

| Variant | In domain | Direction | Score (abs delta log-count, uncalibrated) | Confidence | Abstain reason |
|---|---|---|---:|---|---|
| 1:109274968:G:T | True | increases | 0.7531 | high | - |
| 16:53767042:T:C | True | no meaningful change | 0.0075 | medium | - |

## What this report does NOT claim

- **Not a called peak**: predicted accessibility is CATv1's continuous model output, not
  an ENCODE ATAC/DNase peak call; no peak BED file was consulted.
- **Not a calibrated effect size**: the magnitude is shown but is explicitly uncalibrated;
  direction is the supportable claim.
- **No transcription factor is ever named**.
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   └── results.csv
└── reproducibility/
    └── commands.sh
```

`result.json` records match Contract B exactly: `skill`, `variant_key`, `score`,
`direction`, `confidence`, `in_domain`, `abstain_reason`, `evidence`. `evidence` (when
`in_domain=true`) carries `confirmed_model`, `fold`, `count_pearson_fold0`, `window`,
`predicted_accessibility`, `magnitude_uncalibrated`, `contribution_delta`,
`tf_binding_high_attribution_motif_like_position` (boolean, never a TF name), and
`notes` (the disclaiming sentences for each of the three outputs above).

## Dependencies

- Python 3.10+ (developed and verified on 3.13).
- **Required for live scoring only** (`--demo` needs none of these): `torch`,
  `cherimoya` (`--no-deps`), `bpnet-lite` (`--no-deps`), `tangermeme`,
  `huggingface_hub`. See `--help` output / stderr for the exact install recipe; it is
  also reproduced in Gotchas below.
- Network, for live scoring and for `--propose`/`--list-cell-types`: Hugging Face Hub
  (`programmable-genomics/CATv1`) and Ensembl REST (`rest.ensembl.org`).
- No GPU required. Verified CPU-only on Apple Silicon (M-series, arm64).

## Gotchas

- **`cherimoya`'s real dependency list is incomplete.** `pip install --no-deps
  cherimoya` succeeds, but `import cherimoya` then fails with
  `ModuleNotFoundError: No module named 'bpnetlite'` — an unguarded top-level import in
  `cherimoya/losses.py`. `pip install --no-deps bpnet-lite` fixes it. `triton` is
  imported but wrapped in a `try/except ImportError` inside `cherimoya/cheri.py` and is
  never needed for CPU inference — do not chase a triton wheel that does not exist for
  macOS.
- **The Ensembl REST API can hang on the first request and succeed instantly on a
  retry.** Reproduced twice on this machine, independently, months apart. Always pass an
  explicit timeout (`fetch_window` uses 20s) and always retry once before treating it as
  a real failure — treating the first hang as a permanent outage will make this skill
  needlessly flaky in the field.
- **Batch size does not behave like GPU intuition on CPU.** Measured on this reference
  machine: batch 4–8 gave the *lowest* per-sequence time; batch 256 was ~4x slower per
  sequence. If you change `ATTR_BATCH_SIZE`, re-measure — do not assume "bigger batch is
  faster."
- **A high fuzzy-match score is not a biology check.** The model will want to trust the
  top row of a `--propose` shortlist. Verified live: for `"colorectal"` the top row by
  name-match is a kidney cell line, not colon tissue (see Domain Decisions). Always
  present the full shortlist and require an explicit `--model`/`confirmed_model`; never
  auto-select row 1.
- **`FunctionalClass`-style shortcuts don't apply here, but a parallel trap does**: two
  biosample names can score nearly identically (`"large intestine"` vs `"small
  intestine"`) while being different biology for the variant in question. The shortlist
  score is a name-similarity signal, not a tissue-correctness signal.
- **The count head is `log(count+1)`, not counts.** Recovering a fold-change needs
  `math.expm1` on each side (matching `cherimoya.ExpectedCountsWrapper`), not a naive
  `exp()` of the log-count delta — the two diverge whenever counts are not ≫1.

## Safety

- **Never emits** rare, pathogenic, diagnostic, de novo, or compound heterozygous as a
  claim about a variant.
- **Every `report.md` carries an explicit "What this report does NOT claim" section**
  naming the three overclaim risks (peak call, calibrated magnitude, named TF) so a
  reader cannot miss them even skimming.
- **No silent model selection.** `--propose` never scores; scoring without a confirmed
  model abstains with `"cell-type model proposed but not confirmed"` or `"no --cell-type
  given"`.
- **No fabrication on a failed environment gate.** Missing `torch`/`cherimoya` exits 2
  with the install recipe on stderr — it never falls back to a fabricated score.
- **Disclaimer**: every `report.md` ends with the ClawBio research-tool disclaimer.

## Agent Boundary

The agent dispatches to this skill, relays the shortlist for a human decision at
`--propose`, and explains `report.md` in plain language. The Python script does all
sequence fetching, model inference, and abstention logic; the agent must never itself
pick a `--model` accession from a `--propose` shortlist, must never repeat a
`fold_change_approx` without its uncalibrated caveat, and must never name a
transcription factor that this skill's output did not name.

## Data Attribution

`demo_variants.tsv` uses four literature-known dbSNP variants (rs12740374, rs1421085,
rs6983267, rs4988235); hg38 coordinates and reference alleles were looked up live against
Ensembl's variation and sequence REST APIs on 2026-08-18, not assumed from memory. CATv1
model checkpoints and metadata are from the `programmable-genomics/CATv1` Hugging Face
repository (ENCODE consortium DNase-seq/ATAC-seq experiments). See the root
[`README.md`](../../README.md) for the Corpasome attribution covering the rest of the
`ve-skills` pipeline.
