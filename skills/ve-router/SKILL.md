---
name: ve-router
description: >-
  Parses legacy SnpEff EFF and standardised ANN variant-consequence annotations
  from a VCF and assigns exactly one routing class per variant, showing the
  transcript it picked and every alternative it discarded.
license: MIT
metadata:
  version: 0.1.0
  author: ve-skills team
  domain: genomics
  tags:
    - variant-routing
    - consequence-annotation
    - snpeff
    - vcf-parsing
  inputs:
    - name: input_file
      type: file
      format:
        - vcf
        - json
      description: A VCF (gzipped or plain) with EFF/ANN INFO annotations, or a JSON array of Contract A records carrying a raw info field
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Routing report with class counts and consequence-selection transparency
    - name: result
      type: file
      format:
        - json
      description: Contract A records with class set and the full routing decision
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_input.txt
      description: 10 synthetic/real-derived VCF records covering all five routing classes
  endpoints:
    cli: python skills/ve-router/ve_router.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧭"
    homepage: https://github.com/burgshrimps/ve-skills
    os:
      - macos
      - linux
    install: []
    trigger_keywords:
      - variant routing
      - consequence classification
      - parse EFF annotation
      - splice donor acceptor routing
      - frameshift stop gained routing
---

# ve-router

You are **ve-router**, a specialised variant-effect agent that turns a VCF's raw SnpEff/VEP consequence annotations into exactly one routing class per variant, showing its work.

## Trigger

**Fire this skill when the user says any of:**
- "route these variants by consequence"
- "parse the EFF field" / "parse the legacy SnpEff annotation"
- "classify these variants as truncating, splice, missense, or non-coding"
- "which transcript did you pick for this variant"
- "assign a routing class to this VCF"

**Do NOT fire when:**
- The user wants a pathogenicity, ACMG, or diagnostic call (that is `clinical-variant-reporter`, and out of scope for `ve-skills` entirely).
- The user wants segregation/parent-of-origin (`ve-segregation`) or an allele-frequency provenance check (`ve-frequency`).
- The user wants an actual confidence score for a routed variant (`ve-lof`, `ve-splice`, `ve-missense` consume this skill's output; they do not replace it).
- The input has no VCF INFO annotation field and no VCF at all (e.g. a bare gene list).

## Why This Exists

- **Without it**: every downstream branch (`ve-lof`, `ve-splice`, `ve-missense`) would each have to re-parse legacy `EFF`, re-implement the same worst-consequence-per-record pick, and would silently disagree on which of a variant's several transcript annotations to use.
- **With it**: one shared, tested parser produces one routing class per variant, and the transcript-selection decision — including everything that was NOT picked — is emitted as data instead of being thrown away.
- **Why ve-skills, not an existing ClawBio skill**: `rare-high-impact-variants` cannot parse the legacy `EFF` format this challenge data uses (only modern `ANN`/VEP `CSQ`). No other ClawBio skill picks a single consequence out of a multi-transcript, multi-gene annotation set and reports the discarded alternatives.

## Core Capabilities

1. **Legacy EFF parsing**: `Effect(Impact|FunctionalClass|Codon|AA|AA_len|Gene|BioType|Coding|Transcript|Rank|GT)`, comma-separated, one entry per variant x transcript.
2. **Standardised ANN parsing**: modern SnpEff/VEP `Allele|Annotation|Annotation_Impact|Gene_Name|...`, SO-term effect names.
3. **Build-aware annotator fallback**: when a record carries neither `EFF` nor `ANN`, tries a local `snpEff` binary, then the Ensembl VEP REST API on the host that matches the record's genome build — never the wrong one.
4. **Transparent transcript selection**: for every variant, records which annotation was chosen, the exact rule that fired, and the full set of alternatives that were not chosen.
5. **Routing**: exactly one of `protein_truncating | splice | missense | non_coding | unroutable` per variant.

## Scope

One skill, one task: turn annotation entries into one routing class, honestly. It does not score, rank, or judge a variant's consequence (that is `ve-lof`/`ve-splice`/`ve-missense`), does not call ACMG/AMP categories, does not compute segregation or allele frequency, and does not run its own SnpEff/VEP installation — it only calls out to one if the input file lacks annotation entirely.

## Input Formats

| Format | Extension | Required Fields | Example |
|---|---|---|---|
| VCF (plain or gzipped) | `.vcf`, `.vcf.gz`, or any text/gzip file starting with a VCF header | `#CHROM` header line; `INFO` with `EFF` or `ANN` (or none, to exercise the annotator fallback) | `demo_input.txt` |
| Contract A JSON | `.json` | Array of objects, each with `chrom`, `pos`, `ref`, `alt`; an `info` field (string or object) carrying `EFF`/`ANN` from an upstream skill | see SKILL.md "Workflow" |

Format is auto-detected from the file's first 4 KB (gzip-aware): a `##fileformat=VCF` / `#CHROM` header means VCF, a leading `[` or `{` means JSON. Anything else is rejected before any parsing is attempted.

## Workflow

1. **Detect input format and genome build.** VCF: samples from the `#CHROM` line, build from `##reference` (or `--assembly`). JSON: build from each record's `build` field.
2. **Per variant, find the annotation source in this order**: `EFF` present -> parse directly. Else `ANN` present -> parse directly. Else attempt a local `snpEff` binary (skipped if no binary or no known database for the detected build), then the build-matched Ensembl VEP REST host — skipped entirely in `--demo`.
3. **If no annotation could be obtained at all**, the variant is `unroutable` with the reason `"no consequence annotation present and no annotator reachable"`. Nothing is guessed.
4. **Select one annotation** using the rule in Domain Decisions below — never on `FunctionalClass`.
5. **Classify** the selected annotation's effect name into one of the five routing classes.
6. **Emit** the chosen consequence/gene/transcript, the selection rule text, and the full discarded-alternatives list, for every variant.
7. **Write** `report.md`, `result.json`, `tables/results.csv`, `reproducibility/commands.sh`.

## Domain Decisions

| Decision | Rule | Source |
|---|---|---|
| Selection rule | Highest impact tier wins: `HIGH > MODERATE > LOW > MODIFIER` (the tier SnpEff/VEP assigned to that specific annotation); ties broken by the order the annotation appears in the `EFF`/`ANN` field (first transcript listed) | Task brief: "highest impact tier, then first transcript"; verified against all 68 challenge records |
| Routing NEVER uses `FunctionalClass` | Route on the Effect/Annotation name only | Real challenge record `1:11906068:A:G` has `STOP_LOST(HIGH\|MISSENSE\|...)` — `FunctionalClass=MISSENSE` on a stop-lost record. Routing on `FunctionalClass` would misclassify it as `missense`. |
| `protein_truncating` effects | `FRAME_SHIFT`, `STOP_GAINED`, `START_LOST`, `STOP_LOST` (legacy); `frameshift_variant`, `stop_gained`, `start_lost`, `stop_lost` (SO term) | Task brief mapping table |
| `splice` effects | `SPLICE_SITE_DONOR`, `SPLICE_SITE_ACCEPTOR` (legacy); `splice_donor_variant`, `splice_acceptor_variant` (SO term) | Task brief mapping table |
| `missense` effects | `NON_SYNONYMOUS_CODING` (legacy); `missense_variant` (SO term) | Task brief mapping table |
| `non_coding` effects | Known non-protein-altering context terms: `INTRON`, `UPSTREAM`, `DOWNSTREAM`, `EXON`, `UTR_5_PRIME`, `UTR_3_PRIME`, `SYNONYMOUS_CODING`, `INTERGENIC`, and SO equivalents (see `NONCODING_EFF`/`NONCODING_SO` in `ve_router.py`) | Task brief: "everything in a non-coding context"; effect vocabulary is SnpEff's own documented legacy/ANN term set |
| `unroutable` | Effect name recognised in neither list above, or no annotation obtainable at all | Task brief: "anything you cannot confidently classify" / "do not guess" |
| Build detection | `--assembly` flag wins; else scan `##reference` for `grch37`/`g1k_v37`/`hg19`/`b37` -> `GRCh37`, or `grch38`/`hg38` -> `GRCh38`; never assumed | README trap 3: the GRCh38 Ensembl REST host silently returns wrong-build answers for GRCh37 coordinates with no error |
| VEP REST host | `grch37.rest.ensembl.org` for `GRCh37`, `rest.ensembl.org` for `GRCh38` — different hosts, chosen from the detected build, not hardcoded | README trap 3 |
| `--demo` network isolation | The Ensembl VEP REST fallback is never attempted when `--demo` is set | Hard requirement: `--demo` must work with no network |

## CLI Reference

```bash
python skills/ve-router/ve_router.py --input variants.vcf.gz --output /tmp/ve_router_out
python skills/ve-router/ve_router.py --input variants.vcf --assembly GRCh38 --output /tmp/ve_router_out
python skills/ve-router/ve_router.py --demo --output /tmp/ve_router_demo
```

| Flag | Required | Description |
|---|---|---|
| `--input` | Yes, unless `--demo` | Path to a VCF (plain or `.gz`) or Contract A JSON array |
| `--output` | No (default `ve_router_out`) | Output directory |
| `--demo` | No | Run on the packaged `demo_input.txt`, offline |
| `--assembly` | No | Override build detection, e.g. `GRCh37` or `GRCh38` |

## Demo

```bash
python skills/ve-router/ve_router.py --demo --output /tmp/ve_router_demo
```

Expected: 10 variants routed — 5 `protein_truncating`, 2 `splice`, 1 `missense`, 1 `non_coding`, 1 `unroutable` — with no network access.

## Example Output

```markdown
# ve-router Report

**Input**: `/path/to/challenge1-b37-segregation.vcf.gz`
**Mode**: User-provided local data
**Variants routed**: 68

## Routing class counts

| Class | Count |
|---|---:|
| protein_truncating | 56 |
| splice | 12 |
| missense | 0 |
| non_coding | 0 |
| unroutable | 0 |

## Consequence-selection transparency

39 of 68 variants carried more than one variant x transcript annotation; 131
alternative annotations were discarded in favour of the highest-impact-tier
pick. Every discarded alternative is kept in `result.json` under each
record's `routing.discarded` - nothing is silently collapsed.
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

`result.json` records match Contract A (`variant_key`, `chrom`, `pos`, `ref`, `alt`, `id`, `genotypes`, `gene`, `transcript`, `consequence`, `impact`, `class`) plus a `routing` object: `annotation_source` (`EFF`/`ANN`/`snpEff`/`vep_rest`/`none`), `build`, `selection_rule`, `selected`, `discarded` (full list), `unroutable_reason`.

## Dependencies

- Python 3.10+ standard library only (`argparse`, `csv`, `gzip`, `json`, `re`, `shutil`, `subprocess`, `urllib`). No third-party packages.
- Optional, at runtime, only when a record has neither `EFF` nor `ANN`: a `snpEff` binary on `PATH` with a matching local database, or network access to Ensembl's VEP REST API. Neither is required for `--demo` or for the challenge VCF, which is fully EFF-annotated.

## Gotchas

- **Never route on `FunctionalClass`.** The model will want to trust the label that reads `MISSENSE` next to a variant. On real data, a `STOP_LOST` record carries `FunctionalClass=MISSENSE`. Route on the Effect/Annotation name, always.
- **The HIGH-impact annotation SnpEff already flagged is not automatically the answer for a multi-gene record.** A variant can have five annotations across three genes; the selection rule (impact tier, then first-listed transcript) must run every time, not just when it "looks like" there's ambiguity.
- **Do not silently drop alternates.** It is tempting to keep only the winning annotation in the output to save space. Every discarded annotation must survive into `result.json` — that transparency is the entire point of this skill.
- **`--demo` must never make a network call.** The annotator fallback path is wired to check `allow_network=False` under `--demo`; a demo record with no `EFF`/`ANN` must resolve to `unroutable` immediately, not hang on a timeout.
- **Genome build for the VEP fallback is read from the VCF header, never assumed.** Passing GRCh37 coordinates to `rest.ensembl.org` (the GRCh38 host) returns a response with no error — silently wrong, not absent.

## Safety

- **Local-first**: VCF/JSON parsing and EFF/ANN interpretation are entirely local. The only possible outbound call is the build-matched Ensembl VEP REST lookup, and only when a record has no annotation at all and `--demo` is not set.
- **No verdicts, only evidence**: this skill never emits the words rare, pathogenic, diagnostic, de novo, or compound heterozygous. It assigns a routing bucket and shows its work — nothing more.
- **Abstain over guess**: an effect name outside the documented mapping, or a record no annotator could resolve, is marked `unroutable` with a stated reason, never forced into one of the other four classes.
- **Disclaimer**: every `report.md` ends with the ClawBio research-tool disclaimer.

## Agent Boundary

The agent dispatches to this skill and explains its output in plain language. The Python script (`ve_router.py`) does all parsing, selection, and classification; the agent must not re-derive a routing class by eye from raw `EFF` text, and must not describe an `unroutable` result as if it were a determined class.

## Data Attribution

`demo_input.txt` reuses several `EFF` annotation strings verbatim from the ClawBio Berlin 2026 Challenge 1 VCF, a derived subset of the Corpasome by Manuel Corpas ([DOI 10.6084/m9.figshare.693052](https://figshare.com/articles/dataset/Corpasome/693052)), CC BY 4.0. Genotype values in the demo file are illustrative, not the original cohort's calls.
