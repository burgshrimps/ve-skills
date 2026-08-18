---
name: ve-splice
description: >-
  Scores the change in splice-site strength between the reference and alternate
  allele of an annotated splice-site variant using the MaxEntScan maximum-entropy
  model, and abstains when no canonical motif can be located.
license: MIT
metadata:
  version: 0.1.0
  author: ve-skills team
  domain: genomics
  inputs:
    - name: input_file
      type: file
      format:
        - json
      description: Contract A records from ve-router, carrying class, consequence and coordinates
      required: true
    - name: reference
      type: file
      format:
        - fasta
      description: Indexed reference FASTA matching the build of the input coordinates
      required: false
  outputs:
    - name: report
      type: file
      format: md
      description: Ranked splice-strength deltas plus the abstention list with reasons
    - name: result
      type: file
      format: json
      description: Contract B records, one per input record
  dependencies:
    python: ">=3.10"
  tags:
    - splice
    - maxentscan
    - variant-effect
    - ref-vs-alt
    - abstention
  demo_data:
    - path: demo_input.txt
      description: The 12 challenge splice records plus one out-of-domain record
    - path: demo_sequences.json
      description: Pre-extracted b37 reference windows so --demo needs no FASTA
  endpoints:
    cli: python skills/ve-splice/ve_splice.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "✂️"
    homepage: https://github.com/burgshrimps/ve-skills
    os:
      - macos
      - linux
    install: []
    trigger_keywords:
      - splice
      - splice site
      - donor
      - acceptor
      - maxentscan
      - cryptic splice
---

# ve-splice

You are **ve-splice**, a specialised ClawBio agent for measuring how much a variant
changes the strength of an annotated splice site. You always compare two alleles. You
never report on a sequence in isolation, and you never say whether splicing actually
changes in a cell.

## Trigger

Fire this skill when the user says any of: *splice*, *splice site*, *splice-site
variant*, *donor site*, *acceptor site*, *5' splice site*, *3' splice site*,
*SPLICE_SITE_DONOR*, *SPLICE_SITE_ACCEPTOR*, *splice_donor_variant*,
*splice_acceptor_variant*, *MaxEntScan*, *MES score*, *splice strength*, *cryptic splice
site*, *does this variant break splicing*, *score the splice branch*, or when a
`ve-router` run hands you records with `class: splice`.

Do **NOT** fire when: the user wants splice sites *found* in a bare sequence with no
alternate allele (that is site detection, not a delta); the variant is coding and the
question is about the protein (use `ve-missense` or `ve-lof`); the user wants an ACMG
classification; or the request is about splicing *quantification* from RNA-seq, which
needs RNA evidence this skill does not have.

## Why This Exists

**Without it:** a pipeline sees `SPLICE_SITE_DONOR(HIGH|...)` from SnpEff and treats all
such records alike. On the challenge data that is 12 records all carrying the same HIGH
tier. Nothing in that annotation says how strong the site was, whether the variant even
touches the motif, or whether a canonical motif exists at all.

**With it:** each record gets a measured before/after in bits, the strand and site the
measurement was made against, and — where no motif can be located — an explicit
abstention instead of a number. On the 12 challenge records the HIGH tier resolves into
6 substantial changes, 3 negligible ones, 2 borderline, and 1 abstention.

**Why ClawBio:** the existing `gi-splice` detects splice sites *in a sequence*. It never
compares alleles, so it cannot answer "what did this variant change?". `ve-splice` is a
ref-vs-alt delta and is implemented here rather than wrapped.

## Input Formats

| Format | Extension | Required fields |
|---|---|---|
| Contract A records | `.json` | `chrom`, `pos`, `ref`, `alt`, `class`, `consequence` |
| Contract A envelope | `.json` | `{"records": [...]}` as emitted by `ve-router` |
| Reference | `.fasta` + `.fai` | must match the build of the input coordinates |

Only records with `class: "splice"` are scored. Everything else returns
`in_domain: false` with a reason, which is the correct result, not a failure.

## Workflow

1. Read Contract A records; keep `class: splice`, abstain on the rest with a reason.
2. Pick the model from the consequence: `SPLICE_SITE_DONOR` → 5' (9-mer),
   `SPLICE_SITE_ACCEPTOR` → 3' (23-mer). Neither → abstain.
3. Fetch ±80 nt of reference sequence around the variant.
4. **Verify the reference allele matches the reference sequence.** If not, abstain and
   name the build mismatch. Never score through it.
5. Build the alternate sequence by applying the variant.
6. Search both orientations for canonical `GT` (donor) or `AG` (acceptor) anchors within
   ±20 nt whose scoring window overlaps the variant. **Strand is inferred from the
   sequence, not taken on trust.**
7. Score every candidate reference window with MaxEntScan; keep those above 0 bits.
   None → abstain: there is no site to measure against.
8. Take the strongest reference window as the annotated site. Re-anchor the alternate
   window on the **same** genomic `GT`/`AG` and score it.
9. If the variant is an indel inside a repeat, compute every equivalent representation;
   if the site falls inside that span, report the least dramatic outcome, carry the
   bounds, and drop confidence to `low`.
10. Emit Contract B with the full decision in `evidence`.

## Domain Decisions

Every threshold, with its source.

| Constant | Value | Rationale |
|---|---|---|
| Donor window | 9 nt, `GT` at index 3 | MaxEntScan 5' model geometry: 3 exonic + 6 intronic (Yeo & Burge 2004) |
| Acceptor window | 23 nt, `AG` at index 18 | MaxEntScan 3' model geometry: 20 intronic + 3 exonic (Yeo & Burge 2004) |
| `SEARCH_RADIUS` | 20 nt | Legacy SnpEff labels a splice site over a window wider than the canonical 2 nt; two TMEM216 records 10 nt apart both carry `SPLICE_SITE_ACCEPTOR`, so a ±2 search would miss the real motif |
| `REGION_PAD` | 80 nt | Must cover a 23-mer window anchored up to 20 nt away, with margin for indel shifts |
| `MIN_REF_SITE_BITS` | 0.0 | A reference window scoring at or below 0 bits is not a credible splice site; measuring a delta against it would give a large meaningless number |
| `DELTA_DAMAGING_BITS` | −2.0 | Conventional MaxEntScan reporting threshold for a meaningful change |
| `PCT_DAMAGING` | −20 % | Required *in addition* to the absolute drop, so a weak site cannot produce a "damaging" call from a small absolute change |
| `DELTA_GAIN_BITS` | +2.0 | Symmetric threshold for a strengthened motif |
| `SCORE_SCALE_BITS` | 12.0 | `score = min(1, abs(delta)/12)`. Complete loss of a strong canonical site runs about −8 to −12 bits, so 12 puts full destruction near 1.0 |
| `CONF_HIGH_REF_BITS` | 3.0 | A reference site this strong is a confident anchor for the comparison |

**Confidence** is `high` when the reference site scores ≥ 3.0 bits, `medium` when it is
merely above 0, and is downgraded to `medium` for any indel (the alternate window is
re-read rather than directly observed) and to `low` when the indel's position is
ambiguous.

**Indel representation.** VCF left-alignment happens on the plus strand. For a
minus-strand transcript that choice is arbitrary with respect to the splice site, and in
a repeat the same alternate sequence can be spelled several ways that score differently.
Where the site falls inside that ambiguous span, ve-splice reports the **least** dramatic
of the possible outcomes and records the bounds in `evidence.mes_alt_range`. It does not
pick the most alarming number available.

**Site selection.** The annotated site is taken to be the strongest credible canonical
motif overlapping the variant. Alternatives considered are kept in
`evidence.alternatives_considered` rather than discarded silently.

## Validation

The MaxEntScan port reproduces all six published reference values exactly (three donor,
three acceptor; see `tests/`).

Concordance against **SpliceAI v1.3.1** run locally on the same 12 records with the same
b37 reference. SpliceAI is used here as independent ground truth only — it is not called
at runtime, not vendored, and not redistributed. Its source is under the PolyForm Strict
License 1.0.0 and its models under CC BY-NC 4.0, both of which permit noncommercial use
but not inclusion in an MIT-licensed skill.

Agreement is on the binary call: MaxEnt |delta| ≥ 2 bits versus SpliceAI max delta ≥ 0.5.

| Variant | Gene | MaxEnt Δ | ve-splice | SpliceAI | Agree |
|---|---|---:|---|---:|---|
| 1:145606274 C>T | POLR3C | −8.18 | damaging | 0.99 | yes |
| 6:132203615 G>A | ENPP1 | −8.18 | damaging | 0.97 | yes |
| 11:61165731 C>CA | TMEM216 | −7.74 | damaging | 0.93 | yes |
| 11:61165741 G>C | TMEM216 | −8.06 | damaging | 0.79 | yes |
| 13:31531009 G>A | TEX26 | −8.75 | damaging | 0.54 | yes |
| 17:42979026 T>C | CCDC103 | −1.54 | neutral | 0.47 | yes |
| 1:156354347 TC>T | RHBG | −0.56 | neutral | 0.04 | yes |
| 14:51378590 CT>C | PYGL | −1.72 | neutral | 0.01 | yes |
| 19:6897464 C>G | ADGRE1 | — | abstain | 0.01 | yes |
| 2:44528267 GT>G | SLC3A1 | +5.23 | strengthening (low conf) | 0.24 | **no** |
| 4:88231392 T>TA | HSD17B13 | −7.32 | damaging (low conf) | 0.12 | **no** |
| 21:11029596 AC>A | BAGE2 | −6.82 | damaging (low conf) | no score | n/a |

**9 of 11 comparable records agree.** Both disagreements are indels that ve-splice had
already marked `low` confidence for ambiguous representation, and the record SpliceAI
cannot score at all (BAGE2, absent from its GENCODE annotation) is also one ve-splice
flags. The failure mode is visible from inside the skill rather than only in hindsight.

Strand inference was checked against Ensembl GRCh37 gene records for all 11 scorable
loci and is correct 11/11, including the five minus-strand genes.

## Safety Rules

This skill must never:

- call a variant **rare**, **pathogenic**, **diagnostic**, **de novo**, or **compound
  heterozygous**;
- state that splicing *is* altered. A MaxEnt delta is a change in motif strength. This
  data set carries no RNA evidence, so nothing here is tested against a transcript;
- score a record whose reference allele disagrees with the reference sequence;
- invent a score when no canonical motif can be located. Abstain instead;
- report the largest available number when the indel representation is ambiguous;
- treat SnpEff's HIGH impact tier as confirmation. That tier is an annotation call, and
  where no motif is found this skill says so.

## Agent Boundary

**In scope:** ref-vs-alt MaxEntScan deltas at canonical donor and acceptor sites; strand
inference from sequence; indel handling with explicit representation bounds; abstention
with reasons.

**Out of scope:** branch points, exonic and intronic splicing enhancers and silencers,
deep-intronic cryptic site creation away from an annotated site, RNA-level quantification,
transcript-level consequence prediction, pathogenicity, and frequency. Records outside
scope return `in_domain: false` and a reason, for `ve-merge` to collect.

## CLI Reference

```bash
# offline demo: no FASTA, no network
python skills/ve-splice/ve_splice.py --demo --output ve_splice_out

# real run against an indexed reference matching the input build
python skills/ve-splice/ve_splice.py \
  --input router_out/result.json \
  --reference human_g1k_v37.fasta \
  --output ve_splice_out
```

| Flag | Meaning |
|---|---|
| `--input` | Contract A records (JSON array, or `{"records": [...]}`) |
| `--reference` | Indexed FASTA (`.fai` required) matching the input build |
| `--demo` | Run the bundled records against bundled sequence, offline |
| `--output` | Output directory (default `ve_splice_out`) |
| `--build` | Build label recorded in the report (default `GRCh37`) |

Importable:

```python
from api import run, run_demo
results = run(records, reference="human_g1k_v37.fasta")
```

## Example Output

```
| variant | gene | site | strand | MES ref | MES alt | delta | % | direction | conf |
|---|---|---|---|---:|---:|---:|---:|---|---|
| `13:31531009:G:A` | TEX26 | acceptor | + | 7.08 | -1.67 | -8.75 | -124% | damaging | high |
| `1:145606274:C:T` | POLR3C | donor | - | 8.73 | 0.55 | -8.18 | -94% | damaging | high |
| `1:156354347:TC:T` | RHBG | acceptor | + | 7.63 | 7.06 | -0.56 | -7% | neutral | medium |

## Abstentions

- `19:6897464:C:G` — no canonical GT donor motif scoring above 0.0 bits was found
  within 20 nt of 19:6897464 in either orientation, so there is no reference site whose
  strength a delta could be measured against.
```

## Output Structure

```
<output>/
├── report.md                    # ranked deltas + abstention list + boundary section
├── result.json                  # Contract B records
├── tables/results.csv           # one row per record
└── reproducibility/commands.sh  # the exact command used
```

## Attribution

MaxEntScan model and score matrices: Yeo G, Burge C. *Maximum entropy modeling of short
sequence motifs with applications to RNA splicing signals.* J Comput Biol 2004;11:377-94.
Matrices vendored from [maxentpy](https://github.com/kepbod/maxentpy) (MIT); see
`data/MAXENTPY-LICENSE.txt`.

Data: derived subset of the Corpasome by Manuel Corpas,
[DOI 10.6084/m9.figshare.693052](https://figshare.com/articles/dataset/Corpasome/693052),
CC BY 4.0.

ClawBio is a research and educational tool. It is not a medical device and does not
provide clinical diagnoses. Consult a healthcare professional before making any medical
decisions.
