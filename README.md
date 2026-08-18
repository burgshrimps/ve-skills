# ve-skills — a self-contained variant-effect layer

A series of skills that take a **filtered VCF with pedigree genotypes** and produce a
**prioritised variant list plus an explicit abstention list**, routing each variant to the
predictor that is actually entitled to an opinion about it.

Built for ClawBio Berlin 2026, Challenge 1 ("End the diagnostic odyssey").

**Every `ve-*` skill is implemented by us and runs standalone.** We do not call `gi-splice`,
`gi-chromatin`, `clinical-variant-reporter` or any other existing ClawBio skill at runtime.
No hosted-API dependency, nothing to break in the demo.

---

## The pipeline

```
challenge1 VCF (68, b37) ─┐
outside demo sets ────────┤
                          ▼
              ┌────────────────────┐
              │  ve-segregation    │  pedigree genotypes → parent-of-origin
              └─────────┬──────────┘  30 paternal / 38 maternal
                        ▼
                 ┌──────────────┐
                 │ ve-frequency │   build-aware AF provenance
                 └──────┬───────┘   → rare / common / NO-DATA
                        ▼
                 ┌──────────────┐
                 │  ve-router   │   one class per variant × transcript
                 └──────┬───────┘
      ┌──────────┬──────┴──────┬──────────────┬─────────────┐
      ▼          ▼             ▼              ▼             ▼
   ve-lof   ve-missense   ve-splice    ve-regulatory   unroutable
   56 recs    0 recs       12 recs        0 recs        → abstain
      └──────────┴──────┬──────┴──────────────┴─────────────┘
                        ▼
                 ┌──────────────┐
                 │   ve-merge   │   + segregation evidence
                 └──────┬───────┘
              ranked list + abstention list
```

Record counts are what the **challenge pack actually produces**. The two zero-count
branches are demoed on outside public data (ClinVar, gnomAD), which the challenge rules
allow.

---

## Why this doesn't overlap with ClawBio's 98 skills

Everything predictive already in ClawBio — `gi-splice`, `gi-chromatin`, `gi-enhancer`,
`gi-promoter`, `gi-expression` — takes **a sequence** and asks *"what is in it?"* None
takes a VCF. None does ref-vs-alt.

> **Our skills take a variant and ask "what changed?"**

That is the distinction, and it holds for every branch. We implement the models ourselves
rather than wrapping theirs, so the answer to "didn't ClawBio already have this?" is both
*different question* and *different code*.

Specifically not duplicated, and specifically not called:

| Existing | Why ours is different |
|---|---|
| `gi-splice` | Detects splice sites in a sequence. Never compares alleles. |
| `gi-chromatin` | DeepSEA-style 919-track scan of a sequence. Not a variant delta. |
| `clinical-variant-reporter` | ACMG/AMP classification. We produce evidence, not ACMG codes. |
| `rare-high-impact-variants` | Classifies rare/common/unknown from AF **already in the file**. Cannot fetch, cannot parse legacy `EFF`. |
| `cnv-acmg-classifier` | CNV/SV. We are SNV+indel only. |
| *nothing* | **No skill in the library does pedigree segregation.** Verified by grep across all 98 `SKILL.md` files. |

---

## Shared contracts

These three shapes are the whole reason three people can build in parallel without
blocking each other. **Agree on these before writing code; don't change them unilaterally.**

### Contract A — the variant record (internal currency)

Built up by `ve-segregation` → `ve-frequency` → `ve-router`, consumed by every branch.

```json
{
  "variant_key": "1:11906068:A:G",
  "chrom": "1", "pos": 11906068, "ref": "A", "alt": "G", "id": "rs5065",
  "genotypes": { "ISDBM322015": "0/1", "ISDBM322016": "0/1",
                 "ISDBM322017": "0/0", "ISDBM322018": "0/0" },
  "segregation": {
    "pattern": "paternal",
    "carriers": ["ISDBM322015", "ISDBM322016"],
    "phased": false,
    "rule": "proband carries AND exactly one parent carries"
  },
  "gene": "NPPA",
  "transcript": "NM_006172.3",
  "consequence": "STOP_LOST",
  "impact": "HIGH",
  "freq": {
    "af": null, "source": null, "build": null,
    "class": "NO_DATA",
    "provenance_warning": "INFO/AF present but is GATK cohort AF over 4 samples, not population AF"
  },
  "class": "protein_truncating",
  "routing": {
    "annotation_source": "EFF",
    "build": "GRCh37",
    "selection_rule": "Highest annotation impact tier wins (HIGH > MODERATE > LOW > MODIFIER); ties broken by order in the EFF/ANN field. The Effect name decides the class - FunctionalClass is never used for routing.",
    "selected":  { "effect": "STOP_LOST", "gene": "NPPA", "transcript": "NM_006172.3",
                   "impact": "HIGH", "functional_class": "MISSENSE" },
    "discarded": [ { "effect": "DOWNSTREAM", "gene": "CLCN6", "transcript": "NM_001256959.1",
                     "impact": "MODIFIER", "functional_class": null } ]
  }
}
```

`routing` is written by `ve-router` and is the audit trail for the `class` decision:
which annotation won, under what rule, from what source and build, and every alternative
that was discarded. Branch skills can ignore it; `ve-merge` uses it to show the working.
`annotation_source` is one of `EFF` | `ANN` | `snpEff` | `vep_rest`.

- `segregation.pattern` ∈ `paternal` | `maternal` | `ambiguous` | `excluded`
- `freq.class` ∈ `RARE` | `COMMON` | `NO_DATA`
- `class` (routing) ∈ `protein_truncating` | `missense` | `splice` | `non_coding` | `unroutable`

### Contract B — branch output (every branch returns exactly this)

```json
{
  "skill": "ve-lof",
  "variant_key": "1:11906068:A:G",
  "score": 0.82,
  "direction": "damaging",
  "confidence": "medium",
  "in_domain": true,
  "abstain_reason": null,
  "evidence": { "nmd_escape": false, "last_exon": false, "flags": ["..."] }
}
```

**Why this matters:** when a branch is handed something outside its domain it returns
`in_domain: false` with an `abstain_reason`, e.g.

> `ve-missense` on a frameshift → `"AlphaMissense scores single-amino-acid substitutions
> arising from SNVs; indels are outside its domain."`

`ve-merge` collects every `in_domain: false` into the abstention list. **The abstention
list is generated by construction, not written by hand.** That is the stretch goal of the
challenge and the thing it says it scores highest.

### Contract C — pedigree

`ve-segregation` must **derive** roles, never hardcode them. See the trap below.

```json
{ "proband": "ISDBM322015", "father": "ISDBM322016",
  "mother": "ISDBM322018", "siblings": ["ISDBM322017"] }
```

---

## Verified facts about the challenge data

All independently checked against the file, not taken from the brief. Checksums match the
published SHA-256 values.

| Property | Value |
|---|---|
| Records | 68, GRCh37/b37, all PASS, all autosomal, all biallelic |
| Variant types | **35 SNVs, 33 indels** (max 4 bp), **0 CNVs/SVs** |
| Protein-truncating | **56** — 27 frameshift, 21 stop-gained, 5 start-lost, 3 stop-lost |
| Splice-site | **12** — 7 donor, 4 acceptor, 1 both |
| Missense | **0** as a selected consequence (11 `NON_SYNONYMOUS_CODING` exist, but only on alternate transcripts at MODERATE) |
| Segregation | 30 paternal / 38 maternal, reproduced 68/68 from genotypes |
| Quality | every sample call DP ≥ 10, GQ ≥ 20 |

---

## Four traps that produce silently wrong answers

These are the substance of the project. Each is a real failure mode we can demonstrate.

### 1. The brief's sample roles are wrong

The challenge text says *"ISDBM322015 to ISDBM322018 are son, father, mother, sister."*
**Mother and sister are swapped.**

| Sample | Actual role |
|---|---|
| ISDBM322015 | son |
| ISDBM322016 | father |
| ISDBM322017 | **sister** |
| ISDBM322018 | **mother** |

The TSV column headers say so, and those columns are 68/68 genotype-identical to the VCF
sample columns. Proof by outcome: the data mapping gives **30/38 with 68/68 agreement**
against the `PARENT_OF_ORIGIN_UNPHASED` labels; the brief's mapping gives 11/25/32 and
36/68. **Derive roles from the data. Do not hardcode from the prose.**

### 2. `INFO/AF` is not population frequency

The VCF's INFO keys are `AC, AF, AN, MLEAC, MLEAF, DP, EFF, VQSLOD…`. That `AF` is GATK
**cohort** allele frequency over 4 samples. There are zero occurrences of `AF_TGP`,
`AF_EXAC`, `AF_ESP` or `gnomAD_AF`. Anything reading `AF` as population frequency gets a
confidently wrong answer. `ve-frequency` must detect and name this, not silently use it.

### 3. Build mismatch fails silently in existing tooling

`variant-annotation` hardcodes `VEP_BASE_URL = "https://rest.ensembl.org"` (the GRCh38
host); its `--assembly GRCh37` flag only adds a query parameter that host ignores — GRCh37
needs `grch37.rest.ensembl.org`, a different server. `vcf-annotator` pins
`dataset: gnomad_r4`, which is GRCh38-only. Our data is b37. **Match the build explicitly
or declare that you could not.**

### 4. Consequence belongs to variant × transcript, not variant

31 of 68 records also carry a non-HIGH annotation; 18 of 68 touch more than one gene
(median 2 annotations per record, max 16). The HIGH label is SnpEff's worst-consequence
pick. `ve-router` must emit the routing decision *as data* — which transcript it chose and
what it discarded.

Note the `EFF` field is **legacy SnpEff**, not modern `ANN`. `rare-high-impact-variants`
cannot parse it; `ve-router` owns that parsing.

---

## Skills to build

| Skill | Input → Output | Notes |
|---|---|---|
| **`ve-segregation`** | VCF + pedigree → Contract A `.segregation` | Proband carries AND exactly one parent carries. Unphased — a transmission-consistency label, **not molecular phase**. Must derive roles from data. |
| **`ve-frequency`** | Contract A → `.freq` populated | A **frequency provenance gate**, not another annotator. Is there a usable population-frequency layer at all, in what build, and is any `AF`-looking field actually cohort frequency? Refuses to emit rarity when provenance can't be established. |
| **`ve-router`** | Contract A → `.class` set | Parses legacy `EFF` (and VEP if present), picks one consequence per variant × transcript, emits the decision and the discarded alternatives. |
| **`ve-lof`** | Contract A → Contract B | LOFTEE-style confidence: last exon, NMD escape, low-confidence flags. Our own implementation. Ground truth for validation is free — gnomAD publishes its LOFTEE HC/LC calls. |
| **`ve-missense`** | Contract A → Contract B | AlphaMissense. Prefer the **precomputed score table** (hg19 + hg38) — a lookup, not inference. |
| **`ve-splice`** | Contract A → Contract B | Our own ref-vs-alt splice delta. Prefer **precomputed SpliceAI scores** over running inference; check licensing before committing. |
| **`ve-regulatory`** | Contract A → Contract B | ChromBPNet ref-vs-alt delta, cell-type specified. Heaviest branch — build last. |
| **`ve-merge`** | Contract B[] + Contract A → report | Ranked list + abstention list, with segregation as ranking evidence. Emits evidence, **not** ACMG codes. |

**Model-branch strategy:** for `ve-missense`, `ve-splice` and `ve-regulatory`, prefer a
precomputed score table over running inference. It is faster, needs no GPU, and keeps the
skill self-contained. Ship a small pre-extracted slice as demo data so `--demo` works
offline.

---

## These MUST be ClawBio-format skills

This is not optional and it is not a stylistic preference. A skill that does not follow
this layout **cannot be discovered or run** by `clawbio_list_skills` /
`clawbio_describe_skill` / `clawbio_run_skill`, which is how the hosted BioNeMo agent
reaches everything. If the format is wrong, the demo does not happen.

### Directory layout

```
skills/ve-<name>/
├── SKILL.md          # YAML frontmatter + prose contract (REQUIRED)
├── ve_<name>.py      # CLI: --input, --output, --demo (REQUIRED)
├── api.py            # importable run() entrypoint
├── demo_input.txt    # TINY — library median is 892 bytes, keep under ~10 KB
└── tests/
```

### SKILL.md frontmatter — copy this exactly

**Use the nested `metadata:` form.** The published spec page at
`docs.clawbio.ai/reference/skillmd-spec/` shows a *flat* top-level form — ignore it.
**97 of 97 skills in the repo use the nested form**, and the catalog generator expects it.

```yaml
---
name: ve-segregation
description: >-
  One sentence, active voice, saying exactly what the skill does and nothing else.
license: MIT
metadata:
  version: 0.1.0
  author: <your name>
  domain: genomics
  inputs:
    - name: input_file
      type: file
      format:
        - vcf
        - tsv
      description: Primary input data file
      required: true
  outputs:
    - name: report
      type: file
      format: md
      description: Analysis report
    - name: result
      type: file
      format: json
      description: Machine-readable results
  dependencies:
    python: ">=3.11"
  tags:
    - segregation
    - pedigree
    - inheritance
  demo_data:
    - path: demo_input.txt
      description: Synthetic test data
  endpoints:
    cli: python skills/ve-segregation/ve_segregation.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    homepage: https://github.com/burgshrimps/ve-skills
    os:
      - macos
      - linux
    install: []
    trigger_keywords:
      - segregation
      - parent of origin
      - pedigree
---
```

### SKILL.md body — required sections

Follow the house pattern (see any existing skill, e.g. `rare-high-impact-variants`):

- **`# Title`** then one paragraph: *"You are **X**, a specialised ClawBio agent for …"*
- **`## Trigger`** — "Fire this skill when the user says any of: …" and "Do NOT fire when: …".
  Be loud and literal; models skip subdued descriptions. List exact phrases and synonyms.
- **`## Why This Exists`** — without it / with it / why ClawBio
- **`## Input Formats`** — table of format, extension, required fields
- **`## Workflow`** — numbered steps the agent follows
- **`## Domain Decisions`** — **every threshold and rule, with its source.** This is the
  core of the spec: classification cutoffs, the exact segregation rule, AF thresholds.
- **`## Safety Rules`** — what the skill must never do. Ours all inherit the interpretation
  boundary below: never say rare / pathogenic / diagnostic / de novo / compound het.
- **`## Agent Boundary`** — explicit In Scope / Out of Scope lists
- **`## CLI Reference`** — standard usage plus `--demo`
- **`## Example Output`** — a realistic `report.md` excerpt
- **`## Output Structure`** — the tree below

### Output directory convention

```
<output>/
├── report.md                    # primary markdown report
├── result.json                  # machine-readable results
├── tables/results.csv           # tabular data
└── reproducibility/commands.sh  # exact commands to reproduce
```

### Two hard requirements

1. **`--demo` must work with no network and no user files.** It is the only mode the hosted
   agent is guaranteed to be able to run; `input_path` is often refused for local files.
2. **Every `report.md` ends with the disclaimer:**
   *ClawBio is a research and educational tool. It is not a medical device and does not
   provide clinical diagnoses. Consult a healthcare professional before making any medical
   decisions.*

---

## Build order

The hackathon build window is 2 h 50 min. Eight skills will not happen. Ordered by value
per minute:

1. **`ve-segregation`** — the hour-one deliverable; pure genotype logic, no models, no network
2. **`ve-router`** — nothing downstream works without it; pure parsing
3. **`ve-merge`** — turns the above into the actual deliverable
4. **`ve-lof`** — fires on 56 of 68 records
5. **`ve-frequency`** — the blind spot the brief explicitly names
6. **`ve-splice`** — 12 records
7. **`ve-missense`** / **`ve-regulatory`** — outside demo data only

**Items 1–3 are a complete, honest, end-to-end demo.** The brief says plainly that one
thing working end to end beats four things half-built. Cut from the bottom.

---

## Interpretation boundary (non-negotiable, and it is judged)

Never call anything **rare, pathogenic, diagnostic, de novo, or compound heterozygous**.

`PARENT_OF_ORIGIN_UNPHASED` is a teaching label, not molecular phase. There is no
phenotype, no HPO terms, no valid population-frequency layer, and `EFF` is historical
annotation rather than current clinical evidence.

Every claim must trace to something actually run or read. "I cannot determine this"
scores higher than a confident guess.

Data: derived subset of the Corpasome by Manuel Corpas,
[DOI 10.6084/m9.figshare.693052](https://figshare.com/articles/dataset/Corpasome/693052),
CC BY 4.0. Keep the attribution.

---

## Work split

See [`docs/work-packets.md`](docs/work-packets.md) — three self-contained briefs, one per
person, each pasteable straight into an agentic coding session.
