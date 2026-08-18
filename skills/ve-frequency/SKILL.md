---
name: ve-frequency
description: >-
  A frequency provenance gate for the ve-skills pipeline. Establishes whether a usable
  population-frequency layer exists at all, in which reference build, and whether any
  AF-looking INFO field is actually cohort frequency in disguise — then refuses to emit
  a rarity call when provenance cannot be established.
license: MIT
metadata:
  version: 0.1.0
  author: ve-skills team
  domain: genomics
  tags:
    - allele-frequency
    - provenance
    - build-awareness
    - liftover
    - variant-filtering
  inputs:
    - name: input_file
      type: file
      format:
        - json
        - vcf
      description: Contract A records from ve-segregation, or a VCF to bootstrap from
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Class counts plus an explicit tally of why the gate refused
    - name: result
      type: file
      format:
        - json
      description: Contract A records with .freq populated
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_input.txt
      description: b37 VCF with GATK cohort AF only, plus two records carrying a real population AF
---

# ve-frequency

Populates `.freq` on Contract A. Sits between `ve-segregation` and `ve-router`.

## In plain English

**What it does:** for each variant it answers "how common is this in the general
population?" — or, far more often on this data, it answers **"I cannot tell you, and
here is exactly why."**

**The problem it solves.** Rarity is the first filter anyone applies to a variant list:
something carried by 30% of people is not what is making one patient ill. The obvious
move is to read the `AF` field sitting right there in the VCF.

That field is a trap. In a GATK-called VCF, `AF` is the allele frequency **within the
samples in the file** — here, one family of four. A variant seen in two of eight alleles
gets `AF=0.250`. Read as population frequency, that says "25% of humans carry this,
common, discard it". It actually says "half this family carries it", which for an
inherited disease is the opposite of a reason to discard.

Two more ways the same question goes quietly wrong:

- **The coordinates are in the wrong dialect.** Genome positions were renumbered between
  reference versions. This project's data is b37; gnomAD v4 is GRCh38-only. Ask in the
  wrong numbering and you get a confident "not found" that reads as "ultra-rare".
- **Nobody sequenced that spot.** Coverage varies across the genome. "Not seen in the
  reference database" at a poorly-covered site means "not looked at", not "not there".

None of the three produce an error. They produce a plausible number that is wrong, and a
shortlist built on it.

**What this skill does about it.** It refuses. Before emitting any rarity call it checks
that a genuine population-frequency source exists, that the reference build is known, and
that the field it is about to use is not cohort frequency wearing a population label. If
any check fails the variant gets `NO_DATA` and a warning naming the failure — never a
number. On the challenge data that means most variants come back `NO_DATA`, which is the
correct answer, not a gap.

**What it does not do:** it does not say whether a variant is harmful. Rare is not the
same as causal.

## Trigger

**Fire when:**
- populating `.freq` on Contract A records
- "is this variant rare", "what's the population frequency", "filter by frequency"
- "is INFO/AF usable", "is that cohort AF", "which build is this"
- any step that would filter variants by rarity

**Do NOT fire when:**
- routing consequences → `ve-router`
- scoring impact → `ve-lof` / `ve-missense` / `ve-splice` / `ve-regulatory`
- ranking or merging evidence → `ve-merge`
- inheritance patterns → `ve-segregation`

## Why This Exists

`AF` in the challenge VCF is GATK cohort frequency over four samples; the INFO keys are
`AC, AF, AN, MLEAC, MLEAF, DP, EFF, VQSLOD` with **zero** occurrences of `AF_TGP`,
`AF_EXAC`, `AF_ESP` or `gnomAD_AF`. Anything reading `AF` as population frequency is
confidently wrong on every record. Existing tooling does not catch this: `vcf-annotator`
pins `dataset: gnomad_r4` (GRCh38-only) against b37 data, and `variant-annotation`'s
`--assembly GRCh37` flag only adds a query parameter to the GRCh38 host, which ignores it.

## Core Capabilities

1. **Cohort-AF detection**: names `INFO/AF` as cohort frequency and refuses it.
2. **Build establishment**: from `##reference` and `##contig` lengths, or `None`.
3. **Local liftover**: b37→b38 from a vendored Ensembl chain, minus-strand aware.
4. **Coverage-aware absence**: absence at an uncovered site is `NO_DATA`, not `RARE`.
5. **Refusal accounting**: the report tallies *why* the gate refused, per reason.

## Scope

**One skill, one task.** It establishes frequency provenance and emits `.freq`. It does
not annotate consequence, score impact, or rank variants.

## Input Formats

| Format | Source | Notes |
|---|---|---|
| Contract A JSON | `ve-segregation` output | `{records: [...], header: [...], sample_count: N}` or a bare array |
| VCF (plain or `.gz`) | challenge pack | Bootstraps Contract A directly when run standalone |

## Workflow

1. **Load** Contract A records, or bootstrap them from a VCF.
2. **Establish the build** from the header. Unknown build → every record `NO_DATA`.
3. **Assess AF provenance** per record: a population key is usable; `AF`/`MLEAF` are not.
4. **Lift** b37 coordinates to b38 only if a reference lookup will be consulted.
5. **Classify** into `RARE` / `COMMON` / `NO_DATA`.
6. **Emit** `.freq`, leaving all upstream fields untouched.
7. **Report** counts and a per-reason tally of refusals.

## Output — Contract A `.freq`

```json
"freq": {
  "af": null,
  "source": null,
  "build": "GRCh37",
  "class": "NO_DATA",
  "provenance_warning": "INFO/AF present but is GATK cohort AF over 4 samples, not population AF"
}
```

`class` ∈ `RARE` | `COMMON` | `NO_DATA`

## Example Output

```
ve-frequency: 6 records, build=GRCh37
  RARE     1
  COMMON   1
  NO_DATA  4
  refused (4): INFO/AF present but is GATK cohort AF over 4 samples, not population AF
```

| Variant | Class | AF | Source |
|---|---|---|---|
| `1:11906068:A:G` | NO_DATA | — | — |
| `1:55505647:G:T` | NO_DATA | — | — |
| `1:11906069:C:T` | RARE | 0.0004 | gnomAD_AF |
| `1:11906070:G:A` | COMMON | 0.31 | gnomAD_AF |

## CLI Reference

```bash
python skills/ve-frequency/ve_frequency.py --demo --output /tmp/ve_freq_demo
python skills/ve-frequency/ve_frequency.py --input contract_a.json --output out/
python skills/ve-frequency/ve_frequency.py --input challenge1-b37.vcf.gz --output out/
python skills/ve-frequency/ve_frequency.py --input in.vcf --assembly GRCh37 --output out/
```

## Gotchas

1. **The model will want to use `INFO/AF`. It is cohort frequency.** Over four samples,
   `AF=0.250` means two of eight alleles in *this family*, not 25% of humans. Reading it
   as population AF inverts the filter: family-segregating variants look common and get
   discarded. `AF` and `MLEAF` are on the refusal list by name.

2. **The model will want to treat `NO_DATA` as a failure. It is the finding.** On the
   challenge pack most records are `NO_DATA` because there genuinely is no population
   layer in the file. A pipeline that returned `RARE` for those would be fabricating.
   Downstream must carry the abstention, not drop it.

3. **The model will want to default the build. Refuse instead.** b37 coordinates queried
   against a GRCh38 dataset return "not found", which is indistinguishable from ultra-rare.
   Build is read from `##reference` and `##contig` lengths; when they conflict or are
   absent, `detect_build` returns `None` and every record becomes `NO_DATA`.

4. **The model will want to treat "absent from gnomAD" as rare. Check coverage first.**
   Absence at a site where 32% of samples reached 20x carries no information. Absence at a
   well-covered site does, and is bounded by the rule of three (~3/AN), not reported as zero.

5. **The model will want to skip allele reverse-complement on liftover. 32.5% of chain
   blocks are minus-strand.** A minus-strand lift that keeps the original alleles produces
   a valid-looking variant ID that does not exist.

6. **Never vendor a gnomAD frequency table.** gnomAD is ODbL (share-alike); this repo is
   MIT. Runtime caches are fine; a checked-in extract is a licence conflict. See
   `data/PROVENANCE.md`.

## Dependencies

Python ≥3.10, standard library only. The liftover chain is vendored (285 KB); no
`pyliftover`, `CrossMap`, `bcftools` or reference FASTA required.

## Safety

Reports population frequency provenance only. Makes no clinical claim; `NO_DATA` means
the gate could not establish provenance, never that a variant is absent or benign.
Runs fully offline.

> ClawBio is a research and educational tool. It is not a medical device and does not
> provide clinical diagnoses. Consult a healthcare professional before making any
> medical decisions.

## Agent Boundary

The agent dispatches and explains the refusals. The skill executes. The agent does not
supply `--assembly` to work around a failed build detection, does not reinterpret
`NO_DATA` as rare, and does not read `INFO/AF` itself.

## Data Attribution

Liftover chain from Ensembl (`GRCh37_to_GRCh38.chain.gz`), redistributed under EMBL-EBI
terms, which place no additional restrictions on redistribution. Provenance, SHA-256 and
validation in `data/PROVENANCE.md`.
