---
name: ve-segregation
description: >-
  Determine parent-side segregation for filtered pedigree VCF records and flag
  genes with historical HIGH-effect annotations observed from both parental sides.
license: MIT
metadata:
  version: 0.1.0
  author: ClawBio Berlin 2026 team
  domain: genomics
  inputs:
    - name: input
      type: file
      format:
        - vcf
        - vcf.gz
      description: Filtered pedigree VCF containing GT, DP, GQ, and SnpEff EFF annotations
      required: true
  outputs:
    - name: result
      type: file
      format: json
      description: Machine-readable segregation summary
    - name: variant_segregation
      type: file
      format: tsv
      description: High-effect variant/transcript rows with parental-side calls
    - name: gene_segregation
      type: file
      format: tsv
      description: Gene-level segregation summary sorted with biparental flags first
  dependencies:
    python: ">=3.10"
  tags:
    - variant-effect
    - segregation
    - pedigree
    - rare-disease
    - abstention
---

# ve-segregation

`ve-segregation` is the first skill in the local `ve-*` variant-effect layer.
It reads a filtered pedigree VCF and determines whether each carried variant is
observed on the paternal side, maternal side, or is ambiguous under simple trio
carrier logic.

This skill does **not** infer rarity, pathogenicity, molecular phase, compound
heterozygosity, or diagnosis. It only produces segregation evidence and a lure
flag for genes with historical HIGH-effect annotations from both parental sides.

## CLI

```bash
python ve_segregation.py \
  --input .inputs/challenge1-b37-segregation.vcf.gz \
  --output .outputs/ve_segregation_run \
  --case-sample ISDBM322015 \
  --father-sample ISDBM322016 \
  --mother-sample ISDBM322018 \
  --sibling-sample ISDBM322017
```

## Outputs

The skill writes exactly three files to the output directory:

- `result.json` — run summary, sample mapping, counts, output paths, and
  interpretation boundaries.
- `variant_segregation.tsv` — one row per historical HIGH-effect annotation.
- `gene_segregation.tsv` — one row per gene, sorted with biparental high-effect
  lure genes first.

## Interpretation Boundary

- Segregation side is not molecular phase.
- A biparental high-effect flag is not a compound-heterozygous call.
- Historical SnpEff `EFF` is not modern effect validation.
- No rarity, pathogenicity, or diagnosis is inferred.
