# Data provenance — ve-frequency

## `GRCh37_to_GRCh38.chain.gz`

| | |
|---|---|
| Source | <https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz> |
| Retrieved | 2026-08-18 |
| Size | 285,250 bytes |
| SHA-256 | `351de3cd4a01d9fcffd38881981767b697090d2eba876740891b96d5c546b100` |
| Contig naming | Ensembl style, **no `chr` prefix** (`1`, `17`, `X`) — matches gnomAD variant IDs |
| Licence | EMBL-EBI terms of use: *"places no additional restrictions on the use or redistribution of the data"*. Attribution expected. |

### Why Ensembl and not UCSC

Both chains were validated against the same four ground-truth loci and both got
4/4 exact. Ensembl's was chosen for two reasons:

1. **Redistribution is explicitly permitted.** EMBL-EBI states it places no
   additional restrictions on redistribution. UCSC's terms are less clear for
   vendoring into a third-party MIT repository.
2. **Contig naming matches the target.** Ensembl chains use `1`, not `chr1`,
   which is the same convention as gnomAD variant IDs (`1-55039974-G-T`), so
   one less normalisation step sits between the chain and the query.

### Validation

Verified against two independent sources — Ensembl REST `/map/human/GRCh37/…/GRCh38`
and gnomAD's own `liftover` GraphQL query:

| GRCh37 | GRCh38 | Agreement |
|---|---|---|
| 17:41246481 | 17:43094464 | Ensembl + chain |
| 1:55505647 | 1:55039974 | Ensembl + gnomAD + chain |
| 13:32906729 | 13:32332592 | Ensembl + chain |
| 7:117199644 (indel `ATCT>A`) | 7:117559590 | gnomAD + chain |

`skills/ve-frequency/tests/` re-runs these four as a regression
test. If a chain refresh changes any of them, that is a red flag, not a rebase.

### What is deliberately NOT here

**No gnomAD frequency data is vendored in this directory, and none may be
added.** gnomAD is released under the [ODC Open Database License
(ODbL)](https://opendatacommons.org/licenses/odbl/), which is share-alike:
a derived database must itself be released under ODbL. ClawBio is MIT. A
checked-in `gnomad_af.parquet` — or any extract of gnomAD allele frequencies —
would put the two licences in conflict.

Runtime caches written into the user's own output or cache directory are fine;
those are not redistribution. `tests/test_variant_frequency_resolver.py::test_demo_does_not_vendor_gnomad_frequency_tables`
enforces this.
