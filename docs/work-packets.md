# Work packets

Three packets, one per person. Each is self-contained — **paste your whole packet into your
agentic coding session** (Claude Code, Cursor, whatever you're using) and it has everything
it needs.

Read [`../README.md`](../README.md) first for the pipeline, contracts and traps. Every
packet assumes it.

| Person | Skills | Why this grouping |
|---|---|---|
| **A** | `ve-router`, `ve-merge` | The spine. Pure parsing and assembly — no models, no network. Everything else plugs into it. |
| **B** | `ve-segregation`, `ve-frequency` | VCF-level evidence. Both fire on the challenge data, both have verified ground truth to check against. |
| **C** | `ve-lof`, `ve-splice` | The two branches that actually fire (56 and 12 records). |
| *stretch* | `ve-missense`, `ve-regulatory` | Whoever finishes first. Outside demo data only. |

**`ve-segregation` is the hour-one deliverable — it must not slip.** If Person B is
blocked, someone else picks it up.

**Agree Contracts A and B before anyone writes code.** They are in the README. Nobody
changes them alone.

---

## Common preamble — paste this at the top of every packet

```
We are building `ve-skills`, a self-contained variant-effect layer for the ClawBio
Berlin 2026 hackathon, Challenge 1 ("End the diagnostic odyssey"). Repo:
https://github.com/burgshrimps/ve-skills

Read README.md in the repo root first. It has the pipeline diagram, the shared JSON
contracts, the verified facts about the data, and four documented traps. Everything
below assumes it.

HARD CONSTRAINTS
1. These must be ClawBio-format skills. SKILL.md with NESTED `metadata:` frontmatter
   (97 of 97 repo skills use nested; the published spec page shows a flat form that
   nothing uses — ignore it). Directory layout, body sections and output structure are
   specified in README.md under "These MUST be ClawBio-format skills". A skill that
   doesn't follow it cannot be run by the hosted agent, so the format is not optional.
2. `--demo` must work with NO network and NO user files. Ship a tiny demo input
   (under ~10 KB; library median is 892 bytes).
3. We implement everything ourselves. Do NOT call gi-splice, gi-chromatin,
   clinical-variant-reporter or any other existing ClawBio skill at runtime.
4. Never output the words rare, pathogenic, diagnostic, de novo, or compound
   heterozygous as a claim about a variant. Absence of evidence is not evidence.
5. Every report.md ends with: "ClawBio is a research and educational tool. It is not a
   medical device and does not provide clinical diagnoses. Consult a healthcare
   professional before making any medical decisions."
6. Emit evidence, not verdicts. Where you cannot determine something, say so explicitly
   and put it in the abstention output rather than guessing.

DATA (already downloaded, b37/GRCh37):
  challenge1-b37-segregation.vcf.gz   68 records
  challenge1-b37-segregation.vcf.gz.tbi
  challenge1-b37-segregation.tsv      readable genotype table

VERIFIED FACTS — trust these, they were checked against the file:
  68 records, all PASS, all autosomal, all biallelic
  35 SNVs + 33 indels (max 4 bp). ZERO CNVs/SVs.
  56 protein-truncating (27 frameshift, 21 stop-gained, 5 start-lost, 3 stop-lost)
  12 splice-site (7 donor, 4 acceptor, 1 both)
  0 missense as a selected consequence
  30 paternal / 38 maternal segregation labels
  Sample roles: ISDBM322015=son, 322016=father, 322017=SISTER, 322018=MOTHER
    (the challenge brief says mother/sister in the other order — the brief is WRONG,
     derive roles from the data, never hardcode from prose)
  INFO/AF in this VCF is GATK COHORT frequency over 4 samples, NOT population frequency.
  The consequence field is legacy SnpEff `EFF`, not modern `ANN`.
```

---

# Packet A — `ve-router` + `ve-merge`

> Paste the common preamble above, then this.

```
You are building two skills: ve-router and ve-merge. These are the spine of the
pipeline — pure parsing and assembly, no models and no network. Build ve-router first;
nothing downstream works without it.

## ve-router

INPUT:  Contract A records (or a raw VCF — accept both; if given a VCF, parse it and
        emit Contract A with `class` set).
OUTPUT: Contract A with `class` populated, plus the routing decision as DATA.

WHAT IT DOES
Parse the legacy SnpEff `EFF` INFO field and assign exactly one routing class per
variant. Classes: protein_truncating | missense | splice | non_coding | unroutable.

EFF format is:
  Effect(Impact|FunctionalClass|Codon|AA|AA_len|Gene|BioType|Coding|Transcript|Rank|GT)
e.g.
  STOP_LOST(HIGH|MISSENSE|Tga/Cga|*152R|151|NPPA||CODING|NM_006172.3|3|1)

Map effects to classes:
  protein_truncating <- FRAME_SHIFT, STOP_GAINED, START_LOST, STOP_LOST
  splice             <- SPLICE_SITE_DONOR, SPLICE_SITE_ACCEPTOR
  missense           <- NON_SYNONYMOUS_CODING
  non_coding         <- everything in a non-coding context
  unroutable         <- anything you cannot confidently classify

THE HARD PART, AND THE POINT OF THIS SKILL
Consequence is a property of variant x TRANSCRIPT, not of the variant. In this dataset:
  - median 2 EFF annotations per record, max 16
  - 31 of 68 records also carry a NON-HIGH annotation
  - 18 of 68 records touch more than one gene symbol
So you must pick one, and you must SHOW THE CHOICE. Emit, per variant:
  - the chosen consequence, gene and transcript
  - the rule used to choose (e.g. "highest impact tier, then first transcript")
  - the discarded alternatives, in full
Never silently collapse to one consequence. A reviewer must be able to see what was
dropped. Put the discarded set in the record and summarise the count in report.md.

Also: if a record has NO parseable EFF, class = unroutable with a reason. Do not guess.

VALIDATE AGAINST: routing the 68 challenge records must give
  56 protein_truncating, 12 splice, 0 missense, 0 non_coding, 0 unroutable.
Write that as a test.

## ve-merge

INPUT:  Contract B outputs from all branches + the Contract A records
OUTPUT: report.md with (1) a ranked list and (2) an abstention list

RANKING: combine branch score/confidence with segregation evidence from
Contract A `.segregation`. Keep the ranking rule simple, explicit and written down in
SKILL.md under "Domain Decisions". Show the inputs to each rank, not just the rank.

THE ABSTENTION LIST IS THE MOST IMPORTANT OUTPUT. It is the challenge's stretch goal and
the thing it says it scores highest. Build it BY CONSTRUCTION, not by hand:
  - every branch result with in_domain=false, carrying its abstain_reason
  - every variant with class=unroutable
  - every variant with freq.class=NO_DATA (absence of a frequency is NOT rarity)
  - anything the pipeline never reached
Each entry needs a REASON. "No gnomAD entry" is not a reason; "no build-matched
population-frequency source was available for GRCh37, so rarity could not be
established" is.

Add a standing section to report.md listing what this data cannot support at all:
no phenotype, no HPO terms, no valid population-frequency layer, EFF is historical
annotation rather than current clinical evidence, and the parent-of-origin labels are
unphased transmission labels rather than molecular phase.

Ship both skills with tiny demo inputs and tests. Start with ve-router.
```

---

# Packet B — `ve-segregation` + `ve-frequency`

> Paste the common preamble above, then this.

```
You are building two skills: ve-segregation and ve-frequency. Build ve-segregation
FIRST — it is the hackathon's hour-one deliverable and the whole team depends on it
landing.

## ve-segregation

INPUT:  multi-sample VCF with pedigree genotypes
OUTPUT: Contract A with `.segregation` populated

THE RULE
A variant is assigned a parent of origin when:
  the proband carries the ALT allele, AND exactly one parent carries it.
"Carries" means the genotype contains ALT — i.e. NOT in {0/0, ./., 0|0}.
If both parents carry, or neither does, the pattern is `ambiguous` with a reason.
If the proband does not carry, the pattern is `excluded` with a reason.

CRITICAL — DERIVE THE PEDIGREE, DO NOT HARDCODE IT
The challenge brief states the sample order as "son, father, mother, sister". THIS IS
WRONG: mother and sister are swapped. The correct mapping is
  ISDBM322015=son, ISDBM322016=father, ISDBM322017=sister, ISDBM322018=mother
Get roles from the TSV column headers (SON_/FATHER_/SISTER_/MOTHER_), or accept an
explicit pedigree file (Contract C in the README). Never infer them from column order
and never hardcode them from the prose.

Proof this matters, and worth putting in report.md as a validation section:
  correct mapping  -> 30 paternal / 38 maternal, 68/68 agreement with the
                      PARENT_OF_ORIGIN_UNPHASED labels
  brief's mapping  -> 11 / 25 / 32 ambiguous, only 36/68 agreement

SAFETY — THIS IS JUDGED
These are UNPHASED transmission-consistency labels, NOT molecular phase. Say so in the
report, in SKILL.md Safety Rules, and in the JSON output. Never call anything de novo
or compound heterozygous — this data cannot support either claim.

MAKE THE FILTER TRAIL VISIBLE. The brief asks for the logic, not the two numbers. For
each variant show: proband genotype, each parent genotype, which rule fired, and the
resulting label. Emit as tables/results.csv as well as report.md.

VALIDATE AGAINST: 30 paternal / 38 maternal, 68/68 label agreement. Write it as a test.

## ve-frequency

This is NOT another annotator. ClawBio already has vcf-annotator, variant-annotation and
rare-high-impact-variants. It is a FREQUENCY PROVENANCE GATE. Its job is to answer three
questions and refuse to emit a rarity call when it cannot:
  1. Is there a usable population-frequency layer in this file AT ALL?
  2. What genome build are the coordinates in, and does the frequency source match it?
  3. Is any AF-looking field actually something else wearing that name?

TRAP 1 — COHORT AF MASQUERADING AS POPULATION AF
This VCF's INFO keys are: AC, AF, AN, MLEAC, MLEAF, DP, EFF, VQSLOD, ...
That `AF` is GATK cohort allele frequency across FOUR SAMPLES. It is not population
frequency. There are ZERO occurrences of AF_TGP, AF_EXAC, AF_ESP or gnomAD_AF.
Your skill must DETECT this and emit a provenance_warning naming it. Do not use it.
Getting this right is the single most demoable thing in the project — the challenge
brief cites an audit where 16 of 27 false-actionable calls came from exactly this class
of mistake.

TRAP 2 — SILENT BUILD MISMATCH
Our data is GRCh37/b37. gnomAD r4 is GRCh38-only; gnomAD v2.1.1 is the GRCh37 release.
Ensembl's main REST host (rest.ensembl.org) serves GRCh38 — GRCh37 needs
grch37.rest.ensembl.org, a DIFFERENT HOST. Passing an assembly parameter to the GRCh38
host does nothing; it returns GRCh38 answers with no error. (Two existing ClawBio skills
have exactly this bug — see README trap 3.) Either match the build explicitly or declare
that you could not.

OUTPUT: Contract A `.freq` = {af, source, build, class, provenance_warning}
  class is RARE | COMMON | NO_DATA
NO_DATA is a first-class result, not a null. On the challenge data, ALL 68 records
should come out NO_DATA with a provenance warning — and that is the CORRECT answer.
Write that as a test.

Ship both skills with tiny demo inputs and tests. Start with ve-segregation.
```

---

# Packet C — `ve-lof` + `ve-splice`

> Paste the common preamble above, then this.

```
You are building two skills: ve-lof and ve-splice. These are the two branches that
actually fire on the challenge data — 56 records and 12 records respectively. Build
ve-lof first.

Both take Contract A and return Contract B EXACTLY:
{ skill, variant_key, score, direction, confidence, in_domain, abstain_reason, evidence }

The in_domain / abstain_reason fields are not boilerplate. When handed something outside
your domain, return in_domain=false with a specific, human-readable reason. ve-merge
turns these into the abstention list, which is the challenge's highest-scoring output.

## ve-lof

INPUT:  Contract A records with class=protein_truncating
OUTPUT: Contract B

Assess whether a predicted loss-of-function call is TRUSTWORTHY. A stop-gained in the
last exon is not the same as one in exon 2, and pipelines routinely treat them alike.
Implement LOFTEE-style confidence ourselves (do not shell out to LOFTEE). Flags to
compute, each with its reasoning recorded in `evidence`:
  - variant in the last exon, or the last ~50 bp of the penultimate exon
    -> escapes nonsense-mediated decay -> lower confidence
  - position within the transcript (a truncation at 95% of the CDS removes little)
  - for frameshifts: predicted premature-termination-codon position
  - single-exon genes (no NMD)
  - low-confidence annotations, e.g. non-canonical transcript only

Document every threshold in SKILL.md under "Domain Decisions" with its rationale.
Score = confidence that this is a genuine, consequential LoF. NOT a pathogenicity call —
never use that word.

WHERE THE SCALE COMES FROM: 27 of your 56 records are FRAMESHIFT INDELS, not SNVs.
Handle indels properly; they are the largest class here, not an edge case.

FREE GROUND TRUTH: gnomAD publishes its own LOFTEE HC/LC calls. Pull a small set of
known-HC and known-LC variants as demo data and report CONCORDANCE with them. Showing
agreement (and where you disagree) is far stronger than asserting your score is right.

## ve-splice

INPUT:  Contract A records with class=splice
OUTPUT: Contract B

Score the CHANGE in splice-site strength between ref and alt. This is the distinction
that makes us non-overlapping: ClawBio's gi-splice detects sites in a sequence and never
compares alleles. We do the delta, and we implement it ourselves — do not call gi-splice.

Prefer a PRECOMPUTED SpliceAI score table over running inference: no GPU, no weights, no
reference FASTA, seconds instead of an afternoon. CHECK LICENSING AND AVAILABILITY
BEFORE YOU COMMIT — if it is not cleanly usable, say so immediately and fall back to
running the spliceai package locally, or report the branch as unavailable rather than
faking it. Ship a tiny pre-extracted slice as demo data so --demo works offline.

IMPORTANT DOMAIN NOTE: your 12 records are splice-site variants at canonical intronic
donor/acceptor positions. They are NOT coding-sequence variants — they do not change a
codon. Their HIGH impact rests on a PREDICTION about splicing that this data does not
test. Say that explicitly in the report. It is a legitimate abstention-list entry.

Note also that 7 of the 12 are indels, not SNVs. Precomputed SpliceAI tables generally
cover SNVs only — if that is the case, the 7 indels must return in_domain=false with a
reason, and that is a correct, honest result, not a failure.

Ship both skills with tiny demo inputs and tests. Start with ve-lof.
```

---

## Definition of done (any skill)

- [ ] `SKILL.md` with nested `metadata:` frontmatter, all required body sections
- [ ] `python skills/ve-<name>/ve_<name>.py --demo --output /tmp/x` works with **no network**
- [ ] `report.md`, `result.json`, `tables/`, `reproducibility/` all written
- [ ] Disclaimer line present in `report.md`
- [ ] Contract A / Contract B shapes match the README exactly
- [ ] At least one test asserting a verified number from the README
- [ ] Every threshold documented under "Domain Decisions" with its source
- [ ] No banned words (rare / pathogenic / diagnostic / de novo / compound heterozygous)
      used as a claim about a variant
