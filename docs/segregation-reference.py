#!/usr/bin/env python3
"""Hour-one deliverable seed: reproduce the 30 paternal / 38 maternal split
from genotypes alone, with the filter trail visible.

Sample roles are DERIVED from the data (TSV column headers, verified 68/68
genotype-identical to the VCF sample columns), not from the challenge prose.
"""
import gzip, csv, collections, sys

VCF = 'challenge1-b37-segregation.vcf.gz'
TSV = 'challenge1-b37-segregation.tsv'

ROLES = {'ISDBM322015': 'son', 'ISDBM322016': 'father',
         'ISDBM322017': 'sister', 'ISDBM322018': 'mother'}
SON, FATHER, MOTHER = 'ISDBM322015', 'ISDBM322016', 'ISDBM322018'

def carries(gt):
    """Genotype contains the ALT allele. Unphased; no phase is claimed."""
    return gt not in ('0/0', './.', '0|0')

def load_vcf(path):
    out = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                samples = f[9:]
                continue
            out.append(((f[0], f[1]),
                        {s: f[9 + i].split(':')[0] for i, s in enumerate(samples)}))
    return out

def segregate(recs):
    """Son carries AND exactly one parent carries -> label that parent.
    This is a transmission-consistency label, NOT molecular phase."""
    calls = {}
    for key, g in recs:
        if not carries(g[SON]):
            calls[key] = 'excluded:son-does-not-carry'
            continue
        f, m = carries(g[FATHER]), carries(g[MOTHER])
        if f and not m:
            calls[key] = 'paternal'
        elif m and not f:
            calls[key] = 'maternal'
        else:
            calls[key] = 'ambiguous:both-or-neither-parent-carries'
    return calls

if __name__ == '__main__':
    recs = load_vcf(VCF)
    calls = segregate(recs)
    print('Roles used (derived from data):')
    for sid, role in ROLES.items():
        print(f'  {sid} = {role}')
    print(f'\nRecords: {len(recs)}')
    print('Computed:', dict(collections.Counter(calls.values())))

    with open(TSV) as fh:
        labels = {(r['CHROM'], r['POS']): r['PARENT_OF_ORIGIN_UNPHASED']
                  for r in csv.DictReader(fh, delimiter='\t')}
    agree = sum(1 for k, v in calls.items() if v == labels.get(k))
    print(f'Agreement with PARENT_OF_ORIGIN_UNPHASED: {agree}/{len(recs)}')
    sys.exit(0 if agree == len(recs) else 1)
