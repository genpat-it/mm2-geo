#!/usr/bin/env python3
"""Discordant primary alignments (identical placement, different CIGAR) between stock (geo0) and heuristic
(geo1) PAF files: score loss, severity, truncation and enrichment in GIAB stratifications (chr22)."""
import sys, gzip, numpy as np
d, sdir, conf = sys.argv[1], sys.argv[2], sys.argv[3]
N = 51_000_000
def cov(f, chrom):
    a = np.zeros(N, bool); op = gzip.open if f.endswith('gz') else open
    for l in op(f, 'rt'):
        t = l.split('\t')
        if t[0] == chrom: a[int(t[1]):int(t[2])] = True
    return a
S = {'segdup': cov(f'{sdir}/GRCh38_segdups.bed.gz', 'chr22'),
     'lowmap': cov(f'{sdir}/GRCh38_lowmappabilityall.bed.gz', 'chr22'),
     'TR/homopolymer': cov(f'{sdir}/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz', 'chr22')}
S['outside high-confidence'] = ~cov(conf, 'CM000684.2')
def load(f):
    r = {}
    for l in open(f):
        t = l.rstrip('\n').split('\t'); tg = dict((x[:2], x[5:]) for x in t[12:] if len(x) > 5 and x[2] == ':')
        if tg.get('tp') == 'P': r[t[0]] = dict(key=(t[5], t[7], t[4]), q=(t[2], t[3]), ts=int(t[7]), te=int(t[8]), AS=int(tg['AS']), cg=tg.get('cg', ''))
    return r
for tech in ['hifi', 'ont', 'clr']:
    a = load(f'{d}/{tech}_geo0.paf'); b = load(f'{d}/{tech}_geo1.paf')
    sh = [k for k in a if k in b and a[k]['key'] == b[k]['key']]
    dif = [k for k in sh if a[k]['cg'] != b[k]['cg']]
    lo = [k for k in dif if b[k]['AS'] < a[k]['AS']]
    rel = [(a[k]['AS'] - b[k]['AS']) / max(1, a[k]['AS']) for k in lo]
    span = [k for k in lo if (a[k]['ts'], a[k]['te'], a[k]['q']) != (b[k]['ts'], b[k]['te'], b[k]['q'])]
    n = len(sh); dd = sorted(b[k]['AS'] - a[k]['AS'] for k in lo)
    print(f"\n== {tech}: shared primaries {n:,}; CIGAR differs {len(dif):,} ({100*len(dif)/n:.2f}%); "
          f"lower AS {len(lo):,} ({100*len(lo)/n:.2f}%); median dAS {dd[len(dd)//2] if dd else 0}")
    for lab, (x, y) in [('<1%', (0, .01)), ('1-5%', (.01, .05)), ('>5%', (.05, 9e9))]:
        print(f"   AS loss {lab}: {100*sum(x <= r < y for r in rel)/n:.2f}% of reads")
    print(f"   truncated or split: {100*len(span)/n:.3f}% of reads")
    for name, ks in (('all shared', sh), ('lower AS', lo)):
        print(f"   {name:10s} " + "  ".join(f"{s}={100*sum(S[s][a[k]['ts']:a[k]['te']].mean() >= .5 for k in ks)/max(1,len(ks)):.1f}%" for s in S))
