#!/usr/bin/env python3
"""For reads whose primary alignment differs between minimap2 (a) and another mode (b), count the read bases
placed at their true reference position according to the pbsim3 MAF (first 's' line = reference)."""
import re, sys
def load(f):
    d = {}
    for l in open(f):
        t = l.rstrip('\n').split('\t'); tg = dict((x[:2], x[5:]) for x in t[12:] if len(x) > 5 and x[2] == ':')
        if tg.get('tp') == 'P': d[t[0]] = dict(L=int(t[1]), qs=int(t[2]), qe=int(t[3]), st=t[4], ts=int(t[7]), AS=int(tg['AS']), cg=tg['cg'])
    return d
def amap(r):
    m = {}; q = r['qs'] if r['st'] == '+' else r['L'] - r['qe']; t = r['ts']
    for n, op in re.findall(r'(\d+)([MID])', r['cg']):
        n = int(n)
        if op == 'M':
            for i in range(n): m[q + i] = t + i
            q += n; t += n
        elif op == 'I': q += n
        else: t += n
    return m
a = load(sys.argv[1]); b = load(sys.argv[2]); maf = sys.argv[3]
dif = {k for k in a if k in b and a[k]['cg'] != b[k]['cg']}
truth = {}; block = []
def flush():
    if len(block) >= 2 and block[1][1] in dif:
        ref, rd = block[0], block[1]; m = {}; ri = int(ref[2]); qi = 0
        for x, y in zip(ref[6], rd[6]):
            if x != '-' and y != '-': m[qi] = ri
            if x != '-': ri += 1
            if y != '-': qi += 1
        truth[rd[1]] = (m, rd[4])
for l in open(maf):
    if l.startswith('a'): flush(); block = []
    elif l.startswith('s '): block.append(l.split())
flush()
tot = ca = cb = 0; win = {'minimap2': 0, 'other': 0, 'tie': 0}
for k in sorted(truth):
    T, strand = truth[k]
    if (strand == '-') != (a[k]['st'] == '-'): continue   # orientation must agree with the simulated strand
    ma, mb = amap(a[k]), amap(b[k])
    x = sum(1 for i, r in T.items() if ma.get(i) == r); y = sum(1 for i, r in T.items() if mb.get(i) == r)
    ca += x; cb += y; tot += len(T); win['minimap2' if x > y else 'other' if y > x else 'tie'] += 1
print(f"differing reads {len(dif)}; evaluated {sum(win.values())}; bases at the true position: minimap2 {ca} "
      f"({100*ca/max(1,tot):.2f}%), other {cb} ({100*cb/max(1,tot):.2f}%); closer to the truth: {win}")
