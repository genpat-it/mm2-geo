#!/usr/bin/env python3
"""Rounds per gap, first-pass rejection, joint (d,E) and DP-cell ratio from the MM2_GEO_GAPLOG files."""
import sys, os, numpy as np, pandas as pd
d=sys.argv[1]
for t in ['hifi','ont','clr','wg138k']:
    f1=os.path.join(d,f'{t}_geo1.tsv.gz'); f0=os.path.join(d,f'{t}_geo0.tsv.gz')
    if not os.path.exists(f1): continue
    g=pd.read_csv(f1,sep='\t'); g=g[g.path=='geo']; m=pd.read_csv(f0,sep='\t'); n=len(g)
    vc=g.rounds.value_counts().sort_index()
    print(f"\n== {t}: gaps {n:,}; rounds " + "  ".join(f"{k}:{100*v/n:.3f}%" for k,v in vc.items()))
    print(f"   mean rounds {g.rounds.mean():.4f}; reached minimap2's band {int(g.hit_ceil.sum())}; "
          f"gap-fill DP cells stock/geo {m.cells_total.sum()/g.cells_total.sum():.2f}x")
    g['dbin']=pd.cut(g.d,[-1,5,20,50,200,np.inf],labels=['0-5','5-20','20-50','50-200','>200'])
    tab=g.groupby('dbin',observed=False).agg(pct=('d',lambda s:100*len(s)/n),
        first_pass_rejected=('E1',lambda s:100*(s>=g.loc[s.index,'b0']).mean()),
        mean_rounds=('rounds','mean'),max_rounds=('rounds','max'),
        E50=('Efinal','median'),E90=('Efinal',lambda s:s.quantile(.90)),E99=('Efinal',lambda s:s.quantile(.99)),Emax=('Efinal','max'))
    print(tab.round(3).to_string())

print("\n### Certified mode: how each gap is accepted (share of gaps) and DP cells relative to stock")
print(f"{'set':8s} {'gaps':>11s} {'one-gap':>9s} {'1st band':>9s} {'at b*':>8s} {'b_max':>8s} {'cells stock/cert':>17s}")
for t in ['hifi','ont','clr','wg138k']:
    fc=os.path.join(d,f'{t}_cert.tsv.gz'); f0=os.path.join(d,f'{t}_geo0.tsv.gz')
    if not os.path.exists(fc): continue
    c=pd.read_csv(fc,sep='\t'); c=c[c.path=='geo']; m=pd.read_csv(f0,sep='\t'); n=len(c)
    sc=(c.cells_band==0)&(c.rounds==1)                     # accepted by the exact one-gap shortcut, no DP
    cap=(~sc)&(c.bfinal>=c.ceil)                           # computed at minimap2's own band
    first=(~sc)&(~cap)&(c.bfinal==c.b0)                    # certified at the initial band
    bstar=(~sc)&(~cap)&(c.bfinal>c.b0)                     # realigned once at the certified band b*
    pc=lambda x:100*x.sum()/n
    print(f"{t:8s} {n:11,d} {pc(sc):8.2f}% {pc(first):8.2f}% {pc(bstar):7.3f}% {pc(cap):7.3f}% {m.cells_total.sum()/c.cells_total.sum():16.2f}x")
