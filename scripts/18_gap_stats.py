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
        E50=('Efinal','median'),E99=('Efinal',lambda s:s.quantile(.99)),Emax=('Efinal','max'))
    print(tab.round(3).to_string())
