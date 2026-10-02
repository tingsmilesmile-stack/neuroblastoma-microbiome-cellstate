#!/usr/bin/env python3
# 15_pathway_cox153.py — 153 例：逐通路连续 Cox（Fig4 输入）+ 逐通路存活关联
import os, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(BASE,"results")
g=pd.read_csv(os.path.join(RES,"gsva_ssgsea_scores.tsv"),sep="\t",index_col=0)
cl=pd.read_csv(os.path.join(RES,"sample_clinical_supp.tsv"),sep="\t").set_index("sample")
MI=["TLR_signaling","NOD_like_receptor","Antimicrobial_peptides","LPS_inflammatory",
    "SCFA_butyrate_response","Bile_acid_metabolism","Tryptophan_IDO_kynurenine","Mucin_barrier"]
d=pd.DataFrame({"OS_years":pd.to_numeric(cl["OS_time"],errors="coerce")/365.25,
                "OS_event":pd.to_numeric(cl["OS_event"],errors="coerce")}).loc[g.columns]
rows=[]
for p in MI:
    z=(g.loc[p]-g.loc[p].mean())/g.loc[p].std(ddof=1)
    dd=d.join(z.rename("x")).dropna()
    c=CoxPHFitter().fit(dd,"OS_years","OS_event"); s=c.summary.loc["x"]
    rows.append({"pathway":p,"HR":np.exp(s["coef"]),"CI_lo":s["exp(coef) lower 95%"],
                 "CI_hi":s["exp(coef) upper 95%"],"p":s["p"]})
r=pd.DataFrame(rows); r["q_BH"]=multipletests(r["p"],method="fdr_bh")[1]
r=r.sort_values("p"); r.round(6).to_csv(os.path.join(RES,"main153_pathway_cox.tsv"),sep="\t",index=False)
print(r.round(4).to_string(index=False))
print(f"\n名义 p<0.05: {(r['p']<0.05).sum()}/8 ｜ FDR<0.05: {(r['q_BH']<0.05).sum()}/8")
print("[done] 15_pathway_cox153")
