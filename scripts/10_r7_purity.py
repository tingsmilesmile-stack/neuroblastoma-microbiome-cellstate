#!/usr/bin/env python3
# ---- circularity / purity analysis ----
import os, numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(BASE,"results")
gsva = pd.read_csv(os.path.join(RES,"gsva_ssgsea_scores.tsv"), sep="\t", index_col=0)
axis = pd.read_csv(os.path.join(RES,"main153_mi_axis.tsv"), sep="\t").set_index("sample")
clin = pd.read_csv(os.path.join(RES,"sample_clinical_supp.tsv"), sep="\t").set_index("sample")
INFL=["TLR_signaling","NOD_like_receptor","LPS_inflammatory"]
IMM=["CD8_T_cells","Cytotoxic_NK","Tregs","Macrophage_M1","Macrophage_M2","B_cells","DC","MHC_class_I"]
S=list(gsva.columns); X=axis.loc[S,"MI_axis_z"].values
purity=pd.to_numeric(clin.loc[S,"purity_pct"], errors="coerce")
Z=gsva.loc[IMM].apply(lambda r:(r-r.mean())/r.std(ddof=1), axis=1)
print(f"[info] 纯度可得 {purity.notna().sum()}/{len(S)}（{purity.notna().mean():.1%}）")
rows=[]
for c in IMM:
    y=Z.loc[c].values; ok=purity.notna().values
    # 基线
    b0=sm.OLS(y,sm.add_constant(X)).fit().params[1]
    # 加纯度
    dfm=pd.DataFrame({"y":y,"x":X,"p":purity.values}).dropna()
    m=sm.OLS(dfm["y"],sm.add_constant(dfm[["x","p"]])).fit()
    rows.append({"signature":c,"beta_raw":round(b0,4),
                 "beta_adj_purity":round(m.params["x"],4),"p_adj":m.pvalues["x"],
                 "beta_purity":round(m.params["p"],4),"p_purity":m.pvalues["p"],
                 "n_with_purity":len(dfm)})
r=pd.DataFrame(rows); r.to_csv(os.path.join(RES,"main153_purity_adjust.tsv"),sep="\t",index=False)
print(r.round(4).to_string(index=False))
print("\n→ 加入纯度后 MI 轴系数仍显著的签名数:",
      int(((r['p_adj']<0.05)).sum()), "/8")
