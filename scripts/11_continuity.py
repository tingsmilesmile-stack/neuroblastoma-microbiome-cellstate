#!/usr/bin/env python3
# 11_continuity.py — R8/R9 数字：PC1 方差占比、标签↔PC1 一致性、通路相关结构、SCFA 轴
import os, numpy as np, pandas as pd
from scipy import stats
BASE=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=os.path.join(BASE,"results")
g=pd.read_csv(os.path.join(RES,"gsva_ssgsea_scores.tsv"),sep="\t",index_col=0)
lab=pd.read_csv(os.path.join(RES,"main153_cluster_labels.tsv"),sep="\t").set_index("sample")
MI=["TLR_signaling","NOD_like_receptor","Antimicrobial_peptides","LPS_inflammatory",
    "SCFA_butyrate_response","Bile_acid_metabolism","Tryptophan_IDO_kynurenine","Mucin_barrier"]
Z=g.loc[MI].apply(lambda r:(r-r.mean())/r.std(ddof=1),axis=1)   # 通路 × 样本
X=Z.T.values; lab2=lab.loc[Z.columns,"k2"]

# PC1
Xc=X-X.mean(0); U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
var=(S**2)/np.sum(S**2)*100
pc1=Xc@Vt[0]; 
print(f"① PC1 解释方差 = {var[0]:.1f}% ｜ PC2 = {var[1]:.1f}%")
r_pb=stats.pointbiserialr((lab2=="MI_high").astype(int), pc1)
print(f"   标签(MI_high) ↔ PC1: point-biserial r = {r_pb[0]:.3f}, p = {r_pb[1]:.2e}")

# 通路相关结构
C=Z.T.corr()
infl=["TLR_signaling","NOD_like_receptor","LPS_inflammatory"]
ii=[C.loc[a,b] for i,a in enumerate(infl) for b in infl[i+1:]]
print(f"\n② 炎症组内 pairwise r: {min(ii):.3f} – {max(ii):.3f}")
infl_axis=Z.loc[infl].mean(0)
for p in MI:
    r,pv=stats.pearsonr(infl_axis, Z.loc[p])
    print(f"   {p:28s} vs 炎症轴 r = {r:+.3f} (p={pv:.2e})")

# SCFA 轴
scfa=Z.loc["SCFA_butyrate_response"]; bile=Z.loc["Bile_acid_metabolism"]
print(f"\n③ r(SCFA, 炎症轴) = {stats.pearsonr(scfa,infl_axis)[0]:+.3f}")
print(f"   r(Bile, 炎症轴) = {stats.pearsonr(bile,infl_axis)[0]:+.3f}")
print(f"   r(SCFA, Bile)   = {stats.pearsonr(scfa,bile)[0]:+.3f}")

# SCFA 轴 vs 正式 ADRN/MES 签名 + 纯度
pub=pd.read_csv(os.path.join(RES,"gsva_published_adrn_mes_scores.tsv"),sep="\t",index_col=0)
pubz=pub.apply(lambda r:(r-r.mean())/r.std(ddof=1),axis=1)
if {"VANGRONINGEN_ADRN","VANGRONINGEN_MES"} <= set(pubz.index):
    mes_axis=pubz.loc["VANGRONINGEN_ADRN"]-pubz.loc["VANGRONINGEN_MES"]
    print(f"\n④ r(炎症轴, 正式签名 ADRN−MES) = {stats.pearsonr(infl_axis, mes_axis)[0]:+.3f}")
    print(f"   r(SCFA,  正式签名 ADRN−MES) = {stats.pearsonr(scfa, mes_axis)[0]:+.3f}")
clin=pd.read_csv(os.path.join(RES,"sample_clinical_supp.tsv"),sep="\t").set_index("sample")
pur=pd.to_numeric(clin.loc[Z.columns,"purity_pct"],errors="coerce")
ok=pur.notna()
print(f"   r(SCFA, 纯度) = {stats.pearsonr(scfa[ok], pur[ok])[0]:+.3f} (n={ok.sum()})")
print(f"   r(炎症轴, 纯度) = {stats.pearsonr(infl_axis[ok], pur[ok])[0]:+.3f}")
# 输出
pd.DataFrame({"pathway":MI,
  "r_with_inflammatory_axis":[round(stats.pearsonr(infl_axis,Z.loc[p])[0],3) for p in MI],
  "mean_MIhigh":[round(Z.loc[p][lab2=="MI_high"].mean(),3) for p in MI],
  "mean_MIlow":[round(Z.loc[p][lab2=="MI_low"].mean(),3) for p in MI]}
).to_csv(os.path.join(RES,"main153_pathway_corr.tsv"),sep="\t",index=False)
print("\n[done] 11_continuity")
