#!/usr/bin/env python3
# 07_main153.py — 153 例主分析下游：样本流、免疫关联、生存（含 MYCN/COG/stage）、
# 连续 MI 炎症轴、临床关联。回应 R4/R7/R9/R10/R12。
import os, json, numpy as np, pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(BASE, "results"); DATA = os.path.join(BASE, "data")
FDR = lambda p: multipletests(np.nan_to_num(p, nan=1.0), method="fdr_bh")[1]

gsva = pd.read_csv(os.path.join(RES,"gsva_ssgsea_scores.tsv"), sep="\t", index_col=0)
lab  = pd.read_csv(os.path.join(RES,"main153_cluster_labels.tsv"), sep="\t").set_index("sample")
clin = pd.read_csv(os.path.join(RES,"sample_clinical_supp.tsv"), sep="\t").set_index("sample")
S = list(gsva.columns); assert list(lab.index)==S, "标签顺序不一致"
clin = clin.loc[S]

INFL = ["TLR_signaling","NOD_like_receptor","LPS_inflammatory"]
IMM  = ["CD8_T_cells","Cytotoxic_NK","Tregs","Macrophage_M1","Macrophage_M2",
        "B_cells","DC","MHC_class_I"]
z = lambda v: (v - v.mean())/v.std(ddof=1)
Z = gsva.loc[INFL+IMM].apply(lambda r: z(r), axis=1)          # 每行跨样本 z
MI_axis = Z.loc[INFL].mean(axis=0)                            # 连续 MI 炎症轴
grp = lab["k2"]                                               # MI_high / MI_low

# ---------- ① 样本流 ----------
_pe = json.load(open(os.path.join(DATA,"primary_expr_ids.json")))
uniq_cases = sorted(set(r[2] for r in _pe))
files = sorted(os.listdir(os.path.join(DATA,"expr_all")))
print(f"[info] expr_all 条目 {len(files)}；primary files {len(_pe)}；唯一 case {len(uniq_cases)}")
flow = [
 ("GDC TARGET-NBL 全部文件", 13525, "GDC API 检索（2026-09-16）"),
 ("RNA-seq 基因表达文件（gene expression）", 162, "含 Primary/Recurrent"),
 ("— Primary Tumor", 153, "本分析对象"),
 ("— Recurrent Tumor", 8, "未纳入（非初诊）"),
 ("— Recurrent Blood Marrow", 1, "未纳入（非实体瘤）"),
 ("Primary Tumor 唯一患者数", len(uniq_cases), "无重复 aliquot"),
 ("成功下载并构建矩阵", len(S), "全部示例通过质控"),
 ("临床补充表成功链接", len(clin), "Discovery 343 / Validation 499 记录池"),
 ("最终生存分析样本", int(clin["OS_time"].notna().sum()), "OS 时间可得"),
]
pd.DataFrame(flow, columns=["step","n","note"]).to_csv(
    os.path.join(RES,"main153_sample_flow.tsv"), sep="\t", index=False)
print("① 样本流"); print(pd.DataFrame(flow, columns=["step","n","note"]).to_string(index=False))

# ---------- ② 免疫签名 vs MI 轴 / 分组 ----------
rows=[]
for c in IMM:
    r,p = stats.pearsonr(MI_axis, Z.loc[c]); rs,ps = stats.spearmanr(MI_axis, Z.loc[c])
    g1 = Z.loc[c][grp=="MI_high"]; g0 = Z.loc[c][grp=="MI_low"]
    H,pk = stats.kruskal(g1,g0)
    rows.append({"signature":c,"r_axis":r,"p_axis":p,"rho_axis":rs,
                 "MI_high_median":g1.median(),"MI_low_median":g0.median(),
                 "kruskal_p":pk,"cohens_d":(g1.mean()-g0.mean())/np.sqrt((g1.var()+g0.var())/2)})
imm = pd.DataFrame(rows); imm["FDR"] = FDR(imm["kruskal_p"])
imm.round(6).to_csv(os.path.join(RES,"main153_immune_assoc.tsv"), sep="\t", index=False)
print("\n② 免疫签名 vs MI 炎症轴/分组"); print(imm.round(4).to_string(index=False))
print(f"   FDR<0.05: {(imm['FDR']<0.05).sum()}/8 ｜FDR<0.10: {(imm['FDR']<0.10).sum()}/8")

# ---------- ③ 生存 ----------
d = clin.copy()
d["OS_years"] = pd.to_numeric(d["OS_time"], errors="coerce")/365.25
d["male"] = (d["gender"]=="male").astype(int)
d["MI_high"] = (grp=="MI_high").astype(int)
d["age_years"] = pd.to_numeric(d["age_years"], errors="coerce")
d["mycn_amp"] = pd.to_numeric(d["mycn_amp"], errors="coerce").fillna(0).astype(int)
d["high_risk"] = (d["cog_risk"].astype(str).str.contains("High", case=False, na=False)).astype(int)
d["stage4"] = d["inss_stage"].astype(str).str.contains("4", na=False).astype(int)
dv = d.dropna(subset=["OS_years","OS_event"]).copy()
dv = dv[dv["OS_years"]>=0]
print(f"\n③ 生存：n={len(dv)}｜死亡={int(dv['OS_event'].sum())}｜删失={int((1-dv['OS_event']).sum())}"
      f"｜删失者中位随访={dv.loc[dv['OS_event']==0,'OS_years'].median():.2f} 年")
print(f"   年龄单位=年（连续）；年龄范围 {dv['age_years'].min():.2f}–{dv['age_years'].max():.2f}")

lr = multivariate_logrank_test(dv["OS_years"], grp.loc[dv.index], dv["OS_event"])
km=[]
kmf=KaplanMeierFitter()
for g,sub in dv.groupby(grp.loc[dv.index]):
    kmf.fit(sub["OS_years"], sub["OS_event"])
    ci = kmf.confidence_interval_
    km.append({"subtype":g,"n":len(sub),"events":int(sub["OS_event"].sum()),
               "median_OS_yr":kmf.median_survival_time_,
               "median_CI_lo":ci.iloc[:,0].max(),"median_CI_hi":ci.iloc[:,1].max()})
kmdf=pd.DataFrame(km); kmdf.to_csv(os.path.join(RES,"main153_km_summary.tsv"),sep="\t",index=False)
print(f"   log-rank p = {lr.p_value:.4g}"); print(kmdf.round(3).to_string(index=False))

models = {"M1 subtype+age+sex":["MI_high","age_years","male"],
          "M2 +MYCN":["MI_high","age_years","male","mycn_amp"],
          "M3 +COG risk":["MI_high","age_years","male","high_risk"],
          "M4 +INSS stage4":["MI_high","age_years","male","stage4"],
          "M5 full":["MI_high","age_years","male","mycn_amp","high_risk","stage4"]}
srows=[]
for name,cols in models.items():
    dd = dv[["OS_years","OS_event"]+cols].dropna()
    try:
        cph = CoxPHFitter().fit(dd, "OS_years","OS_event")
    except Exception as e:
        print(f"   [{name}] 不可估计：{type(e).__name__}（协变量共线）")
        srows.append({"model":name,"covariate":"(model)","HR":np.nan,"CI_lo":np.nan,
                      "CI_hi":np.nan,"p":np.nan,"coef":np.nan})
        continue
    s = cph.summary
    for v in cols:
        srows.append({"model":name,"covariate":v,"HR":s.loc[v,"exp(coef)"],
                      "CI_lo":s.loc[v,"exp(coef) lower 95%"],"CI_hi":s.loc[v,"exp(coef) upper 95%"],
                      "p":s.loc[v,"p"],"coef":s.loc[v,"coef"]})
sd=pd.DataFrame(srows); sd.round(6).to_csv(os.path.join(RES,"main153_cox_models.tsv"),sep="\t",index=False)
print("\n   Cox（HR, 95%CI, p）"); print(sd.round(4).to_string(index=False))

# PH 检验（年龄 vs subtype）
from lifelines.statistics import proportional_hazard_test
ph_rows=[]
for v in ["MI_high","age_years","male"]:
    dd = dv[["OS_years","OS_event"]+["MI_high","age_years","male"]].dropna()
    cph = CoxPHFitter().fit(dd,"OS_years","OS_event")
    t = proportional_hazard_test(cph, dd, time_transform="rank")
    ph_rows.append({"covariate":v,"PH_stat":t.summary.loc[v,"test_statistic"],
                    "PH_p":t.summary.loc[v,"p"]})
phdf=pd.DataFrame(ph_rows); phdf.to_csv(os.path.join(RES,"main153_ph_test.tsv"),sep="\t",index=False)
print("\n   PH 假设检验"); print(phdf.round(4).to_string(index=False))

# ---------- ④ 连续 MI 轴 Cox ----------
dd = dv[["OS_years","OS_event"]].join(MI_axis.rename("MI_axis")).dropna()
c = CoxPHFitter().fit(dd,"OS_years","OS_event"); r=c.summary.loc["MI_axis"]
print(f"\n④ 连续 MI 炎症轴 Cox：HR={np.exp(r['coef']):.3f} "
      f"(95%CI {r['exp(coef) lower 95%']:.3f}–{r['exp(coef) upper 95%']:.3f}), p={r['p']:.4g}")

# ---------- ⑤ 临床关联 ----------
crows=[]
for col,test in [("age_years","kruskal"),("mycn_amp","chi2"),("high_risk","chi2"),
                 ("stage4","chi2"),("gender","chi2")]:
    sub = d[[col]].join(grp.rename("g")).dropna()
    if test=="kruskal":
        gs=[v[col].values for _,v in sub.groupby("g")]
        H,p = stats.kruskal(*gs); crows.append({"variable":col,"test":"Kruskal-Wallis","p":p})
    else:
        ct=pd.crosstab(sub["g"],sub[col])
        chi2,p,_,_=stats.chi2_contingency(ct); crows.append({"variable":col,"test":"chi-square","p":p})
cd=pd.DataFrame(crows); cd.round(6).to_csv(os.path.join(RES,"main153_clinical_assoc.tsv"),sep="\t",index=False)
print("\n⑤ 临床关联"); print(cd.round(4).to_string(index=False))
pd.DataFrame({"sample":S,"k2":grp.values,"MI_axis_z":MI_axis.round(4).values}).to_csv(
    os.path.join(RES,"main153_mi_axis.tsv"),sep="\t",index=False)
print("\n[done] 07_main153")
