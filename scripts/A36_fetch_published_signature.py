#!/usr/bin/env python3
"""A36 — 取得并验证 van Groningen 2017 ADRN/MES 正式签名（`published_adrn_mes`）

来源（已核）：
  van Groningen T, et al. "Neuroblastoma is composed of two super-enhancer-associated
  differentiation states." Nat Genet 2017;49(8):1261-1266.
  PMID 28650485 ｜ DOI 10.1038/ng.3899
  补充材料 Table S2 = "List of ADRN and MES signature genes"
  URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3899/
       MediaObjects/41588_2017_BFng3899_MOESM3_ESM.xlsx

产出：
  data/refs/vangroningen2017_adrn_mes_signature.tsv   （gene, group）
  results/published_adrn_mes_coverage.tsv             （覆盖率）
  results/published_adrn_mes_scores.tsv               （GSVA 打分，若 R 可用）
"""
import csv, os, subprocess, sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
RES = os.path.join(BASE, "results")
REFS = os.path.join(DATA, "refs")
XLSX = os.path.join(REFS, "ng3899_MOESM3.xlsx")
OUT_SIG = os.path.join(REFS, "vangroningen2017_adrn_mes_signature.tsv")

import pandas as pd

print("=" * 90)
print("① 从 Table S2 提取签名")
print("=" * 90)
raw = pd.read_excel(XLSX, sheet_name=0, header=None)
print(f"  原始表: {raw.shape}  首行: {str(raw.iloc[0,0])[:80]}")
df = raw.iloc[2:, :2].copy()
df.columns = ["gene", "group"]
df = df.dropna(subset=["gene"])
df["gene"] = df["gene"].astype(str).str.strip()
df["group"] = df["group"].astype(str).str.strip()
df = df[df["group"].isin(["ADRN", "MES"])]
df = df.drop_duplicates(subset=["gene"])
n_adrn = int((df.group == "ADRN").sum())
n_mes = int((df.group == "MES").sum())
print(f"  ADRN {n_adrn} ｜ MES {n_mes} ｜ 合计 {len(df)}")
assert n_adrn > 300 and n_mes > 400, "基因数异常，须复核"
df.to_csv(OUT_SIG, sep="\t", index=False)
print(f"  [out] {os.path.relpath(OUT_SIG, BASE)}")

print()
print("=" * 90)
print("② 与本地矩阵的覆盖核对（20,584 基因）")
print("=" * 90)
with open(os.path.join(RES, "expr_tpm.tsv"), encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t")
    hdr = next(rd)
    n_samples = len(hdr) - 1
    matrix_genes = set()
    for r in rd:
        if r and r[0]:
            matrix_genes.add(r[0].strip())
print(f"  矩阵: {len(matrix_genes)} 基因 × {n_samples} 样本")

rows = []
for g in ["ADRN", "MES"]:
    genes = df[df.group == g]["gene"].tolist()
    present = [x for x in genes if x in matrix_genes]
    missing = [x for x in genes if x not in matrix_genes]
    rows.append(dict(group=g, n_signature=len(genes), n_present=len(present),
                     coverage=len(present) / len(genes), n_missing=len(missing)))
    print(f"  {g}: {len(present)}/{len(genes)} = {len(present)/len(genes):.1%}  "
          f"缺失 {len(missing)}")
    if missing:
        print(f"     缺失样例: {missing[:14]}")

cov = pd.DataFrame(rows)
cov.to_csv(os.path.join(RES, "published_adrn_mes_coverage.tsv"), sep="\t", index=False)
print(f"\n  [out] results/published_adrn_mes_coverage.tsv")

# 阈值判定：GSVA 对低覆盖集不稳定，设 60% 警戒线
for _, r in cov.iterrows():
    if r["coverage"] < 0.60:
        print(f"  ⚠️ {r['group']} 覆盖率 {r['coverage']:.1%} < 60% —— "
              f"打分可能不稳定，须在 Methods 报告")

print()
print("=" * 90)
print("③ 生成 GSVA 打分脚本输入（16 集 + ADRN/MES 两签名）")
print("=" * 90)
sets_dir = os.path.join(DATA, "gene_sets")
os.makedirs(sets_dir, exist_ok=True)
gmt = os.path.join(sets_dir, "vangroningen_adrn_mes.gmt")
with open(gmt, "w", encoding="utf-8") as fh:
    for g in ["ADRN", "MES"]:
        genes = [x for x in df[df.group == g]["gene"].tolist() if x in matrix_genes]
        fh.write("VANGRONINGEN_" + g + "\tna\t" + "\t".join(genes) + "\n")
print(f"  [out] {os.path.relpath(gmt, BASE)}")
print("  下一步：用 03_gsva_ssgsea.R 的同一参数（ssgseaParam, alpha=0.25, normalize=TRUE）对这两个集打分")

print()
print("=" * 90)
print("④ 交叉核对：本签名 vs 本地代理基因集")
print("=" * 90)
LOCAL_ADRN = ["PHOX2B", "DBH", "TH", "CHGA", "CHGB", "NTRK1", "ASCL1", "GATA2"]
LOCAL_MES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2",
             "PDGFRB", "DCN", "LUM", "FAP"]
adrn_sig = set(df[df.group == "ADRN"]["gene"])
mes_sig = set(df[df.group == "MES"]["gene"])
print(f"  本地 ADRN 8 基因在正式 ADRN 签名内: {sum(g in adrn_sig for g in LOCAL_ADRN)}/8 "
      f"→ {[g for g in LOCAL_ADRN if g in adrn_sig]}")
print(f"  本地 MES 10 基因在正式 MES 签名内: {sum(g in mes_sig for g in LOCAL_MES)}/10 "
      f"→ {[g for g in LOCAL_MES if g in mes_sig]}")
print(f"  两签名重叠基因: {len(adrn_sig & mes_sig)} "
      f"→ {sorted(adrn_sig & mes_sig)[:10]}")
print("\n  说明：本地代理集是**手工精选标记基因**，正式签名是**超级增强子关联基因集**")
print("     两者维度差一个量级（8/10 vs 369/485），故 r 与 attenuation 须全报。")
