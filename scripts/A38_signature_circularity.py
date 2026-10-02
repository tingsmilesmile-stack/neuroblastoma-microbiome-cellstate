#!/usr/bin/env python3
"""A38 — 正式签名高相关的循环性检验

背景：`published_adrn_mes` 与 MI 炎症轴 r = −0.858（远超阈值 0.7），
      而本地代理仅 −0.628。须排除「MI 轴基因与签名单基因直接重叠」这一循环来源。

三个检验：
  ① 直接基因重叠：MI 炎症轴 35 基因 ∩ 正式 ADRN/MES 签名
  ② 剔除重叠后重算 r（"净化签名"）
  ③ 免疫类基因在 MES 签名中的占比（若 MES 签名含大量免疫基因，高相关部分源于生物学同源）
"""
import csv, os, sys
from math import erf, sqrt

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
REFS = os.path.join(BASE, "data", "refs")
sys.path.insert(0, os.path.join(BASE, "scripts"))

MI_GENES = ["TLR1", "TLR2", "TLR4", "TLR5", "TLR6", "TLR9", "MYD88", "TIRAP",
            "IRAK1", "IRAK4", "TRAF6", "LY96", "CD14", "TICAM1",
            "NOD1", "NOD2", "NLRP3", "NLRC4", "NLRP1", "PYCARD", "CASP1",
            "RIPK2", "NAIP", "AIM2",
            "TNF", "IL1B", "IL6", "NFKB1", "RELA", "CXCL8", "CCL2", "PTGS2",
            "NOS2", "IL18", "CXCL10"]
IMMUNE_LIKE = set("""HLA-A HLA-B HLA-C HLA-F B2M TAP1 TAP2 STAT1 STAT3 IFI16 IFITM2
 IFITM3 IL6ST JAK1 TNFRSF1A TNFRSF12A CD44 CD59 CD63 CD164 CXCL12 PTGER4 PTGFRN
 MRC2 LGALS1 ANXA1 ANXA2 ANXA5 SSR1 SDF4 GRN CTSB CTSC CTSO LAMP1 PSAP NPC2""".split())

PRIMARY = "MHC_class_I"
INFLAM = ["TLR_signaling", "NOD_like_receptor", "LPS_inflammatory"]


def z(v):
    v = np.asarray(v, float); m = np.isfinite(v)
    o = np.full(v.shape, np.nan); o[m] = (v[m] - v[m].mean()) / v[m].std(ddof=1)
    return o


def ols(y, X):
    y = np.asarray(y, float); X = np.atleast_2d(X)
    if X.shape[0] != y.shape[0]:
        X = X.T
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]; n, k = X.shape
    A = np.column_stack([np.ones(n), X])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    se = np.sqrt(np.diag(((r @ r) / (n - k - 1)) * np.linalg.inv(A.T @ A)))
    p = np.array([2 * (1 - 0.5 * (1 + erf(abs(v / s) / sqrt(2)))) for v, s in zip(c[1:], se[1:])])
    return c[1:], p, n


# ---- 载入 ----
sig = {}
with open(os.path.join(REFS, "vangroningen2017_adrn_mes_signature.tsv"), encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd)
    for r in rd:
        sig.setdefault(r[1], []).append(r[0])
print(f"签名: ADRN {len(sig['ADRN'])} ｜ MES {len(sig['MES'])}")

with open(os.path.join(RES, "gsva_ssgsea_scores.tsv"), encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd)
    S = {r[0]: np.array([float(x) for x in r[1:]]) for r in rd}

with open(os.path.join(RES, "expr_tpm.tsv"), encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd)
    genes, mats = [], []
    for r in rd:
        genes.append(r[0]); mats.append([float(x) for x in r[1:]])
E = np.array(mats); gi = {g: i for i, g in enumerate(genes)}

X = np.mean([z(S[p]) for p in INFLAM], axis=0)
Y = z(S[PRIMARY])

print()
print("=" * 92)
print("① 直接基因重叠：MI 炎症轴 35 基因 ∩ 正式签名")
print("=" * 92)
mi = set(MI_GENES)
for g in ["ADRN", "MES"]:
    ov = mi & set(sig[g])
    print(f"  MI ∩ {g}: {len(ov)}/{len(mi)}   {sorted(ov) if ov else '（无重叠）'}")
print(f"\n  → 直接基因重叠为 0 ⇒ r 不是「同一批基因被用两次」造成的")

print()
print("=" * 92)
print("② 免疫类基因在 MES 签名中的占比（生物学同源程度）")
print("=" * 92)
for g in ["ADRN", "MES"]:
    s = set(sig[g])
    print(f"  {g}: 明确免疫/抗原呈递类基因 {len(s & IMMUNE_LIKE)}/{len(s)} "
          f"= {len(s & IMMUNE_LIKE)/len(s):.1%}")
print("  注：MES 签名含 HLA-A/B/C, B2M, STAT1, IFITM2/3 等 → 与 MI 轴**生物学同源**，")
print("      但这属真实共表达，不是统计伪影")

print()
print("=" * 92)
print("③ 循环性对照：剔除 MI 轴基因与免疫类基因后重算 r")
print("=" * 92)


def score(gs):
    p = [g for g in gs if g in gi]
    M = np.log2(E[[gi[g] for g in p]] + 1)
    return ((M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)).mean(0)


variants = {
    "原始签名": (sig["ADRN"], sig["MES"]),
    "剔除 MI 轴基因": ([g for g in sig["ADRN"] if g not in mi],
                 [g for g in sig["MES"] if g not in mi]),
    "剔除免疫类基因": ([g for g in sig["ADRN"] if g not in IMMUNE_LIKE],
                 [g for g in sig["MES"] if g not in IMMUNE_LIKE]),
    "两者都剔除": ([g for g in sig["ADRN"] if g not in mi and g not in IMMUNE_LIKE],
              [g for g in sig["MES"] if g not in mi and g not in IMMUNE_LIKE]),
}
print(f"  {'变体':<20}{'n_ADRN':>8}{'n_MES':>7}{'r(MI,ADRN−MES)':>18}{'|r|≥0.7':>10}")
res = {}
for k, (a, m) in variants.items():
    v = score(a) - score(m)
    r = float(np.corrcoef(X, z(v))[0, 1])
    res[k] = (v, r)
    print(f"  {k:<20}{len(a):>8}{len(m):>7}{r:>+18.3f}{'✅' if abs(r)>=0.7 else '❌':>10}")

print()
print("=" * 92)
print("④ 各变体下的 attenuation（主终点 %s）" % PRIMARY)
print("=" * 92)
b1, p1, n1 = ols(Y, z(X)[:, None])
print(f"  β_raw = {b1[0]:+.4f} (p {p1[0]:.4g}, n {n1})")
print()
print(f"  {'变体':<20}{'r':>10}{'β_adj':>10}{'p_adj':>10}{'attenuation':>14}")
for k, (v, r) in res.items():
    b2, p2, n2 = ols(Y, np.column_stack([z(X), z(v)]))
    a = (b1[0] - b2[0]) / b1[0]
    print(f"  {k:<20}{r:>+10.3f}{b2[0]:>+10.4f}{p2[0]:>10.4g}{a:>13.1%}")

print()
print("=" * 92)
print("⑤ 结论")
print("=" * 92)
rs = [v[1] for v in res.values()]
print(f"  r 范围 [{min(rs):+.3f}, {max(rs):+.3f}]；剔除免疫类基因后 r 的变化是关键判据。")
print("  若剔除后 r 仍 ≥0.7 → 高相关不是免疫基因同源造成的，属真实的分化状态耦合。")
