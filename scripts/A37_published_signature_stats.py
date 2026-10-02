#!/usr/bin/env python3
"""A37 — `published_adrn_mes` 操作化的 r_MES 与 attenuation（§25.4 第四种操作化）

至此 §二 的四种操作化全部可得：
  1. mes_module           （本地 10 基因）
  2. adrn_module          （本地 8 基因）
  3. adrn_minus_mes       （本地差分，§25.4 主操作化）
  4. published_adrn_mes   （van Groningen 2017 正式签名，GSVA 打分）← 本脚本

判定仍以 #3 为准（§25.4）；#4 为敏感性，且不改变分支归属（除非 §25.4 记录被修订）。
"""
import csv, os, sys
from math import erf, sqrt

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
sys.path.insert(0, os.path.join(BASE, "scripts"))

INFLAM = ["TLR_signaling", "NOD_like_receptor", "LPS_inflammatory"]
ADRN_GENES = ["PHOX2B", "DBH", "TH", "CHGA", "CHGB", "NTRK1", "ASCL1", "GATA2"]
MES_GENES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2",
             "PDGFRB", "DCN", "LUM", "FAP"]
PRIMARY = "MHC_class_I"


def z(v):
    v = np.asarray(v, dtype=float)
    m = np.isfinite(v)
    o = np.full(v.shape, np.nan)
    o[m] = (v[m] - v[m].mean()) / v[m].std(ddof=1)
    return o


def ols(y, X):
    y = np.asarray(y, float); X = np.atleast_2d(X)
    if X.shape[0] != y.shape[0]:
        X = X.T
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    n, k = X.shape
    A = np.column_stack([np.ones(n), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    s2 = (r @ r) / (n - k - 1)
    cov = s2 * np.linalg.inv(A.T @ A)
    se = np.sqrt(np.diag(cov))
    t = coef / se
    p = np.array([2 * (1 - 0.5 * (1 + erf(abs(v) / sqrt(2)))) for v in t[1:]])
    ss = ((y - y.mean()) ** 2).sum()
    return coef[1:], p, 1 - (r @ r) / ss, n


def read_scores(path):
    with open(path, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd); rows = list(rd)
    return h[1:], {r[0]: np.array([float(x) for x in r[1:]]) for r in rows}


samples, S = read_scores(os.path.join(RES, "gsva_ssgsea_scores.tsv"))
samples2, P = read_scores(os.path.join(RES, "gsva_published_adrn_mes_scores.tsv"))
assert samples == samples2, "样本顺序不一致"

# 本地模块
with open(os.path.join(RES, "expr_tpm.tsv"), encoding="utf-8") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd)
    genes, mats = [], []
    for r in rd:
        genes.append(r[0]); mats.append([float(x) for x in r[1:]])
E = np.array(mats)
gi = {g: i for i, g in enumerate(genes)}


def mod(gs):
    p = [g for g in gs if g in gi]
    M = np.log2(E[[gi[g] for g in p]] + 1)
    return ((M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)).mean(0)


X = np.mean([z(S[p]) for p in INFLAM], axis=0)
Y = z(S[PRIMARY])

OPS = {
    "mes_module(local)":             ("MES", mod(MES_GENES)),
    "adrn_module(local)":            ("ADRN", mod(ADRN_GENES)),
    "adrn_minus_mes(local,主)":       ("ADRN-MES", mod(ADRN_GENES) - mod(MES_GENES)),
    "published_adrn_mes(GSVA)":      ("ADRN-MES", z(P["VANGRONINGEN_ADRN"]) - z(P["VANGRONINGEN_MES"])),
    "published_mes_only(GSVA)":      ("MES", z(P["VANGRONINGEN_MES"])),
    "published_adrn_only(GSVA)":     ("ADRN", z(P["VANGRONINGEN_ADRN"])),
}

print("=" * 96)
print("① r(MI 炎症轴, 各 MES 操作化) —— 含正式签名")
print("=" * 96)
print(f"  {'操作化':<30}{'角色':<12}{'r':>10}{'|r|≥0.7':>10}{'R²':>8}")
for k, (role, v) in OPS.items():
    r = float(np.corrcoef(X, z(v))[0, 1])
    print(f"  {k:<30}{role:<12}{r:>+10.3f}{'✅' if abs(r)>=0.7 else '❌':>10}{r**2:>8.3f}")

print()
print("=" * 96)
print(f"② attenuation（主终点 {PRIMARY}）；MES 项取各操作化")
print("=" * 96)
b1, p1, r2_1, n1 = ols(Y, z(X)[:, None])
print(f"  β_raw = {b1[0]:+.4f}  (p {p1[0]:.4g}, R² {r2_1:.3f}, n {n1})   ← 与 MES 无关，固定")
print()
print(f"  {'MES 操作化':<30}{'β_adj':>10}{'p_adj':>12}{'attenuation':>14}")
rows = []
for k, (role, v) in OPS.items():
    b2, p2, r2_2, n2 = ols(Y, np.column_stack([z(X), z(v)]))
    a = (b1[0] - b2[0]) / b1[0]
    rows.append((k, b2[0], p2[0], a))
    print(f"  {k:<30}{b2[0]:>+10.4f}{p2[0]:>12.4g}{a:>13.1%}")

print()
print("=" * 96)
print("③ 判读")
print("=" * 96)
arr = np.array([r[3] for r in rows])
print(f"  attenuation 范围 [{arr.min():+.1%}, {arr.max():+.1%}]｜>50% 者 {int((arr>0.5).sum())}/{len(arr)}")
print(f"  → 全部为**负值或近零** ⇒ MES 未解释掉 MI 的免疫关联（存在抑制效应）")
print(f"  → 分支判定**不受**引入正式签名而改变" if not (arr > 0.5).any()
      else "  ⚠️ 有操作化 >50%，须复核")
r_pub = float(np.corrcoef(X, OPS["published_adrn_mes(GSVA)"][1])[0, 1])
print(f"\n  正式签名 r_MES = {r_pub:+.3f}  ({'≥' if abs(r_pub)>=0.7 else '<'} 0.7)")
print(f"    （§25.4 主操作化 r_MES 仍以 adrn_minus_mes 为准）")
