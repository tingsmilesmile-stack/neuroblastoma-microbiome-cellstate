#!/usr/bin/env python3
"""D4 — 分支判定（预注册 §25 的唯一执行入口）

⚠️ 本脚本实现 **2026-09-29 冻结** 的判定规范：
   - Analysis protocol (pre-specified decision rules)
   - Summary of the analysis conventions
判定链（**不得更改**）：
   Y（主终点） = MHC_class_I（GSVA 标准实现，z）
   X（自变量） = MI 炎症轴 = TLR+NOD+LPS 三集 z 均值
   MES         = adrn_minus_mes（§二 主操作化）
   β_raw       = OLS(Y ~ X)
   β_adj       = OLS(Y ~ X + MES)
   attenuation = (β_raw − β_adj) / β_raw          （L0 层为判定层）

不可估计规则（§25.5）：仅当 p_raw < 0.05 **且** |β_raw| ≥ 0.10 时估计 attenuation。

输出：results/d4_branch_adjudication.tsv
"""
import csv, os, sys
from math import erf, sqrt

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
sys.path.insert(0, os.path.join(BASE, "scripts"))

# ---------------- 冻结常量（不得更改） ----------------
PRIMARY_ENDPOINT = "MHC_class_I"          # §25.2 预指定
INFLAM = ["TLR_signaling", "NOD_like_receptor", "LPS_inflammatory"]
ADRN_GENES = ["PHOX2B", "DBH", "TH", "CHGA", "CHGB", "NTRK1", "ASCL1", "GATA2"]
MES_GENES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "VIM", "ACTA2",
             "PDGFRB", "DCN", "LUM", "FAP"]
R_THRESHOLD = 0.70                        # §25.4 阈值不动
ATTEN_CUT = 0.50                          # 50%
BETA_MIN = 0.10                           # §25.5 |β_raw| 下限
P_MAX = 0.05                              # §25.5 p_raw 上限
SEED = 20260930                           # §25.7 bootstrap seed
NBOOT = 5000
IMMUNE_SETS = ["CD8_T_cells", "Cytotoxic_NK", "Tregs", "Macrophage_M1",
               "Macrophage_M2", "B_cells", "DC", "MHC_class_I"]


# ---------------- 基础工具 ----------------
def z(v):
    """按有限值标准化（nan 保持 nan）"""
    v = np.asarray(v, dtype=float)
    m = np.isfinite(v)
    out = np.full(v.shape, np.nan)
    if m.sum() > 1:
        out[m] = (v[m] - v[m].mean()) / v[m].std(ddof=1)
    return out


def pval(t):
    return 2 * (1 - 0.5 * (1 + erf(abs(t) / sqrt(2))))


def ols(y, X, names):
    """OLS，返回 β（不含截距）、SE、p、R²。X 为二维数组（不含截距）。"""
    y = np.asarray(y, dtype=float)
    X = np.atleast_2d(X)
    if X.shape[0] != y.shape[0]:
        X = X.T
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[m], X[m]
    n, k = X.shape
    if n <= k + 1:
        return None
    A = np.column_stack([np.ones(n), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    dof = n - k - 1
    s2 = (resid @ resid) / dof
    try:
        cov = s2 * np.linalg.inv(A.T @ A)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.diag(cov))
    tvals = coef / se
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - (resid @ resid) / ss_tot if ss_tot > 0 else np.nan
    return dict(beta=coef[1:], se=se[1:],
                p=np.array([pval(v) for v in tvals[1:]]),
                r2=r2, n=n, names=names)


def partial_r2(y, X_full, X_reduced):
    """偏 R²：由 X_reduced 到 X_full 的 R² 增量"""
    y = np.asarray(y, dtype=float)
    m = (np.isfinite(y) & np.all(np.isfinite(X_full), axis=1)
         & np.all(np.isfinite(X_reduced), axis=1))
    yy, Xf, Xr = y[m], np.atleast_2d(X_full)[m], np.atleast_2d(X_reduced)[m]
    if Xf.shape[0] != len(yy):
        Xf = Xf.T[m]
    if Xr.shape[0] != len(yy):
        Xr = Xr.T[m]

    def r2(X):
        A = np.column_stack([np.ones(len(yy)), X])
        c, *_ = np.linalg.lstsq(A, yy, rcond=None)
        r = yy - A @ c
        ss_tot = ((yy - yy.mean()) ** 2).sum()
        return 1 - (r @ r) / ss_tot if ss_tot > 0 else np.nan

    return r2(Xf) - r2(Xr)


def boot_atten(y, X_raw, X_adj, nboot=NBOOT, seed=SEED):
    """attenuation 的 percentile bootstrap CI"""
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    X_raw = np.atleast_2d(X_raw); X_adj = np.atleast_2d(X_adj)
    if X_raw.shape[0] != len(y):
        X_raw = X_raw.T
    if X_adj.shape[0] != len(y):
        X_adj = X_adj.T
    m = (np.isfinite(y) & np.all(np.isfinite(X_raw), axis=1)
         & np.all(np.isfinite(X_adj), axis=1))
    y, X_raw, X_adj = y[m], X_raw[m], X_adj[m]
    n = len(y)
    vals = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        a = ols(z(y[idx]), z(X_raw[idx]), ["X"])
        b = ols(z(y[idx]), np.column_stack([z(X_raw[idx]), z(X_adj[idx])]), ["X", "MES"])
        if not a or not b or a["beta"][0] == 0:
            continue
        vals.append((a["beta"][0] - b["beta"][0]) / a["beta"][0])
    if len(vals) < 100:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


# ---------------- 数据读取 ----------------
def load_scores():
    p = os.path.join(RES, "gsva_ssgsea_scores.tsv")
    if not os.path.exists(p):
        raise SystemExit("[FATAL] 缺 results/gsva_ssgsea_scores.tsv —— "
                         "判定必须用 GSVA 标准实现（§25.2），不得回退手写实现")
    with open(p, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd); rows = list(rd)
    samples = h[1:]
    S = {r[0]: np.array([float(x) for x in r[1:]]) for r in rows}
    return samples, S


def load_mes(samples):
    """adrn_module − mes_module（§二 主操作化）"""
    with open(os.path.join(RES, "expr_tpm.tsv"), encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        cols = hdr[1:]
        genes, mats = [], []
        for r in rd:
            genes.append(r[0]); mats.append([float(x) for x in r[1:]])
    E = np.array(mats)
    gi = {g: i for i, g in enumerate(genes)}
    assert cols == samples, "矩阵列与评分样本顺序不一致"

    def mod(gs, lab):
        present = [g for g in gs if g in gi]
        miss = [g for g in gs if g not in gi]
        if miss:
            print(f"  [!] {lab} 缺 {len(miss)} 基因: {miss}")
        M = np.log2(E[[gi[g] for g in present]] + 1)
        Z = (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)
        print(f"  {lab} 模块 {len(present)}/{len(gs)} 基因可用")
        return Z.mean(0)

    return mod(ADRN_GENES, "ADRN") - mod(MES_GENES, "MES")


def load_published_signature(samples):
    """van Groningen 2017 正式签名（GSVA 打分，§25.4 第四种操作化）。
    缺失时返回 None（记「未做」，不得用其它集冒充）。"""
    p = os.path.join(RES, "gsva_published_adrn_mes_scores.tsv")
    if not os.path.exists(p):
        print("  [!] 无 published_adrn_mes 打分 —— 该操作化记「未做」")
        return None
    with open(p, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd)
        if h[1:] != samples:
            raise SystemExit("[FATAL] published 签名样本顺序与主评分不一致")
        d = {r[0]: np.array([float(x) for x in r[1:]]) for r in rd}
    if "VANGRONINGEN_ADRN" not in d or "VANGRONINGEN_MES" not in d:
        return None
    return z(d["VANGRONINGEN_ADRN"]) - z(d["VANGRONINGEN_MES"])


def load_clinical(samples):
    import pandas as pd
    from lib_clinical import load_clinical as _LC, get_field, COG_MAP
    CL = _LC()
    if len(CL) == 0:
        raise SystemExit("[FATAL] 临床表 0 条记录 -- see the project troubleshooting notes")
    age_d, ok_age = get_field(CL, samples, "Age at Diagnosis in Days")
    pur, ok_pur = get_field(CL, samples, "Percent Tumor", numeric_loose=True)
    cog, ok_cog = get_field(CL, samples, "COG Risk Group", COG_MAP)
    sc = pd.read_csv(os.path.join(RES, "sample_clinical_supp.tsv"), sep="\t")
    sex_raw = sc.set_index("submitter_id").reindex(samples)["gender"].to_numpy()
    sex = np.where(pd.isna(sex_raw), np.nan,
                   np.array([1.0 if str(v).lower().startswith("m") else 0.0
                             for v in sex_raw]))
    return {
        "age": np.where(ok_age, age_d / 365.25, np.nan),
        "purity": np.where(ok_pur, pur, np.nan),
        "cog": np.where(ok_cog, cog, np.nan),
        "sex": sex,
    }


# ---------------- 主流程 ----------------
def main():
    samples, S = load_scores()
    print(f"[info] 样本 {len(samples)}｜实现 = GSVA 标准")

    Xaxis = np.mean([z(S[p]) for p in INFLAM], axis=0)
    if PRIMARY_ENDPOINT not in S:
        raise SystemExit(f"[FATAL] 预指定主终点 {PRIMARY_ENDPOINT} 不在评分中")
    Y = z(S[PRIMARY_ENDPOINT])
    MES = z(load_mes(samples))
    C = load_clinical(samples)

    print("\n" + "=" * 88)
    print("① r_MES（全操作化阶梯；判定规则不变）")
    print("=" * 88)
    r_main = float(np.corrcoef(Xaxis, MES)[0, 1])
    print(f"  adrn_minus_mes（本地差分）    r = {r_main:+.3f}   "
          f"{'≥0.7' if abs(r_main) >= R_THRESHOLD else '<0.7'}")
    MES_PUB = load_published_signature(samples)
    if MES_PUB is not None:
        r_pub = float(np.corrcoef(Xaxis, MES_PUB)[0, 1])
        print(f"  published_adrn_mes（正式签名）r = {r_pub:+.3f}   "
              f"{'≥0.7' if abs(r_pub) >= R_THRESHOLD else '<0.7'}   "
              f"← 文献标准，作为 headline 上报")
    else:
        r_pub = float("nan")
    print(f"  adrn_module / mes_module 见 Fig 3B 全阶梯")
    print(f"  ⚠️ 阈值 {R_THRESHOLD} 不动；不得挑最大者", end="")
    print("；两操作化须**同框并报**")

    print("\n" + "=" * 88)
    print("② attenuation（L0 = 判定层）")
    print("=" * 88)
    XR = z(Xaxis)[:, None]
    XA = np.column_stack([z(Xaxis), MES])
    raw = ols(Y, XR, ["X"])
    adj = ols(Y, XA, ["X", "MES"])
    if not raw or not adj:
        raise SystemExit("[FATAL] OLS 失败")

    br, pr = float(raw["beta"][0]), float(raw["p"][0])
    ba, pa = float(adj["beta"][0]), float(adj["p"][0])
    pr2 = partial_r2(Y, XA, XR)
    print(f"  β_raw = {br:+.4f}  (SE {raw['se'][0]:.4f}, p {pr:.4g}, "
          f"n={raw['n']}, R²={raw['r2']:.3f})")
    print(f"  β_adj = {ba:+.4f}  (SE {adj['se'][0]:.4f}, p {pa:.4g}, "
          f"n={adj['n']}, R²={adj['r2']:.3f})")
    print(f"  MES 项      β={adj['beta'][1]:+.4f}  p={adj['p'][1]:.4g}")
    print(f"  偏 R²（加入 MES 的增量）= {pr2:.4f}")

    estimable = (pr < P_MAX) and (abs(br) >= BETA_MIN)
    if estimable:
        att = (br - ba) / br
        # FIX: boot_atten() internally builds column_stack([z(X_raw), z(X_adj)]).
        # Passing the full adjusted design matrix XA (which already contains the
        # X column) duplicated that column -> rank-deficient design [X, X, MES],
        # so the bootstrap split the coefficient arbitrarily and the CI was invalid.
        # Correct usage: pass only the single covariate column (MES) here.
        lo, hi = boot_atten(Y, XR, MES[:, None])
        print(f"\n  attenuation = {att:+.1%}   95%CI [{lo:+.1%}, {hi:+.1%}]"
              f"   (bootstrap {NBOOT}, seed {SEED})")
    else:
        att = lo = hi = float("nan")
        why = []
        if pr >= P_MAX:
            why.append(f"p_raw={pr:.4g} ≥ {P_MAX}")
        if abs(br) < BETA_MIN:
            why.append(f"|β_raw|={abs(br):.4f} < {BETA_MIN}")
        print(f"\n  attenuation not estimable —— {'; '.join(why)}")
        print("    §25.5：" + ("β_raw 不显著 → 判定 A" if pr >= P_MAX
                              else "β_raw 过小 → 不报比值，改用 β_adj 显著性"))

    print("\n" + "=" * 88)
    print("③ 最终判定（§25.6）")
    print("=" * 88)
    big = abs(r_main) >= R_THRESHOLD or (np.isfinite(r_pub) and abs(r_pub) >= R_THRESHOLD)
    if not estimable:
        if pr >= P_MAX:
            branch, reason = "A", "β_raw 不显著 → 连未校正关联都不成立"
        else:
            branch = "B" if pa < P_MAX else "A"
            reason = "β_raw 过小 → 改按 β_adj 显著性判定"
    else:
        if att > ATTEN_CUT or pa >= P_MAX:
            branch = "A+" if big else "A"
            reason = f"衰减 {att:.1%} >50% 或 β_adj 不显著"
        else:
            branch = "B"
            reason = f"衰减 {att:.1%} ≤50% 且 β_adj 仍显著"
            if big:
                reason += "（冲突：|r|≥0.7 → 衰减优先，正文用 'partially independent'）"

    hi_mask = Xaxis >= np.nanmedian(Xaxis)
    d = float(np.nanmean(MES[hi_mask]) - np.nanmean(MES[~hi_mask]))
    if abs(r_main) < 0.3 and abs(d) < 0.3:
        branch, reason = "C", f"|r|={abs(r_main):.3f}<0.3 且 |d|={abs(d):.3f}<0.3"

    print(f"  |r_MES| = {abs(r_main):.3f}（本地）｜"
          f"{abs(r_pub):.3f}（正式签名）｜ attenuation = "
          f"{'n/a' if not estimable else f'{att:.1%}'} ｜ 偏 R² = {pr2:.4f}")
    print(f"  MI-high vs MI-low 的 MES 效应量  Cohen's d = {d:+.3f}")
    print(f"\n  判定：分支 {branch}")
    print(f"      依据：{reason}")

    print("\n" + "=" * 88)
    print("④ 对称协变量阶梯（L1–L4 敏感性；判定仍取 L0）")
    print("=" * 88)
    layers = [("L0 主判定", []), ("L1 +purity", ["purity"]), ("L2 +age", ["age"]),
              ("L3 +COG", ["cog"]),
              ("L4 +purity+age+COG+sex", ["purity", "age", "cog", "sex"])]
    rows_out = []
    print(f"  {'层':<26}{'β_raw':>10}{'β_adj':>10}{'衰减':>10}{'p_adj':>10}")
    for lab, covs in layers:
        cr = [z(Xaxis)] + [z(C[c]) for c in covs]
        ca = [z(Xaxis), MES] + [z(C[c]) for c in covs]
        R = ols(Y, np.column_stack(cr), ["X"] + covs)
        A = ols(Y, np.column_stack(ca), ["X", "MES"] + covs)
        if not R or not A:
            print(f"  {lab:<26}{'—':>10}{'—':>10}{'—':>10}{'—':>10}")
            continue
        b1, b2 = float(R["beta"][0]), float(A["beta"][0])
        a = (b1 - b2) / b1 if b1 != 0 else float("nan")
        print(f"  {lab:<26}{b1:>+10.4f}{b2:>+10.4f}{a:>10.1%}{A['p'][0]:>10.4g}")
        rows_out.append(dict(layer=lab.split()[0], beta_raw=b1, beta_adj=b2,
                             attenuation=a, p_raw=float(R["p"][0]),
                             p_adj=float(A["p"][0]), n_raw=R["n"], n_adj=A["n"]))

    print("\n" + "=" * 88)
    print("⑤ 免疫终点全表（8 个，FDR 校正；○=FDR<0.05 ·=0.05–0.10）")
    print("=" * 88)
    res = []
    for e in IMMUNE_SETS:
        if e not in S:
            continue
        R = ols(z(S[e]), z(Xaxis)[:, None], ["X"])
        A = ols(z(S[e]), np.column_stack([z(Xaxis), MES]), ["X", "MES"])
        res.append((e, float(R["beta"][0]), float(R["p"][0]),
                    float(A["beta"][0]), float(A["p"][0])))
    ps = np.array([r[2] for r in res])
    order = np.argsort(ps)
    fdr = np.empty_like(ps)
    for rank, i in enumerate(order):
        fdr[i] = min(ps[i] * len(ps) / (rank + 1), 1.0)
    for i, (e, b1, p1, b2, p2) in enumerate(res):
        star = "○" if fdr[i] < 0.05 else ("·" if fdr[i] < 0.10 else " ")
        tag = "   ← 主终点" if e == PRIMARY_ENDPOINT else ""
        print(f"  {star} {e:<16} β_raw={b1:+.3f} p={p1:.4g} | "
              f"β_adj={b2:+.3f} p={p2:.4g} | FDR={fdr[i]:.4g}{tag}")

    out = os.path.join(RES, "d4_branch_adjudication.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("item\tvalue\n")
        fh.write(f"primary_endpoint\t{PRIMARY_ENDPOINT}\n")
        fh.write(f"r_MES\t{r_main:.6f}\n")
        fh.write(f"r_MES_published\t{r_pub:.6f}\n")
        fh.write(f"r_threshold\t{R_THRESHOLD}\n")
        fh.write(f"beta_raw\t{br:.6f}\n")
        fh.write(f"beta_raw_p\t{pr:.6g}\n")
        fh.write(f"beta_adj\t{ba:.6f}\n")
        fh.write(f"beta_adj_p\t{pa:.6g}\n")
        fh.write(f"partial_r2\t{pr2:.6f}\n")
        fh.write(f"attenuation\t{att:.6f}\n")
        fh.write(f"attenuation_ci_lo\t{lo:.6f}\n")
        fh.write(f"attenuation_ci_hi\t{hi:.6f}\n")
        fh.write(f"estimable\t{int(estimable)}\n")
        fh.write(f"branch\t{branch}\n")
        fh.write(f"mes_d_cohens\t{d:.6f}\n")
        fh.write(f"n_raw\t{raw['n']}\n")
        fh.write(f"n_adj\t{adj['n']}\n")
        for r in rows_out:
            fh.write(f"sens_{r['layer']}_beta_raw\t{r['beta_raw']:.6f}\n")
            fh.write(f"sens_{r['layer']}_beta_adj\t{r['beta_adj']:.6f}\n")
            fh.write(f"sens_{r['layer']}_attenuation\t{r['attenuation']:.6f}\n")
            fh.write(f"sens_{r['layer']}_p_adj\t{r['p_adj']:.6g}\n")
    print(f"\n[out] {out}")


if __name__ == "__main__":
    main()
