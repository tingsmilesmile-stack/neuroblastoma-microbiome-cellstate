#!/usr/bin/env python3
"""D1 — 153 例「临床补充表 + GDC 临床」合并（严格解析）+ 协变量覆盖率重测

# Purpose: build the merged clinical/sample table for the analysis cohort.
--------------------------------
* D1 要把样本量从被截断的 100 例扩到全部 153 例原发瘤；
* 同时「临床补充数据清洗合并」（MYCN / INSS / COG / ploidy / Percent Tumor / EFS）
  必须在 153 例上重做；
* 旧覆盖率表基于 100 例，须在 153 例上**重测**：覆盖率 <80% 的协变量降级为敏感性分析。

数据来源
--------
1. 表达矩阵列名 = 153 例 submitter_id（唯一，无重复 aliquot）→ 作为**样本定标**，
   保证协变量与矩阵逐列对齐，不靠二次映射。
2. `data/clinical_supplement/*.xlsx`（由 `00b_fetch_clinical_supplement.py` 落盘，含 SHA256）：
   TARGET 临床补充表，**行稀疏**，
   ⚠️ 历史上曾放 `/tmp` 被 macOS 清理 → `lib_clinical` 静默返回 0 条(documented in the project troubleshooting notes).
   必须用 `lib_clinical.py` 按单元格引用归位（严禁 dict(zip(header,row)) 位置对齐）。
   字段：Age at Diagnosis in Days / INSS Stage / MYCN status / Ploidy /
        COG Risk Group / Percent Tumor / First Event / Event Free Survival Time in Days
3. `data/nbl_clinical.tsv`（GDC 官方临床导出）：OS 时间 / 生存状态 / 性别。

输出
----
* results/sample_clinical_supp.tsv    — 153 行 × 协变量宽表（S1/Table S1 基础）
* results/covariate_coverage_153.tsv  — 协变量覆盖率 + <80% 降级判定
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_clinical import load_clinical, get_field, COG_MAP, MYCN_MAP, INSS_MAP, PLOIDY_MAP

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
DATA = os.path.join(BASE, "data")

COVERAGE_FLOOR = 0.80   # MEMO：<80% 降级为敏感性分析


# ---------- 样本定标：以表达矩阵列名为准 ----------
def matrix_samples():
    with open(os.path.join(RES, "expr_tpm.tsv"), encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    return hdr[1:]        # 第一列是 gene_name


samples = matrix_samples()
print(f"[info] 表达矩阵样本数: {len(samples)}")
assert len(set(samples)) == len(samples), "矩阵列名存在重复，需先处理重复 aliquot"

# ---------- 临床补充表（严格解析）----------
CL = load_clinical()
print(f"[info] 临床补充表记录总数: {len(CL)}")

# 逐字段取：全部走 lib_clinical 的 get_field，禁止位置对齐
age_d,  ok_age  = get_field(CL, samples, "Age at Diagnosis in Days")
purity, ok_pur  = get_field(CL, samples, "Percent Tumor", numeric_loose=True)
cog,    ok_cog  = get_field(CL, samples, "COG Risk Group", COG_MAP)
mycn,   ok_mycn = get_field(CL, samples, "MYCN status", MYCN_MAP)
inss,   ok_inss = get_field(CL, samples, "INSS Stage", INSS_MAP)
ploidy, ok_plo  = get_field(CL, samples, "Ploidy", PLOIDY_MAP)
efs_t,  ok_efs  = get_field(CL, samples, "Event Free Survival Time in Days")

# `First Event` **不是二值**，实际为 6 类：
#   Censored / Relapse / Progression / Event / Death / Second Malignant Neoplasm
# EFS 事件 = 非 Censored（含 Relapse / Progression / Death / SMN 等全部首事件）
# ⚠️ 踩过的坑：曾把 "Relapse"/"Progression" 等误判为非事件，导致事件数被压到 16/73。
first_event = [str(CL.get(s, {}).get("First Event", "")).strip() for s in samples]
efs_event = np.array([0.0 if v.lower().startswith("censor") else 1.0
                      for v in first_event])
ok_efs_event = np.array([bool(v) for v in first_event])   # 153/153 均有首事件记录

# ---------- GDC 官方临床（OS / 性别）----------
gdc = pd.read_csv(os.path.join(DATA, "nbl_clinical.tsv"), sep="\t")
gdc = gdc.rename(columns={
    "demographic.days_to_death": "days_to_death",
    "demographic.gender": "gender",
    "demographic.vital_status": "vital_status",
    "diagnoses.0.age_at_diagnosis": "gdc_age_days",
    "diagnoses.0.days_to_last_follow_up": "days_to_last_follow_up",
})


def os_time(r):
    if r["vital_status"] == "Dead" and not pd.isna(r["days_to_death"]):
        return r["days_to_death"]
    return r["days_to_last_follow_up"]


gdc["OS_time"] = gdc.apply(os_time, axis=1)
gdc["OS_event"] = (gdc["vital_status"] == "Dead").astype(int)
gdc = gdc[["submitter_id", "gender", "vital_status", "OS_time", "OS_event",
           "gdc_age_days"]]

out = pd.DataFrame({
    "sample": samples,
    "submitter_id": samples,
    "age_days": age_d,
    "age_years": age_d / 365.25,
    "cog_risk": cog,          # 0=Low 1=Intermediate 2=High
    "mycn_amp": mycn,         # 1=Amplified 0=Not
    "inss_stage": inss,
    "ploidy": ploidy,         # 1=Hyperdiploid 0=Diploid
    "purity_pct": purity,
    "efs_time": efs_t,
    "efs_event": efs_event,
    "efs_first_event": first_event,
})
out = out.merge(gdc, on="submitter_id", how="left")
out.to_csv(os.path.join(RES, "sample_clinical_supp.tsv"), sep="\t", index=False)
print(f"[out] results/sample_clinical_supp.tsv  {out.shape[0]} 行 × {out.shape[1]} 列")

# ---------- 协变量覆盖率重测（153 例口径）----------
cov = {
    "age_years": ok_age,
    "age_days": ok_age,
    "cog_risk": ok_cog,
    "mycn_amp": ok_mycn,
    "inss_stage": ok_inss,
    "ploidy": ok_plo,
    "purity_pct": ok_pur,
    "efs_time": ok_efs,
    "efs_event": ok_efs_event,
    "os_time": out["OS_time"].notna().to_numpy(),
    "os_event": out["OS_event"].notna().to_numpy(),
    "gender": out["gender"].notna().to_numpy(),
}
rows = []
for k, m in cov.items():
    n = int(np.asarray(m).sum())
    rate = n / len(samples)
    rows.append({
        "covariate": k,
        "n_available": n,
        "n_total": len(samples),
        "coverage": round(rate, 4),
        "verdict": "OK" if rate >= COVERAGE_FLOOR else "降级为敏感性(<80%)",
    })
cvd = pd.DataFrame(rows).sort_values("coverage", ascending=False)
cvd.to_csv(os.path.join(RES, "covariate_coverage_153.tsv"), sep="\t", index=False)

print(f"\n[覆盖率] 153 例口径（阈值 {COVERAGE_FLOOR:.0%}）")
for _, r in cvd.iterrows():
    print(f"  {r['covariate']:<12} {r['n_available']:>4}/{r['n_total']}  "
          f"{r['coverage']:.1%}  {r['verdict']}")

# 事件数（用于建模预算：每协变量 >=10 事件）
print(f"\n[事件数] OS 事件 = {int(out['OS_event'].sum())} / 153"
      f"    EFS 事件 = {int(np.nansum(out['efs_event']))} / "
      f"{int(ok_efs_event.sum())}")
print("[done] D1_merge_supp_clinical")
