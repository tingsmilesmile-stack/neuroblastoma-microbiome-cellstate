#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_supplementary_v2.py
Rebuilds the entire suite of Supplementary Tables for the 153-sample cohort.
Outputs:
- submission/Supplementary_Tables.xlsx (Sheets: S1_Sample_Flow, S2_MI_Pathway_Gene_Sets, S3_Alpha_Sensitivity, S4_Pathway_OS_Cox, S5_Immune_Associations, S6_Differential_Expression, S7_GDC_Sample_IDs, S8_Immune_Signatures)
- Individual CSVs in submission/ to maintain consistency with the old file structure.
"""

import os
import json
import pandas as pd
import openpyxl

# Set paths
RES = "results"
SUB = "submission"
os.makedirs(SUB, exist_ok=True)

# ----------------- Sheet S1: Sample Flow -----------------
print("Building S1...")
df_s1 = pd.read_csv(os.path.join(RES, "main153_sample_flow.tsv"), sep="\t")

# ----------------- Sheet S2: MI Pathway Gene Sets -----------------
print("Building S2...")
# Gene set definitions and rationale from 03_gsva_ssgsea.R
MICROBE_SETS = {
    "TLR_signaling": ("TLR1, TLR2, TLR4, TLR5, TLR6, TLR9, MYD88, TIRAP, IRAK1, IRAK4, TRAF6, LY96, CD14, TICAM1", "TLR activation by bacterial LPS, lipoteichoic acid, flagellin, and CpG DNA"),
    "NOD_like_receptor": ("NOD1, NOD2, NLRP3, NLRC4, NLRP1, PYCARD, CASP1, RIPK2, NAIP, AIM2", "Intracellular sensing of bacterial peptidoglycan and inflammasome assembly"),
    "Antimicrobial_peptides": ("DEFB1, DEFB4A, DEFA1, DEFA5, DEFA6, CAMP, LCN2, REG3A, REG3G, LTF, BPI, S100A8, S100A9", "Innate host bactericidal and barrier-defense effector proteins"),
    "LPS_inflammatory": ("TNF, IL1B, IL6, NFKB1, RELA, CXCL8, CCL2, PTGS2, NOS2, IL18, CXCL10", "Host-cell transcriptional inflammatory response cascade triggered by LPS"),
    "SCFA_butyrate_response": ("FFAR2, FFAR3, HCAR2, SLC5A8, SLC16A1, HDAC1, HDAC3, ACADM", "Sensing and response to bacterial short-chain fatty acids (SCFAs)"),
    "Bile_acid_metabolism": ("NR1H4, NR1H3, CYP7A1, CYP8B1, SLC10A1, SLCO1B1, ABCB11, FABP6, NR0B2, GPBAR1", "Sensing and feedback regulation of host and microbial bile acids"),
    "Tryptophan_IDO_kynurenine": ("IDO1, IDO2, TDO2, KYNU, KMO, AHR, ARNT, QPRT, HAAO, AFMID", "Tryptophan catabolism to immunosuppressive kynurenines"),
    "Mucin_barrier": ("MUC1, MUC2, MUC4, MUC5AC, TFF3, CLDN1, CLDN3, CLDN4, OCLN, TJP1", "Epithelial mucous layer and tight-junction physical barrier proteins")
}

rows_s2 = []
cov_df = pd.read_csv(os.path.join(RES, "gsva_gene_coverage.tsv"), sep="\t")
for name, (genes, rationale) in MICROBE_SETS.items():
    cov_row = cov_df[cov_df["gene_set"] == name]
    found = int(cov_row["n_found"].iloc[0]) if not cov_row.empty else 0
    defined = int(cov_row["n_defined"].iloc[0]) if not cov_row.empty else 0
    missing = cov_row["missing"].iloc[0] if not cov_row.empty else ""
    if pd.isna(missing) or not missing:
        missing = "None (100% coverage)"
    rows_s2.append({
        "Pathway Name": name.replace("_", " "),
        "Defined Genes": genes,
        "Biological Rationale": rationale,
        "Genes Present in Matrix": f"{found} / {defined}",
        "Missing Genes in TARGET": missing
    })
df_s2 = pd.DataFrame(rows_s2)

# ----------------- Sheet S3: Alpha Sensitivity -----------------
print("Building S3...")
df_s3 = pd.read_csv(os.path.join(RES, "main153_alpha_sensitivity.tsv"), sep="\t")
df_s3.columns = ["ssGSEA alpha parameter", "r (Inflammatory axis, MHC class I)", "r (SCFA, Inflammatory axis)", "r (Inflammatory axis, ADRN/MES)"]

# ----------------- Sheet S4: Pathway OS Cox -----------------
print("Building S4...")
df_s4 = pd.read_csv(os.path.join(RES, "main153_pathway_cox.tsv"), sep="\t")

# ----------------- Sheet S5: Immune Associations -----------------
print("Building S5...")
df_s5 = pd.read_csv(os.path.join(RES, "main153_immune_assoc.tsv"), sep="\t")

# ----------------- Sheet S6: Differential Expression (top 2000) -----------------
print("Building S6...")
df_deg_all = pd.read_csv(os.path.join(RES, "main153_deg_deseq2.tsv"), sep="\t")
df_deg_sig = df_deg_all[(df_deg_all["padj"] < 0.05) & (df_deg_all["log2FoldChange"].abs() > 1)].copy()
df_deg_sig["abs_lfc"] = df_deg_sig["log2FoldChange"].abs()
df_deg_sig = df_deg_sig.sort_values(by=["padj", "abs_lfc"], ascending=[True, False]).drop(columns=["abs_lfc"])
df_s6 = df_deg_sig.head(2000) # Keep top 2000 for supplementary sheets (FDR sorted)

# ----------------- Sheet S7: GDC Sample IDs -----------------
print("Building S7...")
# Load primary_expr_ids.json
with open("data/primary_expr_ids.json") as f:
    gdc_list = json.load(f) # list of list [file_id, file_name, case_id]
# We also have main153_mi_axis.tsv for cluster assignments
df_ax = pd.read_csv(os.path.join(RES, "main153_mi_axis.tsv"), sep="\t")
df_gdc = pd.DataFrame(gdc_list, columns=["GDC_File_ID", "GDC_File_Name", "TARGET_Submitter_ID"])
# Merge to add subtype assignment
df_s7 = pd.merge(df_gdc, df_ax, left_on="TARGET_Submitter_ID", right_on="sample", how="inner").drop(columns=["sample"])
df_s7.columns = ["GDC File ID", "GDC File Name", "TARGET Submitter ID", "Subtype Group Assignment (Median Split)", "Continuous Inflammatory-axis Score (z)"]

# ----------------- Sheet S8: Immune Signatures -----------------
print("Building S8...")
IMMUNE_SETS = {
    "CD8_T_cells": "CD8A, CD8B, GZMK, GZMA, CD3D, CD3E",
    "Cytotoxic_NK": "NKG7, KLRD1, KLRK1, GNLY, PRF1, GZMB, NCR1",
    "Tregs": "FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18",
    "Macrophage_M1": "NOS2, IL1B, TNF, CXCL9, CXCL10, CD80",
    "Macrophage_M2": "CD163, MRC1, MSR1, CD68, IL10, CCL22, ARG1",
    "B_cells": "CD19, MS4A1, CD79A, CD79B, IGHM",
    "DC": "ITGAX, CD1C, CLEC9A, LILRA4, IRF8",
    "MHC_class_I": "HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2"
}
rows_s8 = []
for name, genes in IMMUNE_SETS.items():
    cov_row = cov_df[cov_df["gene_set"] == name]
    found = int(cov_row["n_found"].iloc[0]) if not cov_row.empty else 0
    defined = int(cov_row["n_defined"].iloc[0]) if not cov_row.empty else 0
    missing = cov_row["missing"].iloc[0] if not cov_row.empty else ""
    if pd.isna(missing) or not missing:
        missing = "None (100% coverage)"
    rows_s8.append({
        "Signature Name": name.replace("_", " "),
        "Marker Genes": genes,
        "Functional Role": "Antigen presentation" if name == "MHC_class_I" else f"Marker-based cellular deconvolution signature for {name.replace('_', ' ')}",
        "Genes Present in Matrix": f"{found} / {defined}",
        "Missing Genes in TARGET": missing
    })
df_s8 = pd.DataFrame(rows_s8)

# ----------------- Write unified xlsx -----------------
xlsx_path = os.path.join(SUB, "Supplementary_Tables.xlsx")
print(f"Writing unified XLSX to {xlsx_path}...")
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    df_s1.to_excel(writer, sheet_name="S1_Sample_Flow", index=False)
    df_s2.to_excel(writer, sheet_name="S2_MI_Pathway_Gene_Sets", index=False)
    df_s3.to_excel(writer, sheet_name="S3_Alpha_Sensitivity", index=False)
    df_s4.to_excel(writer, sheet_name="S4_Pathway_OS_Cox", index=False)
    df_s5.to_excel(writer, sheet_name="S5_Immune_Associations", index=False)
    df_s6.to_excel(writer, sheet_name="S6_Differential_Expression", index=False)
    df_s7.to_excel(writer, sheet_name="S7_GDC_Sample_IDs", index=False)
    df_s8.to_excel(writer, sheet_name="S8_Immune_Signatures", index=False)

# ----------------- Write individual CSVs (old compatibility) -----------------
# 1. S1_sample_subtype_clinical.csv (153 clinical data with axis score)
print("Writing S1_sample_subtype_clinical.csv...")
df_clin_raw = pd.read_csv(os.path.join(RES, "sample_clinical_supp.tsv"), sep="\t")
df_s1_csv = pd.merge(df_clin_raw, df_ax, on="sample", how="inner")
df_s1_csv.to_csv(os.path.join(SUB, "S1_sample_subtype_clinical.csv"), index=False)

# 2. S2_pathway_ssGSEA.csv (raw ssGSEA scores of 8 pathways for 153 samples)
print("Writing S2_pathway_ssGSEA.csv...")
s_raw = pd.read_csv(os.path.join(RES, "gsva_ssgsea_scores.tsv"), sep="\t")
# filter to 8 MI pathways
s_raw_8 = s_raw[s_raw["gene_set"].isin(MICROBE_SETS.keys())]
s_raw_8.to_csv(os.path.join(SUB, "S2_pathway_ssGSEA.csv"), index=False)

# 3. S3_KM_summary.csv
print("Writing S3_KM_summary.csv...")
df_s3_csv = pd.read_csv(os.path.join(RES, "main153_km_summary.tsv"), sep="\t")
df_s3_csv.to_csv(os.path.join(SUB, "S3_KM_summary.csv"), index=False)

# 4. S4_pathway_Cox.csv
print("Writing S4_pathway_Cox.csv...")
df_s4.to_csv(os.path.join(SUB, "S4_pathway_Cox.csv"), index=False)

# 5. S5_immune_tests.csv
print("Writing S5_immune_tests.csv...")
df_s5.to_csv(os.path.join(SUB, "S5_immune_tests.csv"), index=False)

# 6. S6_significant_DEG.csv (full FDR<0.05, |log2FC|>1 list)
print("Writing S6_significant_DEG.csv...")
df_deg_sig.to_csv(os.path.join(SUB, "S6_significant_DEG.csv"), index=False)

print("SUCCESS: All supplementary tables built!")
