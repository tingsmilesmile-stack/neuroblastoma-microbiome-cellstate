#!/usr/bin/env Rscript
# =============================================================================
# 03_gsva_ssgsea.R - standard ssGSEA scoring (GSVA) for all gene sets
# =============================================================================
# Rationale: pathway scores are computed with the standard GSVA ssGSEA
#   implementation (all parameters reported), not an in-house approximation.
# 本脚本：用 GSVA 2.6.6 的 ssgseaParam() 标准实现重算，报告全部参数。
#
# NOTE: the GSVA 2.x API differs from GSVA 1.x (verified experimentally):
#     旧写法 gsva(expr, geneSets, method="ssgsea", ...) 已废除；
#     新写法：先 ssgseaParam(exprData, geneSets, ...)，再 gsva(param)。
#
# 输入：results/expr_tpm.tsv（基因 × 样本，TPM，20,584 × 153）
# 输出：results/gsva_ssgsea_scores.tsv      （基因集 × 样本）
#       results/gsva_run_params.tsv          （运行参数与版本，供 Methods 直接引用）
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(jsonlite)
})

# NOTE: run this script from the repository root.
RES <- "results"

# ---- 基因集定义（必须与 02_pathways_ssgsea.py 完全一致，避免口径分叉）----
MICROBE_SETS <- list(
  TLR_signaling = c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9","MYD88","TIRAP",
                    "IRAK1","IRAK4","TRAF6","LY96","CD14","TICAM1"),
  NOD_like_receptor = c("NOD1","NOD2","NLRP3","NLRC4","NLRP1","PYCARD","CASP1",
                        "RIPK2","NAIP","AIM2"),
  Antimicrobial_peptides = c("DEFB1","DEFB4A","DEFA1","DEFA5","DEFA6","CAMP",
                             "LCN2","REG3A","REG3G","LTF","BPI","S100A8","S100A9"),
  LPS_inflammatory = c("TNF","IL1B","IL6","NFKB1","RELA","CXCL8","CCL2","PTGS2",
                       "NOS2","IL18","CXCL10"),
  SCFA_butyrate_response = c("FFAR2","FFAR3","HCAR2","SLC5A8","SLC16A1","HDAC1",
                             "HDAC3","ACADM"),
  Bile_acid_metabolism = c("NR1H4","NR1H3","CYP7A1","CYP8B1","SLC10A1","SLCO1B1",
                           "ABCB11","FABP6","NR0B2","GPBAR1"),   # TGR5 -> GPBAR1
  Tryptophan_IDO_kynurenine = c("IDO1","IDO2","TDO2","KYNU","KMO","AHR","ARNT",
                                "QPRT","HAAO","AFMID"),
  Mucin_barrier = c("MUC1","MUC2","MUC4","MUC5AC","TFF3","CLDN1","CLDN3","CLDN4",
                    "OCLN","TJP1")
)
IMMUNE_SETS <- list(
  CD8_T_cells = c("CD8A","CD8B","GZMK","GZMA","CD3D","CD3E"),
  Cytotoxic_NK = c("NKG7","KLRD1","KLRK1","GNLY","PRF1","GZMB","NCR1"),
  Tregs = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  Macrophage_M1 = c("NOS2","IL1B","TNF","CXCL9","CXCL10","CD80"),
  Macrophage_M2 = c("CD163","MRC1","MSR1","CD68","IL10","CCL22","ARG1"),
  B_cells = c("CD19","MS4A1","CD79A","CD79B","IGHM"),          # IGHM included
  DC = c("ITGAX","CD1C","CLEC9A","LILRA4","IRF8"),
  MHC_class_I = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2")
)
ALL_SETS <- c(MICROBE_SETS, IMMUNE_SETS)

# ---- 读入矩阵 ----
m <- as.matrix(read.delim(file.path(RES, "expr_tpm.tsv"), row.names = 1, check.names = FALSE))
cat(sprintf("[input] %d genes x %d samples\n", nrow(m), ncol(m)))

# ---- 基因集命中统计（每个基因集"实际参与评分 / 总数"）----
cov <- do.call(rbind, lapply(names(ALL_SETS), function(nm) {
  g <- ALL_SETS[[nm]]
  hit <- intersect(g, rownames(m))
  data.frame(gene_set = nm, n_defined = length(g), n_found = length(hit),
             missing = paste(setdiff(g, rownames(m)), collapse = ";"),
             stringsAsFactors = FALSE)
}))
write.table(cov, file.path(RES, "gsva_gene_coverage.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("[coverage] 基因集数:", nrow(cov), " 命中不足的集:",
    sum(cov$n_found < cov$n_defined), "\n")
print(cov[cov$n_found < cov$n_defined, c("gene_set","n_defined","n_found","missing")])

# ---- 标准 ssGSEA（GSVA 2.6.6 新 API）----
ALPHA <- 0.25; MIN_SIZE <- 1; MAX_SIZE <- Inf; NORMALIZE <- TRUE
param <- ssgseaParam(exprData = m, geneSets = ALL_SETS,
                     alpha = ALPHA, normalize = NORMALIZE,
                     minSize = MIN_SIZE, maxSize = MAX_SIZE)
t0 <- Sys.time()
sc <- gsva(param)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[gsva] 完成：%d 基因集 x %d 样本，耗时 %.1f 秒\n",
            nrow(sc), ncol(sc), elapsed))

# ---- 输出 ----
out <- data.frame(gene_set = rownames(sc), sc, check.names = FALSE)
write.table(out, file.path(RES, "gsva_ssgsea_scores.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

params <- data.frame(
  item = c("R_version","GSVA_version","Bioconductor_version","method",
           "alpha","normalize","minSize","maxSize","kcdf",
           "input_matrix","input_scale","n_genes","n_samples","runtime_sec",
           "gene_set_source","date"),
  value = c(paste(R.version$major, R.version$minor, sep = "."),
            as.character(packageVersion("GSVA")),
            as.character(BiocManager::version()),
            "ssGSEA (Barbie et al. 2009) via GSVA::ssgseaParam",
            as.character(ALPHA), as.character(NORMALIZE),
            as.character(MIN_SIZE), as.character(MAX_SIZE),
            "Gaussian (TPM, log2 not applied — GSVA handles ranks)",
            "results/expr_tpm.tsv", "TPM (unstranded)",
            as.character(nrow(m)), as.character(ncol(m)),
            sprintf("%.1f", elapsed),
            "curated (see gsva_gene_coverage.tsv)",
            format(Sys.time(), "%Y-%m-%d %H:%M")),
  stringsAsFactors = FALSE)
write.table(params, file.path(RES, "gsva_run_params.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("[done] 03_gsva_ssgsea.R\n")
print(params)
