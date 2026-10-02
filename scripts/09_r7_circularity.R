#!/usr/bin/env Rscript
# ---- circularity / purity analysis ----
# ③ MI 轴内部去重（leave-one-set-out）；④ 与免疫集的重叠明细。
suppressPackageStartupMessages({library(GSVA)})
BASE <- "."; RES <- file.path(BASE,"results")
setwd(BASE)
MIC <- list(
  TLR_signaling=c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9","MYD88","TIRAP","IRAK1","IRAK4","TRAF6","LY96","CD14","TICAM1"),
  NOD_like_receptor=c("NOD1","NOD2","NLRP3","NLRC4","NLRP1","PYCARD","CASP1","RIPK2","NAIP","AIM2"),
  Antimicrobial_peptides=c("DEFB1","DEFB4A","DEFA1","DEFA5","DEFA6","CAMP","LCN2","REG3A","REG3G","LTF","BPI","S100A8","S100A9"),
  LPS_inflammatory=c("TNF","IL1B","IL6","NFKB1","RELA","CXCL8","CCL2","PTGS2","NOS2","IL18","CXCL10"),
  SCFA_butyrate_response=c("FFAR2","FFAR3","HCAR2","SLC5A8","SLC16A1","HDAC1","HDAC3","ACADM"),
  Bile_acid_metabolism=c("NR1H4","NR1H3","CYP7A1","CYP8B1","SLC10A1","SLCO1B1","ABCB11","FABP6","NR0B2","GPBAR1"),
  Tryptophan_IDO_kynurenine=c("IDO1","IDO2","TDO2","KYNU","KMO","AHR","ARNT","QPRT","HAAO","AFMID"),
  Mucin_barrier=c("MUC1","MUC2","MUC4","MUC5AC","TFF3","CLDN1","CLDN3","CLDN4","OCLN","TJP1"))
IMM <- list(
  CD8_T_cells=c("CD8A","CD8B","GZMK","GZMA","CD3D","CD3E"),
  Cytotoxic_NK=c("NKG7","KLRD1","KLRK1","GNLY","PRF1","GZMB","NCR1"),
  Tregs=c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  Macrophage_M1=c("NOS2","IL1B","TNF","CXCL9","CXCL10","CD80"),
  Macrophage_M2=c("CD163","MRC1","MSR1","CD68","IL10","CCL22","ARG1"),
  B_cells=c("CD19","MS4A1","CD79A","CD79B","IGHM"),
  DC=c("ITGAX","CD1C","CLEC9A","LILRA4","IRF8"),
  MHC_class_I=c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2"))
INFL <- c("TLR_signaling","NOD_like_receptor","LPS_inflammatory")

m <- as.matrix(read.delim(file.path(RES,"expr_tpm.tsv"), row.names=1, check.names=FALSE))
score <- function(sets){ p <- ssgseaParam(exprData=m, geneSets=sets, alpha=0.25, normalize=TRUE, minSize=1)
  s <- gsva(p); t(scale(t(s))) }        # 每集 z

# ① 直接重叠
ov <- data.frame(microbe_gene=unlist(MIC))
ov2 <- data.frame(immune_gene=unlist(IMM))
overlap <- intersect(unlist(MIC), unlist(IMM))
cat("① MI 集 ∩ 免疫集 直接重叠基因:", paste(overlap, collapse=", "), "｜共", length(overlap), "\n")
cat("   占免疫集比例:", sprintf("%d/%d = %.1f%%", length(overlap), length(unique(unlist(IMM))),
    100*length(overlap)/length(unique(unlist(IMM)))), "\n")

Z0 <- score(c(MIC, IMM))
axis0 <- colMeans(Z0[INFL,,drop=FALSE]); mhc0 <- Z0["MHC_class_I",]
r0 <- cor(axis0, mhc0)
cat(sprintf("   基线 r(MI 炎症轴, MHC-I) = %.3f\n", r0))

# ② 剔除 4 个重叠基因后重算
strip <- function(L, drop) lapply(L, function(g) setdiff(g, drop))
M2 <- strip(MIC, overlap); I2 <- strip(IMM, overlap)
Z1 <- score(c(M2, I2))
axis1 <- colMeans(Z1[INFL,,drop=FALSE]); mhc1 <- Z1["MHC_class_I",]
r1 <- cor(axis1, mhc1)
cat(sprintf("② 剔除重叠基因后 r(MI 炎症轴, MHC-I) = %.3f（Δ=%.4f）\n", r1, r1-r0))
cat("   免疫集剔除后剩余基因数:", sapply(I2, length), "\n")

# ③ leave-one-set-out：逐一去掉炎症轴的一个集
loo <- sapply(INFL, function(g){ zz <- score(c(MIC[setdiff(INFL,g)], IMM))
  cor(colMeans(zz[setdiff(INFL,g),,drop=FALSE]), zz["MHC_class_I",]) })
cat("③ 炎症轴 leave-one-set-out r:\n"); print(round(loo,3))

# ④ 每个免疫签名与 MI 轴的 r（基线 & 剔除后）
res <- data.frame(
  signature=names(IMM),
  r_with_MI_axis = round(sapply(names(IMM), function(s) cor(axis0, Z0[s,])),3),
  r_with_MI_axis_excl_overlap = round(sapply(names(IMM), function(s) cor(axis1, Z1[s,])),3))
write.table(res, file.path(RES,"main153_overlap_exclusion.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(res)
cat(sprintf("\n[summary] 基线 r=%.3f｜剔除后 r=%.3f（|Δ|=%.4f）｜重叠基因 %d 个\n",
    r0, r1, abs(r1-r0), length(overlap)))
cat("[done] 09_r7_circularity\n")
