#!/usr/bin/env Rscript
# 12_alpha_sensitivity.R — 冻结要求：凡报 ssGSEA 关联须同时给 alpha ∈ {0.25,0.50,1.00}
suppressPackageStartupMessages(library(GSVA))
BASE <- "."; RES <- file.path(BASE,"results")
MIC <- list(
  TLR_signaling=c("TLR1","TLR2","TLR4","TLR5","TLR6","TLR9","MYD88","TIRAP","IRAK1","IRAK4","TRAF6","LY96","CD14","TICAM1"),
  NOD_like_receptor=c("NOD1","NOD2","NLRP3","NLRC4","NLRP1","PYCARD","CASP1","RIPK2","NAIP","AIM2"),
  LPS_inflammatory=c("TNF","IL1B","IL6","NFKB1","RELA","CXCL8","CCL2","PTGS2","NOS2","IL18","CXCL10"),
  SCFA_butyrate_response=c("FFAR2","FFAR3","HCAR2","SLC5A8","SLC16A1","HDAC1","HDAC3","ACADM"))
MHC <- list(MHC_class_I=c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2"))
m <- as.matrix(read.delim(file.path(RES,"expr_tpm.tsv"), row.names=1, check.names=FALSE))
set <- c(MIC, MHC)
rows <- list()
for (a in c(0.25, 0.50, 1.00)){
  sc <- gsva(ssgseaParam(exprData=m, geneSets=set, alpha=a, normalize=TRUE, minSize=1))
  z <- t(scale(t(sc)))
  infl <- colMeans(z[c("TLR_signaling","NOD_like_receptor","LPS_inflammatory"),,drop=FALSE])
  r_imm  <- cor(infl, z["MHC_class_I",])
  r_scfa <- cor(z["SCFA_butyrate_response",], infl)
  # 与正式 ADRN/MES 签名
  pub <- read.delim(file.path(RES,"gsva_published_adrn_mes_scores.tsv"), row.names=1, check.names=FALSE)
  pz <- t(scale(t(as.matrix(pub))))
  mes <- pz["VANGRONINGEN_ADRN",] - pz["VANGRONINGEN_MES",]
  r_mes <- cor(infl, mes)
  rows[[length(rows)+1]] <- data.frame(alpha=a, r_MIaxis_MHCclassI=round(r_imm,3),
    r_SCFA_inflaxis=round(r_scfa,3), r_MIaxis_ADRNMES=round(r_mes,3))
  cat(sprintf("[alpha=%.2f] r(MI轴,MHC-I)=%.3f  r(SCFA,炎症轴)=%+.3f  r(MI轴,ADRN-MES)=%+.3f\n",
      a, r_imm, r_scfa, r_mes))
}
out <- do.call(rbind, rows)
write.table(out, file.path(RES,"main153_alpha_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("[done] 12_alpha_sensitivity\n")
