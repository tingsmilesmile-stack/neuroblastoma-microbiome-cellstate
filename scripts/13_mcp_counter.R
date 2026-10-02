#!/usr/bin/env Rscript
# ---- circularity / purity analysis ----
suppressPackageStartupMessages({library(MCPcounter)})
BASE <- "."; RES <- file.path(BASE,"results")
m <- as.matrix(read.delim(file.path(RES,"expr_tpm.tsv"), row.names=1, check.names=FALSE))
# MCP-counter 用 HUGO 符号；去掉版本号、去重
rn <- sub("\\..*$", "", rownames(m)); keep <- !duplicated(rn)
m <- m[keep,]; rownames(m) <- rn[keep]
cat(sprintf("[info] %d genes x %d samples\n", nrow(m), ncol(m)))
genes <- read.table(file.path(BASE,"data/refs/mcpcounter_genes.txt"), sep="\t",
                    header=TRUE, stringsAsFactors=FALSE, colClasses="character", check.names=FALSE)
res <- MCPcounter.estimate(m, featuresType="HUGO_symbols", genes=genes)
write.table(data.frame(population=rownames(res), res, check.names=FALSE),
            file.path(RES,"main153_mcpcounter.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
axis <- read.delim(file.path(RES,"main153_mi_axis.tsv"))
stopifnot(identical(axis$sample, colnames(res)))
lab <- axis$k2
clin <- read.delim(file.path(RES,"sample_clinical_supp.tsv"), row.names=1, check.names=FALSE)
pur <- suppressWarnings(as.numeric(clin[axis$sample,"purity_pct"]))
rows <- list()
for (p in rownames(res)){
  v <- res[p,]
  r <- cor(v, axis$MI_axis_z)
  wt <- wilcox.test(v[lab=="MI_high"], v[lab=="MI_low"])
  # 纯度校正后的 MI 系数
  ok <- !is.na(pur)
  m1 <- summary(lm(scale(v)[,1] ~ scale(axis$MI_axis_z) + scale(pur)))$coefficients
  rows[[length(rows)+1]] <- data.frame(population=p, r_with_MIaxis=round(r,3),
    p_wilcox=signif(wt$p.value,3), beta_adj_purity=round(m1[2,1],3),
    p_adj_purity=signif(m1[2,4],3))
}
out <- do.call(rbind, rows); out$FDR <- signif(p.adjust(out$p_wilcox,"BH"),3)
write.table(out, file.path(RES,"main153_mcpcounter_vs_MIaxis.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(out, row.names=FALSE)
cat(sprintf("\n[summary] 与 MI 轴 |r|≥0.5 的群体: %d/%d；FDR<0.05: %d/%d\n",
    sum(abs(out$r_with_MIaxis)>=0.5), nrow(out), sum(out$FDR<0.05), nrow(out)))
cat("[done] 13_mcp_counter\n")
