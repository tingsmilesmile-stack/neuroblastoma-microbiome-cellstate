#!/usr/bin/env Rscript
# 17_deg_continuous.R — 连续炎症轴关联的差异表达（v2 主分析）
# 动机：v2 放弃中位数二分，主分析改连续；本脚本以轴 z 分数为连续自变量拟合 DESeq2。
# 设计：~ axis_z（连续），检验基因表达随轴单调变化的强度。
# 输出：results/v2_deg_continuous.tsv
suppressPackageStartupMessages({library(DESeq2)})
setwd("/Users/ting/Documents/数据/01肿瘤微生物")
RES <- "results"

cnt <- read.delim(file.path(RES, "main153_counts_unstranded.tsv"),
                  row.names = 1, check.names = FALSE)
# 第一列是 gene_name，需移出计数矩阵
gname <- cnt[[1]]
cnt <- cnt[, -1, drop = FALSE]
counts <- as.matrix(cnt)
storage.mode(counts) <- "integer"
cat(sprintf("[input] %d genes x %d samples\n", nrow(counts), ncol(counts)))

ax <- read.delim(file.path(RES, "main153_mi_axis.tsv"), stringsAsFactors = FALSE)
samples <- ax$sample
stopifnot(all(samples %in% colnames(counts)))
counts <- counts[, samples]
axis_z <- as.numeric(scale(ax$MI_axis_z))

coldata <- data.frame(row.names = samples, axis_z = axis_z)

# 独立过滤：rowSums >= 10（与 08_deg153.R 一致）
keep <- rowSums(counts) >= 10
cat(sprintf("[filter] rowSums>=10 : %d -> %d (removed %d)\n",
            nrow(counts), sum(keep), nrow(counts) - sum(keep)))
counts <- counts[keep, ]
gname <- gname[keep]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ axis_z)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, name = "axis_z", alpha = 0.05)
res$gene_name <- gname
out <- as.data.frame(res)
out <- data.frame(gene_id = rownames(out), gene_name = out$gene_name,
                  baseMean = out$baseMean, log2FoldChange = out$log2FoldChange,
                  lfcSE = out$lfcSE, stat = out$stat, pvalue = out$pvalue, padj = out$padj)
write.table(out, file.path(RES, "v2_deg_continuous.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("[done] tested=%d | padj<0.05 : %d\n",
            sum(!is.na(out$padj)), sum(out$padj < 0.05, na.rm = TRUE)))

s <- out[!is.na(out$padj), ]
cat(sprintf("  padj<0.05 & |log2FC|>0.5 : %d\n", sum(s$padj < 0.05 & abs(s$log2FoldChange) > 0.5)))
cat(sprintf("  padj<0.05 & |log2FC|>1.0 : %d\n", sum(s$padj < 0.05 & abs(s$log2FoldChange) > 1)))
# 分类计数（每 1 SD 轴变化）
up <- s[s$padj < 0.05 & s$log2FoldChange > 0.2, ]
dn <- s[s$padj < 0.05 & s$log2FoldChange < -0.2, ]
cat(sprintf("  per-SD |log2FC|>0.2 : up=%d down=%d\n", nrow(up), nrow(dn)))
cat("  top10 up:", paste(head(up[order(-up$log2FoldChange), "gene_name"], 10), collapse = ", "), "\n")
cat("  top10 down:", paste(head(dn[order(dn$log2FoldChange), "gene_name"], 10), collapse = ", "), "\n")
