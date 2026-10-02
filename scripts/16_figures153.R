#!/usr/bin/env Rscript
# 16_figures153.R -- 重制 6 主图 + 2 补充图（153 例队列）
# 叙事：连续炎症轴（禁止 "subtype"/"MI-high"/"MI-low" 出现在图中）
# 输出：figures/v2/*.pdf + *.png(300dpi)
suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap); library(survival); library(survminer)
  library(ggrepel); library(patchwork); library(RColorBrewer); library(grid)
})
set.seed(42)
OUT <- "figures/v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
HI <- "#D7263D"; LO <- "#1B6CA8"; MID <- "#F2F2F2"
MI_SETS <- c("TLR_signaling","NOD_like_receptor","Antimicrobial_peptides","LPS_inflammatory",
             "SCFA_butyrate_response","Bile_acid_metabolism","Tryptophan_IDO_kynurenine","Mucin_barrier")
IMM_SETS <- c("CD8_T_cells","Cytotoxic_NK","Tregs","Macrophage_M1","Macrophage_M2","B_cells","DC","MHC_class_I")
PRETTY <- c(TLR_signaling="TLR signaling", NOD_like_receptor="NOD-like receptor",
            Antimicrobial_peptides="Antimicrobial peptides", LPS_inflammatory="LPS/inflammatory",
            SCFA_butyrate_response="SCFA/butyrate", Bile_acid_metabolism="Bile-acid metabolism",
            Tryptophan_IDO_kynurenine="Tryptophan/IDO-kynurenine", Mucin_barrier="Mucin barrier",
            CD8_T_cells="CD8+ T cells", Cytotoxic_NK="Cytotoxic NK", Tregs="Regulatory T cells",
            Macrophage_M1="M1 macrophages", Macrophage_M2="M2 macrophages", B_cells="B cells",
            DC="Dendritic cells", MHC_class_I="MHC class I (antigen-presentation/expression signature)")
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

## ---------- 数据 ----------
s   <- read.delim("results/gsva_ssgsea_scores.tsv", row.names = 1, check.names = FALSE)
ax  <- read.delim("results/main153_mi_axis.tsv", stringsAsFactors = FALSE)  # sample,k2,MI_axis_z
samples <- ax$sample
stopifnot(all(samples %in% colnames(s)))
Z <- as.matrix(s[c(MI_SETS, IMM_SETS), samples])
ax$k2lab <- ifelse(ax$k2 == "MI_high", "Higher inflammatory-axis activity", "Lower inflammatory-axis activity")
ax$k2lab <- factor(ax$k2lab, levels = c("Lower inflammatory-axis activity", "Higher inflammatory-axis activity"))

## ---------- Fig1: 通路热图，样本按连续炎症轴排序 ----------
X <- t(scale(t(Z[MI_SETS, ])))
ord <- order(ax$MI_axis_z)
ann <- data.frame(`Inflammatory-axis activity` = ax$MI_axis_z[ord], check.names = FALSE,
                  row.names = samples[ord])
pal <- colorRampPalette(c(LO, MID, HI))(100)
pheatmap(X[, ord], color = pal, breaks = seq(-2, 2, length.out = 101), cluster_cols = FALSE,
         cluster_rows = TRUE, annotation_col = ann, show_colnames = FALSE, border_color = NA,
         labels_row = unname(PRETTY[MI_SETS]),
         annotation_colors = list(`Inflammatory-axis activity` = c(LO, MID, HI)),
         main = "Microbiome–host interaction pathway activity (ssGSEA z-score)",
         fontsize_row = 9, filename = file.path(OUT, "Fig1_pathway_heatmap.pdf"), width = 6, height = 3.2)
pheatmap(X[, ord], color = pal, breaks = seq(-2, 2, length.out = 101), cluster_cols = FALSE,
         cluster_rows = TRUE, annotation_col = ann, show_colnames = FALSE, border_color = NA,
         labels_row = unname(PRETTY[MI_SETS]),
         annotation_colors = list(`Inflammatory-axis activity` = c(LO, MID, HI)),
         main = "Microbiome–host interaction pathway activity (ssGSEA z-score)",
         fontsize_row = 9, filename = file.path(OUT, "Fig1_pathway_heatmap.png"), width = 6, height = 3.2, res = 300)

## ---------- Fig2: PCA 连续着色 ----------
pca <- prcomp(t(X), center = TRUE, scale. = FALSE)
pv  <- round(100 * summary(pca)$importance[2, 1:2], 1)
d2  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], axis = ax$MI_axis_z)
r_pc1 <- cor(d2$PC1, d2$axis)
p2 <- ggplot(d2, aes(PC1, PC2, colour = axis)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_colour_gradient2(low = LO, mid = MID, high = HI, midpoint = 0,
                         name = "Inflammatory-axis\nactivity (z)") +
  labs(title = "Principal component analysis of pathway activity",
       subtitle = sprintf("PC1 = %.1f%% of variance; r(PC1, inflammatory axis) = %.3f; n = %d",
                          pv[1], r_pc1, nrow(d2)),
       x = sprintf("PC1 (%.1f%% variance)", pv[1]), y = sprintf("PC2 (%.1f%% variance)", pv[2])) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save2(p2, "Fig2_PCA", 5.4, 4.2)

## ---------- Fig3: KM + 删失 + risk table ----------
cl <- read.delim("results/sample_clinical_supp.tsv", check.names = FALSE)
cl <- cl[cl$sample %in% samples, c("sample", "OS_time", "OS_event")]
cl <- merge(cl, ax[, c("sample", "k2lab")], by = "sample")
cl <- cl[is.finite(cl$OS_time) & is.finite(cl$OS_event), ]
fit <- survfit(Surv(OS_time, OS_event) ~ k2lab, data = cl)
lr  <- survdiff(Surv(OS_time, OS_event) ~ k2lab, data = cl)
p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
g3 <- ggsurvplot(fit, data = cl, pval = sprintf("Log-rank p = %.3f", p_lr), pval.coord = c(0.5, 0.15),
                 risk.table = TRUE, censor = TRUE, conf.int = TRUE, palette = c(LO, HI),
                 legend.title = "", legend.labs = levels(cl$k2lab),
                 xlab = "Time (years)", ylab = "Overall survival probability",
                 title = "Overall survival by inflammatory-axis activity",
                 ggtheme = theme_bw(base_size = 11), tables.height = 0.28)
cairo_pdf(file.path(OUT, "Fig3_KM.pdf"), width = 5.6, height = 5.4); print(g3); dev.off()
png(file.path(OUT, "Fig3_KM.png"), width = 5.6, height = 5.4, units = "in", res = 300); print(g3); dev.off()

## ---------- Fig4: 逐通路 Cox 森林图（p + q） ----------
cox <- read.delim("results/main153_pathway_cox.tsv", stringsAsFactors = FALSE)
cox <- cox[cox$pathway %in% MI_SETS, ]
cox$lab <- factor(unname(PRETTY[cox$pathway]), levels = rev(unname(PRETTY[MI_SETS])))
cox$sig <- ifelse(cox$q_BH < 0.05, "q<0.05", "ns")
p4 <- ggplot(cox, aes(HR, lab)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.18, colour = "grey30") +
  geom_point(aes(colour = sig), size = 2.6) +
  scale_colour_manual(values = c("q<0.05" = HI, "ns" = "grey45"), name = NULL) +
  geom_text(aes(x = max(cox$CI_hi) * 1.02, label = sprintf("p=%.3f, q=%.2f", p, q_BH)),
            hjust = 0, size = 3, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(title = "Per-pathway continuous Cox regression for overall survival",
       subtitle = "Hazard ratio per 1 SD of pathway activity (95% CI)", x = "Hazard ratio (95% CI)", y = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save2(p4, "Fig4_pathway_forest", 6.4, 3.4)

## ---------- Fig5: 免疫签名箱线图 ----------
imm <- read.delim("results/main153_immune_assoc.tsv", stringsAsFactors = FALSE)
long <- do.call(rbind, lapply(IMM_SETS, function(g)
  data.frame(sig = g, score = Z[g, ], grp = ax$k2lab)))
long$sig <- factor(unname(PRETTY[long$sig]), levels = unname(PRETTY[IMM_SETS]))
ann5 <- setNames(sprintf("%s\nKruskal p = %.3g; FDR = %.3g", unname(PRETTY[imm$signature]), imm$kruskal_p, imm$FDR),
                 unname(PRETTY[imm$signature]))
p5 <- ggplot(long, aes(grp, score, fill = grp)) +
  geom_boxplot(outlier.size = 0.5, width = 0.62) +
  facet_wrap(~sig, ncol = 4, scales = "free_y", labeller = labeller(sig = ann5)) +
  scale_fill_manual(values = c(LO, HI), guide = "none") +
  labs(title = "Immune and antigen-presentation signatures by inflammatory-axis activity",
       subtitle = "ssGSEA z-scores; boxes = median and interquartile range, whiskers = 1.5×IQR; open points = outliers",
       x = NULL, y = "ssGSEA z-score") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, size = 8), strip.text = element_text(size = 7.5))
save2(p5, "Fig5_immune_box", 9.2, 5.4)

## ---------- Fig6: 火山图 ----------
deg <- read.delim("results/main153_deg_deseq2.tsv", stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$padj) & is.finite(deg$log2FoldChange), ]
deg$q <- deg$padj
deg$sig <- ifelse(deg$q < 0.05 & abs(deg$log2FoldChange) > 1, "FDR<0.05 & |log2FC|>1", "ns")
top <- deg[deg$sig != "ns", ]
top <- rbind(head(top[order(-top$log2FoldChange), ], 12), head(top[order(top$log2FoldChange), ], 12))
p6 <- ggplot(deg, aes(log2FoldChange, -log10(q), colour = sig)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_colour_manual(values = c("FDR<0.05 & |log2FC|>1" = HI, "ns" = "grey75"), name = NULL) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_text_repel(data = top, aes(label = gene_name), size = 2.6, max.overlaps = 30, segment.size = 0.2) +
  labs(title = "Differential expression between higher and lower inflammatory-axis activity",
       subtitle = "DESeq2 on raw counts (unstranded)", x = "log2 fold-change", y = expression(-log[10]("adjusted p-value ("*italic(q)*")"))) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save2(p6, "Fig6_volcano", 6.2, 5.0)

## ---------- FigS1: 样本流 ----------
fl <- read.delim("results/main153_sample_flow.tsv", stringsAsFactors = FALSE)
fl$step <- factor(fl$step, levels = rev(fl$step))
pS1 <- ggplot(fl, aes(n, step, fill = grepl("最终|Primary Tumor 唯一|成功下载|临床补充", step))) +
  geom_col(width = 0.62) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c("TRUE" = HI, "FALSE" = "grey60"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Sample flow: TARGET-NBL primary-tumour RNA-seq cohort",
       x = "Number of files / patients", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
save2(pS1, "FigS1_sample_flow", 7.0, 4.0)

## ---------- FigS2: 共识聚类诊断 ----------
mats <- lapply(2:6, function(k) as.matrix(read.delim(sprintf("results/main153_consensus_k%d.tsv", k), row.names = 1, check.names = FALSE)))
hm <- lapply(seq_along(mats), function(i) {
  m <- mats[[i]]; o <- hclust(as.dist(1 - m), method = "ward.D2")$order
  df <- data.frame(x = rep(seq_len(ncol(m)), each = nrow(m)), y = rep(seq_len(nrow(m)), ncol(m)), v = as.vector(m))
  ggplot(df, aes(x, y, fill = v)) + geom_raster() +
    scale_fill_gradientn(colors = c(MID, HI), limits = c(0, 1), name = "consensus") +
    coord_fixed() + labs(title = paste0("k = ", i + 1), x = NULL, y = NULL) +
    theme_void(base_size = 9) + theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0.25, "cm"))
})
st <- read.delim("results/main153_stability.tsv", stringsAsFactors = FALSE)
pPAC <- ggplot(st, aes(k, PAC)) + geom_line(colour = HI) + geom_point(size = 2, colour = HI) +
  scale_x_continuous(breaks = 2:6) + labs(title = "PAC vs k", x = "k", y = "PAC (lower = more stable)") +
  theme_bw(base_size = 9)
pSil <- ggplot(st, aes(k, mean_silhouette)) + geom_line(colour = LO) + geom_point(size = 2, colour = LO) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") + scale_x_continuous(breaks = 2:6) +
  labs(title = "Mean silhouette vs k", x = "k", y = "Mean silhouette width") + theme_bw(base_size = 9)pS2 <- wrap_plots(c(hm, list(pPAC, pSil)), ncol = 4) +
  plot_annotation(title = "Consensus-clustering diagnostics (k = 2 to 6)",
                  subtitle = "Ward.D2, Euclidean, 1,000 resamples, 80% subsampling" +
                    theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 10)))
ggsave(file.path(OUT, "FigS2_consensus.pdf"), pS2, width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(OUT, "FigS2_consensus.png"), pS2, width = 10, height = 7, dpi = 300)

cat(sprintf("Fig2: PC1=%.1f%% | r(PC1,axis)=%.3f\n", pv[1], r_pc1))
cat(sprintf("Fig3: log-rank p=%.3f | events=%d\n", p_lr, sum(cl$OS_event)))
cat(sprintf("Fig5: FDR range %.3g - %.3g | max Kruskal p = %.3g\n", min(imm$FDR), max(imm$FDR), max(imm$kruskal_p)))
cat(sprintf("Fig6: DEG (FDR<0.05 & |log2FC|>1) = %d\n", sum(deg$sig != "ns")))
cat("DONE\n")
