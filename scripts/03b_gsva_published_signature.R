#!/usr/bin/env Rscript
# =============================================================================
# 03b_gsva_published_signature.R — 用同一参数对 van Groningen 2017 ADRN/MES 签名打分
# =============================================================================
# 目的：§25.4 的 `published_adrn_mes` 操作化（D4 判定量之一）。
# 参数与 03_gsva_ssgsea.R **完全一致**（ssgseaParam, alpha=0.25, normalize=TRUE），
# 保证与炎症轴分数可比。
#
# 输入：results/expr_tpm.tsv（20,584 × 153）
#       data/gene_sets/vangroningen_adrn_mes.gmt（ADRN 353 / MES 470 命中基因）
# 输出：results/gsva_published_adrn_mes_scores.tsv
#       results/gsva_published_adrn_mes_params.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(BiocManager)
})

# NOTE: run this script from the repository root.
RES <- "results"

# ---- 读入签名（GMT）----
gmt_path <- "data/gene_sets/vangroningen_adrn_mes.gmt"
lines <- readLines(gmt_path)
sets <- list()
for (ln in lines) {
  f <- strsplit(ln, "\t")[[1]]
  f <- f[f != ""]
  sets[[f[1]]] <- f[-c(1, 2)]
}
cat(sprintf("[input] 签名数: %d（%s）\n", length(sets),
            paste(sprintf("%s=%d", names(sets), sapply(sets, length)), collapse = ", ")))

m <- as.matrix(read.delim(file.path(RES, "expr_tpm.tsv"), row.names = 1, check.names = FALSE))
cat(sprintf("[input] 矩阵: %d genes x %d samples\n", nrow(m), ncol(m)))

cov <- do.call(rbind, lapply(names(sets), function(nm) {
  g <- sets[[nm]]
  hit <- intersect(g, rownames(m))
  data.frame(gene_set = nm, n_defined = length(g), n_found = length(hit),
             coverage = round(length(hit) / length(g), 4), stringsAsFactors = FALSE)
}))
print(cov)

# ---- 与主分析完全相同的参数 ----
ALPHA <- 0.25; MIN_SIZE <- 1; MAX_SIZE <- Inf; NORMALIZE <- TRUE
param <- ssgseaParam(exprData = m, geneSets = sets,
                     alpha = ALPHA, normalize = NORMALIZE,
                     minSize = MIN_SIZE, maxSize = MAX_SIZE)
t0 <- Sys.time()
sc <- gsva(param)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[gsva] 完成：%d 集 x %d 样本，耗时 %.1f 秒\n", nrow(sc), ncol(sc), el))

out <- data.frame(gene_set = rownames(sc), sc, check.names = FALSE)
write.table(out, file.path(RES, "gsva_published_adrn_mes_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

params <- data.frame(
  item = c("source_paper", "PMID", "DOI", "supplementary_table", "signature_counts",
           "R_version", "GSVA_version", "Bioconductor_version",
           "method", "alpha", "normalize", "minSize", "maxSize",
           "input_matrix", "n_genes", "n_samples", "runtime_sec", "date"),
  value = c("van Groningen T, et al. Nat Genet 2017;49(8):1261-1266",
            "28650485", "10.1038/ng.3899",
            "Supplementary Table 2 (MOESM3)",
            paste(sprintf("%s=%d", cov$gene_set, cov$n_found), collapse = "; "),
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(packageVersion("GSVA")),
            as.character(BiocManager::version()),
            "ssGSEA (Barbie et al. 2009) via GSVA::ssgseaParam",
            as.character(ALPHA), as.character(NORMALIZE),
            as.character(MIN_SIZE), as.character(MAX_SIZE),
            "results/expr_tpm.tsv",
            as.character(nrow(m)), as.character(ncol(m)),
            sprintf("%.1f", el), format(Sys.time(), "%Y-%m-%d %H:%M")),
  stringsAsFactors = FALSE)
write.table(params, file.path(RES, "gsva_published_adrn_mes_params.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[done] 03b_gsva_published_signature.R\n")
print(params)
