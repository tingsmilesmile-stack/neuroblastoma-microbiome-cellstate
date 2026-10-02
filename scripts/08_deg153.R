#!/usr/bin/env Rscript
# 08_deg153.R — DESeq2 / edgeR differential expression for MI_high vs MI_low (n=153)
# Revision, reviewer comment #11 (IN VIVO, MANUSCRIPT NO. 8500-L):
#   replace Mann-Whitney with a count-based model (DESeq2 primary, edgeR robustness)
#   and document the origin of the "19,938 -> 15,800" gene filtering.
#
# Run:  Rscript scripts/08_deg153.R
# Reads (read-only): data/expr_all/<file_id>/*augmented_star_gene_counts.tsv
#                    data/primary_expr_ids.json
#                    results/main153_cluster_labels.tsv
# Writes: results/main153_counts_unstranded.tsv
#         results/main153_deg_deseq2.tsv
#         results/main153_deg_edger_top.tsv
#         results/main153_deg153_log.txt        (key numbers, machine-greppable)
#         results/main153_enrich_go_bp_*.tsv    (only if packages available)
#         results/main153_enrich_status.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(DESeq2)
  library(edgeR)
})

BASE <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
# fallback if --file not available (e.g. sourced)
if (is.na(BASE) || !dir.exists(file.path(BASE, "data"))) BASE <- getwd()
DAT <- file.path(BASE, "data")
RES <- file.path(BASE, "results")
dir.create(RES, showWarnings = FALSE)

set.seed(42)
FORMAT  <- "%s"
say <- function(...) cat(sprintf(...), "\n", sep = "")
log_lines <- character(0)
rec <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); log_lines <<- c(log_lines, s) }

say("=== 08_deg153.R ===  R %s", R.version.string)
say("BASE = %s", BASE)

# ---------------------------------------------------------------- 1. labels ----
lab <- fread(file.path(RES, "main153_cluster_labels.tsv"), sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE)
stopifnot(all(c("sample", "k2") %in% names(lab)))
lab <- lab[, c("sample", "k2")]
colnames(lab) <- c("sample", "condition")
lab$condition <- factor(lab$condition, levels = c("MI_low", "MI_high"))  # MI_low = reference
nhigh <- sum(lab$condition == "MI_high")
nlow  <- sum(lab$condition == "MI_low")
rec("[labels] samples=%d  MI_high=%d  MI_low=%d  (reference=MI_low)", nrow(lab), nhigh, nlow)
if (nrow(lab) != 153 || nhigh != 76 || nlow != 77) {
  stop(sprintf("Group sizes do not match the label file contract (expected 153 / 76 / 77, got %d / %d / %d).",
               nrow(lab), nhigh, nlow))
}

# ------------------------------------------------------- 2. file -> sample map ----
man <- fromJSON(file.path(DAT, "primary_expr_ids.json"), simplifyVector = FALSE)
fid2case <- setNames(vapply(man, `[[`, "", 3), vapply(man, `[[`, "", 1))
rec("[manifest] primary tumor expr ids = %d (unique file_id=%d, unique case=%d)",
    length(man), length(unique(names(fid2case))), length(unique(unname(fid2case))))

find_counts <- function(fid) {
  d <- file.path(DAT, "expr_all", fid)
  f <- list.files(d, pattern = "augmented_star_gene_counts\\.tsv$", full.names = TRUE)
  f[1]
}

# ----------------------------------------------- 3. build counts matrix (long) ----
t0 <- Sys.time()
long_list <- vector("list", length(lab$sample))
for (i in seq_along(lab$sample)) {
  case <- lab$sample[i]
  fid  <- names(fid2case)[match(case, fid2case)]
  if (is.na(fid)) stop(sprintf("no file_id for case %s", case))
  f <- find_counts(fid)
  if (is.na(f)) stop(sprintf("no counts file for file_id %s (case %s)", fid, case))
  dt <- fread(f, sep = "\t", skip = "gene_id", header = TRUE, showProgress = FALSE)
  need <- c("gene_id", "gene_name", "unstranded")
  if (!all(need %in% names(dt)))
    stop(sprintf("unexpected columns in %s: %s", f, paste(names(dt), collapse = ",")))
  dt <- dt[grepl("^ENSG", gene_id), .(gene_id, gene_name, unstranded)]
  dt[, gene_id := sub("\\.\\d+", "", gene_id)]          # drop only the .NN version, keeps _PAR_Y
  dt[, sample := case]
  long_list[[i]] <- dt
  if (i %% 25 == 0 || i == length(lab$sample)) say("  loaded %d/%d : %s", i, length(lab$sample), case)
}
long <- rbindlist(long_list, use.names = TRUE)
rm(long_list); gc(verbose = FALSE)
rec("[counts] long rows=%d  distinct gene_id=%d  (%.1fs)",
    nrow(long), uniqueN(long$gene_id), as.numeric(difftime(Sys.time(), t0, units = "secs")))

wide <- dcast(long, gene_id ~ sample, value.var = "unstranded", fun.aggregate = sum)
gene_name <- long[, .(gene_name = gene_name[1]), by = gene_id]
setkey(wide, gene_id); setkey(gene_name, gene_id)
ann <- gene_name[wide$gene_id]

gene_ids <- wide$gene_id
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- gene_ids
mat <- mat[, lab$sample, drop = FALSE]           # enforce label order
mat[is.na(mat)] <- 0L
storage.mode(mat) <- "integer"

rec("[counts] matrix = %d genes x %d samples", nrow(mat), ncol(mat))
rec("[counts] total counts over matrix = %s", format(sum(mat), big.mark = ","))
rec("[counts] zero-count genes = %d", sum(rowSums(mat) == 0))

fwrite(data.table(gene_id = gene_ids, gene_name = ann$gene_name, as.data.table(mat)),
       file.path(RES, "main153_counts_unstranded.tsv"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------- 4. DESeq2 ----
say("--- DESeq2 ---")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData   = data.frame(condition = lab$condition,
                                                     row.names = lab$sample),
                              design    = ~ condition)
n_pre <- nrow(dds)
keep <- rowSums(counts(dds)) >= 10              # DESeq2's documented default pre-filter
dds <- dds[keep, ]
n_post <- nrow(dds)
rec("[DESeq2] filter rowSums(counts)>=10 : before=%d  after=%d  removed=%d",
    n_pre, n_post, n_pre - n_post)

dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "MI_high", "MI_low"))
res <- res[order(res$padj), ]

shrunk_note <- "lfcShrink(apeglm) skipped: apeglm not installed"
if (requireNamespace("apeglm", quietly = TRUE)) {
  res <- lfcShrink(dds, coef = "condition_MI_high_vs_MI_low", type = "apeglm", quiet = TRUE)
  res <- res[order(res$padj), ]
  shrunk_note <- "lfcShrink(apeglm) applied"
} else if (requireNamespace("ashr", quietly = TRUE)) {
  res <- lfcShrink(dds, contrast = c("condition", "MI_high", "MI_low"), type = "ashr", quiet = TRUE)
  res <- res[order(res$padj), ]
  shrunk_note <- "lfcShrink(ashr) applied (apeglm unavailable)"
}
rec("[DESeq2] %s", shrunk_note)

out <- data.frame(gene_id   = rownames(res),
                  gene_name = ann$gene_name[match(rownames(res), gene_ids)],
                  baseMean  = res$baseMean,
                  log2FoldChange = res$log2FoldChange,
                  lfcSE     = res$lfcSE,
                  stat      = res$stat,
                  pvalue    = res$pvalue,
                  padj      = res$padj,
                  row.names = NULL, check.names = FALSE)
out <- out[order(out$padj), ]
fwrite(out, file.path(RES, "main153_deg_deseq2.tsv"), sep = "\t", na = "NA", quote = FALSE)

fdr05  <- sum(out$padj < 0.05, na.rm = TRUE)
sig    <- out[!is.na(out$padj) & out$padj < 0.05 & abs(out$log2FoldChange) > 1, ]
n_up   <- sum(sig$log2FoldChange > 0)
n_down <- sum(sig$log2FoldChange < 0)
rec("[DESeq2] genes tested (padj not NA)=%d  (NA from independent filtering=%d)",
    sum(!is.na(out$padj)), sum(is.na(out$padj)))
rec("[DESeq2] padj<0.05 : %d", fdr05)
rec("[DESeq2] padj<0.05 & |log2FC|>1 : %d  (up=%d  down=%d)", nrow(sig), n_up, n_down)
rec("[DESeq2] top up   : %s", paste(head(sig$gene_name[sig$log2FoldChange > 0][order(-sig$log2FoldChange[sig$log2FoldChange > 0])], 15), collapse = ", "))
rec("[DESeq2] top down : %s", paste(head(sig$gene_name[sig$log2FoldChange < 0][order(sig$log2FoldChange[sig$log2FoldChange < 0])], 15), collapse = ", "))

sf <- sizeFactors(dds)
rec("[DESeq2] sizeFactors: min=%.3f median=%.3f max=%.3f", min(sf), median(sf), max(sf))
rec("[DESeq2] dispersion fit: %s", ifelse(is.null(dds@dispersionFunction), "parametric fallback", "parametric"))

# outlier flag: samples with Cook's distance outliers
cooks <- assays(dds)[["cooks"]]
if (!is.null(cooks)) {
  n_out <- sum(apply(cooks, 2, function(z) any(z > qf(0.99, df1 = 1, df2 = ncol(dds) - 2), na.rm = TRUE)))
  rec("[DESeq2] samples with >=1 Cook's-distance outlier gene (refit applied): %d", n_out)
}

# ------------------------------------------------------------- 5. edgeR ----
say("--- edgeR (robustness) ---")
y <- DGEList(counts = mat, group = lab$condition)
keep_e <- filterByExpr(y, group = lab$condition)
y <- y[keep_e, , keep.lib.sizes = FALSE]
rec("[edgeR] filterByExpr : before=%d  after=%d  removed=%d", nrow(mat), sum(keep_e), sum(!keep_e))
y <- calcNormFactors(y)
y <- estimateDisp(y, design = model.matrix(~ lab$condition))
fit <- glmQLFit(y, design = model.matrix(~ lab$condition))
qlf <- glmQLFTest(fit, coef = 2)
tt  <- topTags(qlf, n = nrow(y))$table
tt <- tt[order(tt$FDR), ]
tt$gene_id <- rownames(tt)
tt$gene_name <- ann$gene_name[match(tt$gene_id, gene_ids)]
tt <- tt[, c("gene_id", "gene_name", "logFC", "logCPM", "F", "PValue", "FDR")]
fwrite(head(tt, 200), file.path(RES, "main153_deg_edger_top.tsv"), sep = "\t", na = "NA", quote = FALSE)

efdr  <- sum(tt$FDR < 0.05, na.rm = TRUE)
esig  <- tt[!is.na(tt$FDR) & tt$FDR < 0.05 & abs(tt$logFC) > 1, ]
e_up  <- sum(esig$logFC > 0); e_dn <- sum(esig$logFC < 0)
rec("[edgeR] FDR<0.05 : %d", efdr)
rec("[edgeR] FDR<0.05 & |log2FC|>1 : %d  (up=%d  down=%d)", nrow(esig), e_up, e_dn)
rec("[edgeR] norm factors: min=%.3f median=%.3f max=%.3f",
    min(y$samples$norm.factors), median(y$samples$norm.factors), max(y$samples$norm.factors))

# agreement with DESeq2
agree <- intersect(out$gene_id[!is.na(out$padj) & out$padj < 0.05],
                   tt$gene_id[!is.na(tt$FDR)  & tt$FDR  < 0.05])
rec("[overlap] genes significant in BOTH DESeq2(padj<0.05) and edgeR(FDR<0.05): %d", length(agree))

# ------------------------------------------------------------- 6. enrichment ----
say("--- enrichment ---")
enrich_status <- data.frame(database = c("GO_BP", "Hallmark"), status = "not_run", note = "")
has_cp <- requireNamespace("clusterProfiler", quietly = TRUE)
has_org <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
has_msig <- requireNamespace("msigdbr", quietly = TRUE)
enrich_done <- FALSE

if (has_cp && has_org && nrow(sig) > 0) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  run_go <- function(genes, tag) {
    genes <- unique(genes[!is.na(genes) & genes != ""])
    if (length(genes) < 5) { rec("[GO %s] too few genes (%d), skipped", tag, length(genes)); return(NULL) }
    eg <- suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = org.Hs.eg.db, drop = TRUE))
    rec("[GO %s] input symbols=%d  mapped Entrez=%d", tag, length(genes), nrow(eg))
    if (nrow(eg) < 5) return(NULL)
    ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, readable = TRUE)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) { rec("[GO %s] no enriched terms at BH<0.05", tag); return(NULL) }
    df <- as.data.frame(ego)
    fwrite(df, file.path(RES, sprintf("main153_enrich_go_bp_%s.tsv", tag)), sep = "\t", quote = FALSE)
    rec("[GO %s] terms at BH<0.05 = %d ; top: %s", tag, nrow(df),
        paste(head(df$Description, 5), collapse = " | "))
    df
  }
  up_df   <- run_go(sig$gene_name[sig$log2FoldChange > 0], "up")
  down_df <- run_go(sig$gene_name[sig$log2FoldChange < 0], "down")
  both_df <- run_go(sig$gene_name, "all")
  enrich_done <- !is.null(up_df) || !is.null(down_df) || !is.null(both_df)
  enrich_status$status[enrich_status$database == "GO_BP"] <- ifelse(enrich_done, "done", "ran_no_terms")
  enrich_status$note[enrich_status$database == "GO_BP"] <-
    sprintf("clusterProfiler=%s org.Hs.eg.db=%s; genes: up=%d down=%d",
            packageVersion("clusterProfiler"), packageVersion("org.Hs.eg.db"), n_up, n_down)
} else {
  enrich_status$note[enrich_status$database == "GO_BP"] <-
    sprintf("clusterProfiler=%s org.Hs.eg.db=%s", has_cp, has_org)
}

if (has_msig) {
  enrich_status$status[enrich_status$database == "Hallmark"] <- "available_not_run"
  enrich_status$note[enrich_status$database == "Hallmark"] <- "msigdbr present; not run in this pass"
} else {
  enrich_status$status[enrich_status$database == "Hallmark"] <- "not_run"
  enrich_status$note[enrich_status$database == "Hallmark"] <- "msigdbr NOT installed"
}
fwrite(enrich_status, file.path(RES, "main153_enrich_status.tsv"), sep = "\t")
rec("[enrich] GO_BP done=%s ; Hallmark=%s", enrich_done,
    ifelse(has_msig, "available", "NOT RUN (msigdbr missing)"))

# top 30 up/down by name (always recorded, used if enrichment unavailable)
top_up   <- head(out$gene_name[!is.na(out$padj) & out$padj < 0.05 & out$log2FoldChange > 1][order(-out$log2FoldChange[!is.na(out$padj) & out$padj < 0.05 & out$log2FoldChange > 1])], 30)
top_down <- head(out$gene_name[!is.na(out$padj) & out$padj < 0.05 & out$log2FoldChange < -1][order(out$log2FoldChange[!is.na(out$padj) & out$padj < 0.05 & out$log2FoldChange < -1])], 30)
rec("[top30 up]   %s", paste(top_up, collapse = ", "))
rec("[top30 down] %s", paste(top_down, collapse = ", "))

# ------------------------------------------------------------- 7. versions ----
say("--- versions ---")
vers <- c(R = R.version.string,
          DESeq2 = as.character(packageVersion("DESeq2")),
          edgeR = as.character(packageVersion("edgeR")),
          limma = as.character(packageVersion("limma")),
          clusterProfiler = ifelse(has_cp, as.character(packageVersion("clusterProfiler")), "NA"),
          org.Hs.eg.db = ifelse(has_org, as.character(packageVersion("org.Hs.eg.db")), "NA"),
          msigdbr = ifelse(has_msig, as.character(packageVersion("msigdbr")), "NOT INSTALLED"),
          apeglm = ifelse(requireNamespace("apeglm", quietly = TRUE), as.character(packageVersion("apeglm")), "NOT INSTALLED"),
          ashr = ifelse(requireNamespace("ashr", quietly = TRUE), as.character(packageVersion("ashr")), "NOT INSTALLED"))
for (k in names(vers)) rec("[version] %s = %s", k, vers[[k]])
rec("[seed] set.seed(42)")

writeLines(log_lines, file.path(RES, "main153_deg153_log.txt"))
say("=== done ===")
