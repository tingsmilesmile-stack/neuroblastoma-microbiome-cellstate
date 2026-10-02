#!/usr/bin/env Rscript
# 06_cluster153.R — 153 例主分析：标准共识聚类 + 完整稳定性证据（回应 R9）
# 输入：results/gsva_ssgsea_scores.tsv（GSVA 标准实现，8 MI 通路 + 8 免疫，153 例）
# 输出：results/main153_*（标签、PAC 曲线、共识矩阵、簇内共识、silhouette、bootstrap ARI）
suppressPackageStartupMessages({
  library(ConsensusClusterPlus); library(cluster)
})
BASE <- "."
RES  <- file.path(BASE, "results")
MI_PATHWAYS <- c("TLR_signaling","NOD_like_receptor","Antimicrobial_peptides",
                 "LPS_inflammatory","SCFA_butyrate_response","Bile_acid_metabolism",
                 "Tryptophan_IDO_kynurenine","Mucin_barrier")
set.seed(42)

raw <- read.delim(file.path(RES,"gsva_ssgsea_scores.tsv"), row.names=1, check.names=FALSE)
stopifnot(all(MI_PATHWAYS %in% rownames(raw)))
X <- t(scale(t(as.matrix(raw[MI_PATHWAYS, , drop=FALSE]))))   # 每通路跨样本 z
X <- t(X)                                                    # 样本 × 通路
cat("[info] 样本", nrow(X), "｜通路", ncol(X), "\n")

# ---- 标准共识聚类（ConsensusClusterPlus, Ward.D2, 1000×0.8 重抽样）----
res <- ConsensusClusterPlus(t(X), maxK=6, reps=1000, pItem=0.8, pFeature=1,
                            clusterAlg="hc", distance="euclidean",
                            innerLinkage="ward.D2", finalLinkage="ward.D2",
                            seed=42, plot=NULL, verbose=FALSE)

pac_of <- function(M){ v <- M[upper.tri(M)]; mean(v > 0.1 & v < 0.9) }
icl <- calcICL(res, plot=NULL)
icl_cl <- as.data.frame(icl$clusterConsensus)
names(icl_cl) <- c("k","cluster","clusterConsensus")
stab <- list(); cons_k2 <- NULL
for (k in 2:6){
  M <- res[[k]]$consensusMatrix; rownames(M) <- colnames(M) <- colnames(t(X))
  pac <- pac_of(M)
  lab <- res[[k]]$consensusClass
  # 簇内共识
  cc <- icl_cl[icl_cl$k==k, , drop=FALSE]
  # silhouette（以 (1-共识) 为距离）
  d <- as.dist(1 - M)
  sil <- silhouette(lab, d); mean_sil <- mean(sil[,3]); min_sil <- min(sil[,3])
  # bootstrap ARI（重抽样后对全体样本聚类 vs 参考标签）
  aris <- numeric(100)
  for (b in 1:100){
    idx <- sample(nrow(X), floor(0.8*nrow(X))); Xb <- X[idx, , drop=FALSE]
    hb <- hclust(dist(Xb), method="ward.D2")
    lb <- cutree(hb, k)
    a <- table(lab[idx], lb)
    # 手算 ARI
    n <- sum(a); s <- 0
    for (i in 1:nrow(a)) for (j in 1:ncol(a)) s <- s + choose(a[i,j],2)
    sa <- sum(sapply(rowSums(a), choose, 2)); sb <- sum(sapply(colSums(a), choose, 2))
    expct <- sa*sb/choose(n,2); mx <- (sa+sb)/2; ari <- (s-expct)/(mx-expct)
    aris[b] <- ari
  }
  stab[[length(stab)+1]] <- data.frame(
    k=k, PAC=round(pac,4),
    mean_cluster_consensus=round(mean(cc$clusterConsensus),4),
    min_cluster_consensus=round(min(cc$clusterConsensus),4),
    mean_silhouette=round(mean_sil,4), min_silhouette=round(min_sil,4),
    bootstrap_ARI_median=round(median(aris),4),
    bootstrap_ARI_q025=round(quantile(aris,.025),4),
    bootstrap_ARI_q975=round(quantile(aris,.975),4))
  if (k==2){ cons_k2 <- M; k2lab <- lab }
  write.table(M, file.path(RES, sprintf("main153_consensus_k%d.tsv", k)), sep="\t", quote=FALSE)
}
stab <- do.call(rbind, stab)
write.table(stab, file.path(RES,"main153_stability.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n=== 稳定性（PAC 越低越稳；ARI 越高越稳）===\n"); print(stab)

# 选 k（PAC 最小）
kbest <- stab$k[which.min(stab$PAC)]
cat("[info] PAC 最小 k =", kbest, "\n")

# ---- final labels (k = 2 and kbest; the main analysis uses k = 2) ----
mean_act <- rowMeans(X[, MI_PATHWAYS])
n2 <- as.integer(factor(res[[2]]$consensusClass)); act2 <- tapply(mean_act, n2, mean)
lab2 <- ifelse(n2 == names(act2)[which.max(act2)], "MI_high", "MI_low")
out <- data.frame(sample=rownames(X), MI_axis_z=round(mean_act,4), k2=lab2)
for (k in 3:6){
  nk <- as.integer(factor(res[[k]]$consensusClass)); actk <- tapply(mean_act, nk, mean)
  ord <- order(actk); nm <- character(length(actk)); nm[ord] <- paste0("C", seq_along(actk))
  out[[paste0("k",k)]] <- nm[nk]
}
write.table(out, file.path(RES,"main153_cluster_labels.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n=== MI-high/low 计数 ===\n"); print(table(out$k2))

# ---- per-pathway direction under k = 2 ----
dircheck <- data.frame(pathway=MI_PATHWAYS,
  mean_MIhigh=round(sapply(MI_PATHWAYS, function(p) mean(X[lab2=="MI_high",p])),3),
  mean_MIlow =round(sapply(MI_PATHWAYS, function(p) mean(X[lab2=="MI_low", p])),3))
dircheck$diff <- round(dircheck$mean_MIhigh - dircheck$mean_MIlow, 3)
write.table(dircheck, file.path(RES,"main153_pathway_direction.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n=== 各通路 MI-high − MI-low（z 差）===\n"); print(dircheck)
cat("[done] 06_cluster153\n")
