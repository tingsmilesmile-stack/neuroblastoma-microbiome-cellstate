# Host-transcriptome microbiome-interaction scores in neuroblastoma are largely explained by tumour cell state

Analysis code for:

> **Host-Transcriptome Microbiome-Interaction Scores in Neuroblastoma Are Largely Explained by
> Tumour Cell State: A TARGET Cohort Analysis**

All analyses use the publicly available **GDC TARGET-NBL** cohort. No primary data are redistributed here.

## Requirements

- **R** >= 4.6.1: GSVA (>= 2.6.6), DESeq2, edgeR, ConsensusClusterPlus, survival, survminer,
  clusterProfiler, org.Hs.eg.db, ggplot2, pheatmap, patchwork, ggrepel, RColorBrewer
- **Python** >= 3.9: pandas, numpy, scipy, statsmodels

## Data access

Download from the GDC TARGET-NBL project: <https://portal.gdc.cancer.gov/projects/TARGET-NBL>.
The exact file and case identifiers used are listed in Additional file 7 of the manuscript.

## How to run

Run the scripts **from the repository root** (all internal paths are relative), in this order:

| # | Script | Purpose |
|---|---|---|
| 1 | `D1_merge_supp_clinical.py` | build the clinical/sample table (n = 153) |
| 2 | `03_gsva_ssgsea.R` | ssGSEA scoring of the eight pathway gene sets (GSVA) |
| 3 | `A36_fetch_published_signature.py` | retrieve the published ADRN/MES signature |
| 4 | `03b_gsva_published_signature.R` | score the published ADRN/MES signature |
| 5 | `A37_published_signature_stats.py` | published-signature association statistics |
| 6 | `06_cluster153.R` | consensus-clustering diagnostics (PAC, silhouette, bootstrap ARI) |
| 7 | `07_main153.py` | inflammatory axis, immune associations, survival analysis |
| 8 | `08_deg153.R` | DESeq2 / edgeR differential expression and GO-BP enrichment |
| 9 | `09_r7_circularity.R` | gene-overlap / circularity analysis |
| 10 | `10_r7_purity.py` | tumour-purity adjustment |
| 11 | `11_continuity.py` | continuity analysis (PCA, pairwise correlations) |
| 12 | `12_alpha_sensitivity.R` | sensitivity to the ssGSEA weighting parameter alpha |
| 13 | `13_mcp_counter.R` | orthogonal MCP-counter deconvolution |
| 14 | `15_pathway_cox153.py` | per-pathway Cox regression |
| 15 | `A38_signature_circularity.py` | circularity of the published-signature correlation |
| 16 | `D4_branch_adjudication.py` | adjudication of the cell-state confounding |
| 17 | `build_supplementary_v2.py` | assemble the supplementary tables |
| 18 | `16_figures153.R` | generate all figures |

## Reproducibility notes

- Score definitions and decision rules were fixed **before** the association analyses were performed.
- Random seeds: `42` (clustering) and `20260930` (bootstrap), stated in the scripts.
- Result tables cited in the manuscript are provided as additional files.

## License

Released under the **MIT License** — see `LICENSE`. Copyright (c) 2026 Ting Li and Xiaojuan Wu.
