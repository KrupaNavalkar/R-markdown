# README: GSE147507_NGS_DESeq2_masked_signatures.Rmd

**Type:** R Markdown (knits to HTML, Word, or PDF)  
**Lines:** ~374  
**Companion render script:** `render_GSE147507_NGS_DESeq2_masked_signatures.R`

---

## Purpose

This R Markdown document is the **production-ready, reproducible report** version of the GSE147507 DESeq2 analysis. It produces a self-contained HTML report with:

1. **Raw count QC** — log2 boxplot of all 69 samples.
2. **Metadata cleaning** — structured phenotype table from the corrected CSV.
3. **DESeq2 normalization** — size-factor normalization and variance-stabilizing transformation (VST) across all samples.
4. **Masked signature scoring** — five proprietary gene expression scores are calculated in hidden chunks; visible outputs label them only as **Score 1–5** to protect proprietary gene identities. Score-to-name mapping:
   - Score 1
   - Score 2 
   - Score 3 
   - Score 4 
   - Score 5 
5. **AUC analysis** — ROC-AUC (infected vs. mock) per series and per score; outputs a CSV.
6. **DESeq2 differential expression** — five SARS-CoV-2 vs. Mock contrasts (Series 1, 2, 5, 6, 7) with volcano plots, MA plots, and heatmaps using original gene names.
7. **PCA** — across all samples using VST-transformed counts.
8. **Export** — all tables and plots saved to disk.

---

## Inputs

Supplied via YAML `params` block (defaults shown):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `raw_counts_file` | `GSE147507_RawReadCounts_Human.tsv` | Tab-delimited raw count matrix (21,797 genes × 69 samples) |
| `metadata_file` | `GSE147507_all_pheno_corrected.csv` | Corrected sample phenotype table |
| `output_dir` | `.` (current directory) | Root directory for all output files/subdirectories |
| `min_n_per_group_auc` | `3` | Minimum samples per group required to compute AUC |

---

## Expected Outputs

All paths are relative to `output_dir`:

| Output | Location | Description |
|--------|----------|-------------|
| HTML report | `GSE147507_NGS_DESeq2_masked_signatures.html` | Self-contained report (all plots embedded) |
| Raw count boxplot | `plots/GSE147507_raw_count_distribution.png` | QC plot |
| Per-series score plots | `plots/GSE147507_Series*_masked_scores.png` | One PNG per series (10 total) |
| AUC table | `GSE147507_masked_score_auc_by_series.csv` | AUC per score × series |
| DESeq2 results | `deg_outputs/GSE147507_DESeq2_results_original_gene_names.csv` | Full results for all 5 contrasts |
| Volcano plots | `deg_outputs/GSE147507_volcano_*_original_gene_names.png` | One per contrast |
| MA plots | `deg_outputs/GSE147507_MA_*_original_gene_names.png` | One per contrast |
| PCA plot | `deg_outputs/GSE147507_PCA_vst_original_gene_names.png` | VST-based PCA |
| Heatmaps | `deg_outputs/GSE147507_heatmap_top_DESeq2_*_original_gene_names.png` | Top 50 DEGs per contrast |
| Cleaned metadata | `GSE147507_cleaned_metadata.csv` | Parsed sample phenotype table |
| Normalized counts | `GSE147507_DESeq2_normalized_counts_original_gene_names.csv` | Size-factor normalized counts |
| Score + expression table | `GSE147507_expression_with_masked_scores.csv` | All samples with raw expression + 5 scores |
| Long score table | `GSE147507_masked_scores_long.csv` | Tidy format: one row per sample × score |

---

## Dependencies

**Required CRAN packages:**
```
dplyr, forcats, ggplot2, knitr, pROC, purrr, readr, stringr, tibble, tidyr, rmarkdown
```

**Required Bioconductor packages:**
```
DESeq2
```

**Optional (enhance output):**
```
pheatmap   — row-clustered heatmaps (skipped gracefully if absent)
ggrepel    — non-overlapping gene labels on volcano plots (skipped gracefully if absent)
```

Install all at once using the render script, or manually:
```r
install.packages(c("dplyr","forcats","ggplot2","knitr","pROC","purrr",
                   "readr","stringr","tibble","tidyr","rmarkdown",
                   "pheatmap","ggrepel"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

---

## How to Run

### Option A — Use the render script (recommended)

This script installs missing packages automatically, then renders the report:

```r
source("render_GSE147507_NGS_DESeq2_masked_signatures.R")
```

Run from the directory containing both the `.Rmd` and the input data files, or set your working directory first.

### Option B — Render directly from R

```r
setwd("/path/to/GSE147507/")
rmarkdown::render(
  input       = "GSE147507_NGS_DESeq2_masked_signatures.Rmd",
  output_file = "GSE147507_NGS_DESeq2_masked_signatures.html",
  envir       = new.env(parent = globalenv())
)
```

### Option C — Override parameters

To point to different input files or a different output directory:
```r
rmarkdown::render(
  input  = "GSE147507_NGS_DESeq2_masked_signatures.Rmd",
  params = list(
    raw_counts_file = "path/to/counts.tsv",
    metadata_file   = "path/to/metadata.csv",
    output_dir      = "path/to/output"
  )
)
```

### Option D — Knit in RStudio

Open the `.Rmd` in RStudio and click **Knit → Knit to HTML**.

---

## Validation / Testing

### 1. Dependency check (runs automatically)
The `setup` chunk checks for all required packages at render time and stops with an informative error if any are missing:
```
Error: Missing CRAN package(s): pROC, stringr
```
Fix by installing the listed packages, then re-render.

### 2. Input file validation
The `import-counts` and `import-metadata` chunks call `stop()` with descriptive messages if input files are missing or required columns are absent:
```
Error: Missing raw counts file: GSE147507_RawReadCounts_Human.tsv
Error: Missing metadata columns: geo_accession, characteristics_ch1
```

### 3. Check sample counts
After the metadata import chunk, a summary table is printed:
```
series     cell_line   treatment_status   n
Series1    NHBE        Infected           3
Series1    NHBE        Mock               3
...
```
Confirm: 69 samples total across 10 series, with ≥2 samples per group for the 5 contrasts used in DESeq2.

### 4. Verify DESeq2 output
After the `deseq2-deg-original-gene-names` chunk, a significance summary is printed:
```
comparison_id       significance_class        n
Series1_NHBE        Higher in comparison    XXX
Series1_NHBE        Higher in reference     XXX
Series1_NHBE        Not significant        XXXX
```
For SARS-CoV-2 vs. Mock in responsive lines (NHBE, Calu-3), expect hundreds to thousands of significant DEGs (padj < 0.05, |log2FC| > 1).

### 5. Verify AUC table
Open `GSE147507_masked_score_auc_by_series.csv`. For strongly responsive comparisons:
- Score 5 (VIRUS_ratio1) should show AUC near 1.0 for Series 1 (NHBE) and Series 7 (Calu-3).
- AUC = `NA` is expected for series with < 3 samples per group.

### 6. Confirm output files exist
After a successful render, check that the output files listed in the table above were created in the expected subdirectories (`plots/`, `deg_outputs/`).

### 7. Cross-validate with the exploratory script
Compare AUC values printed in the title of plots from `GSE147507_DeSeq2.R` with the AUC CSV from this report. Values should agree closely (note: the `.R` script uses raw counts while this report uses DESeq2 size-factor normalized counts, so minor differences are expected).
