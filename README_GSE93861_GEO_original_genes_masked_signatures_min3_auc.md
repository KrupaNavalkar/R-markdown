# README: GSE93861_GEO_original_genes_masked_signatures_min3_auc.Rmd

**Author:** Krupa Navalkar  
**Date:** 2020-05-19

---

## Purpose

This R Markdown workflow performs a complete gene expression analysis of the GEO dataset **GSE93861** (Ebola virus disease, platform GPL6480). It:

- Downloads the dataset directly from GEO and applies log2 transformation if needed
- Harmonizes sample labels: Ebola-derived samples → `CASE`, healthy samples → `Control group`
- Collapses multi-probe microarray data to unique gene symbols (highest mean expression probe retained)
- Averages technical replicates and bins CASE samples into post-infection day windows (Days 7–13, 14–20, 21–33, 34–58, 59–180, 181–270)
- Computes five gene-combination signature scores (LAB, VIRUS, BACT, TRIAGE, RAPID) with signature names and component genes intentionally **masked** in rendered output
- Calculates AUC (CASE vs. Control, and per day window) using a safe AUC helper that returns `NA` when either group has fewer than 3 samples
- Runs genome-wide differential expression (limma) using original GEO gene symbols
- Produces volcano, MA, heatmap (top 50 DEGs), and PCA (top 500 DEGs) plots

---

## Dependencies

### R packages

| Source | Packages |
|--------|----------|
| CRAN | `data.table`, `dplyr`, `forcats`, `ggplot2`, `knitr`, `pROC`, `purrr`, `readr`, `stringr`, `tibble`, `tidyr`, `rmarkdown` |
| Bioconductor | `Biobase`, `GEOquery`, `limma` |
| Optional | `pheatmap` (better heatmap), `ggrepel` (gene labels on volcano plot) |

Install missing packages:

```r
# CRAN
install.packages(c("data.table", "dplyr", "forcats", "ggplot2", "knitr", "pROC",
                   "purrr", "readr", "stringr", "tibble", "tidyr", "rmarkdown",
                   "pheatmap", "ggrepel"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biobase", "GEOquery", "limma"))
```

### Environment

- R ≥ 4.0
- Internet access (GEO download on first run; results cached by `GEOquery` in `~/GEOquery_cache/` by default)

---

## Inputs

| Input | Source | Notes |
|-------|--------|-------|
| GEO dataset GSE93861 | Downloaded automatically via `GEOquery::getGEO()` | Requires internet on first run |
| `gse_accession` param | `"GSE93861"` (default) | Changeable via YAML params |
| `gpl_platform` param | `"GPL6480"` (default) | Selects the correct platform if multiple exist |
| `output_dir` param | `"."` (default) | Directory for all outputs |

No local data files are required.

---

## How to Run

**Option 1 — RStudio:** Open the `.Rmd` file and click **Knit**.

**Option 2 — R console:**

```r
rmarkdown::render(
  input      = "GSE93861_GEO_original_genes_masked_signatures_min3_auc.Rmd",
  output_file = "GSE93861_GEO_original_genes_masked_signatures.html",
  envir      = new.env(parent = globalenv())
)
```

**Option 3 — Custom output directory:**

```r
rmarkdown::render(
  input      = "GSE93861_GEO_original_genes_masked_signatures_min3_auc.Rmd",
  output_file = "GSE93861_GEO_original_genes_masked_signatures.html",
  params     = list(output_dir = "/path/to/results"),
  envir      = new.env(parent = globalenv())
)
```

---

## Expected Outputs

All files are written to `output_dir` (default: working directory).

| File | Description |
|------|-------------|
| `GSE93861_GEO_original_genes_masked_signatures.html` | Self-contained HTML report |
| `GSE93861_all_metadata_from_GEO.csv` | Raw sample metadata from GEO |
| `GSE93861_analysis_matrix_original_gene_names_masked_signatures.csv` | Analysis-ready expression + signature matrix |
| `GSE93861_signature_scores_long_masked_labels.csv` | Signature scores in long format (masked labels) |
| `GSE93861_signature_case_control_auc_masked.csv` | AUC table: CASE vs. Control group |
| `GSE93861_signature_day_window_auc_masked.csv` | AUC table: each day window vs. Control group |
| `GSE93861_signature_plot_objects_masked_labels.rds` | R list of all signature ggplot objects |
| `deg_outputs/GSE93861_DEG_CASE_vs_Control_group_original_gene_names.csv` | Full limma DEG results |
| `plots/Score*_case_control.png` | Boxplots per signature (CASE vs. Control) |
| `plots/Score*_day_windows.png` | Boxplots per signature by day window |
| `plots/Score*_time_course.png` | Time-course line plots per signature |
| `deg_outputs/GSE93861_DEG_summary_original_gene_names.png` | Bar chart of DEG counts by class |
| `deg_outputs/GSE93861_volcano_*.png` | Volcano plot |
| `deg_outputs/GSE93861_MA_*.png` | MA plot |
| `deg_outputs/GSE93861_heatmap_*.png` | Heatmap of top 50 DEGs |
| `deg_outputs/GSE93861_PCA_*.png` | PCA of top 500 DEGs |

---

## Validation and Checks

### Built-in checks (run automatically)
- The setup chunk halts with an informative `stop()` if any required package is missing
- A `stop()` is raised if the `GENE_SYMBOL` column is absent from feature metadata
- A `stop()` is raised if the limma design matrix does not contain both `CASE` and `Control.group` columns
- The `safe_auc()` helper returns `NA` (rather than erroring) when either group has fewer than 3 samples — AUC values of `NA` in the output tables indicate this condition

### Manual checks after running
1. **Sample counts** — the rendered report prints a `count(case_status)` table; expect both `CASE` and `Control group` rows
2. **Probe collapse** — a summary table shows `n_genes_after_probe_collapse`; a reasonable value for GPL6480 is ~15,000–20,000 genes
3. **Log2 flag** — the download summary reports `log2_transform_applied`; verify this matches the known data scale for GSE93861 (should be `TRUE` for raw intensity data)
4. **Signature NAs** — if any signature column is entirely `NA`, the relevant gene was not found after probe collapse; check gene names in the expression matrix
5. **AUC range** — all non-NA AUC values in the output CSVs should be between 0 and 1
6. **DEG table** — open `GSE93861_DEG_CASE_vs_Control_group_original_gene_names.csv` and confirm `adj.P.Val` values are in [0, 1] and `logFC` values span a reasonable range (approximately −5 to +5 for microarray data)
7. **Plots directory** — confirm `plots/` contains 15 PNG files (3 per signature × 5 signatures) and `deg_outputs/` contains 4 PNG files plus the DEG CSV
