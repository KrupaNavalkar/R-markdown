# README: GSE70311_GEO_original_genes_masked_signatures.Rmd

**Date:** 2020-05-19  
**Dataset:** GSE70311 (Illumina HumanHT-12 v4, GPL10558)

---

## Purpose

This R Markdown script reproduces a full GEO analysis workflow for **GSE70311**, a sepsis/SIRS longitudinal microarray dataset. It:

- Downloads raw expression and metadata directly from GEO
- Harmonizes class labels: **CASE** (Sepsis) vs **Control group** (SIRS)
- Annotates Illumina probes to gene symbols using `illuminaHumanv4.db`
- Computes five proprietary gene expression **signatures** (Score 1–5) with calculation code hidden in the rendered HTML output (labels and gene components are masked)
- Runs genome-wide **differential expression analysis** (DEG) via `limma` using original GEO gene symbols
- Produces publication-quality plots: volcano, MA, heatmap, PCA, and per-signature trajectory/boxplots
- Exports all results to CSV and PNG files

---

## Dependencies

### R Packages

| Source | Packages |
|---|---|
| CRAN | `dplyr`, `forcats`, `ggplot2`, `knitr`, `pROC`, `purrr`, `readr`, `stringr`, `tibble`, `tidyr`, `rmarkdown` |
| Bioconductor | `AnnotationDbi`, `Biobase`, `GEOquery`, `illuminaHumanv4.db`, `limma` |
| Optional (enhance output) | `pheatmap` (better heatmap), `ggrepel` (volcano gene labels) |

Install missing packages before running:

```r
# CRAN
install.packages(c("dplyr","forcats","ggplot2","knitr","pROC","purrr",
                   "readr","stringr","tibble","tidyr","rmarkdown",
                   "pheatmap","ggrepel"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi","Biobase","GEOquery","illuminaHumanv4.db","limma"))
```

### Environment

- R ≥ 4.0
- Internet access (for GEO download on first run)
- ~500 MB disk space for GEO data cache and outputs

---

## Inputs

| Input | Source | Notes |
|---|---|---|
| GSE70311 expression matrix | Downloaded automatically from NCBI GEO | ~47,000 Illumina probes × ~300 samples |
| Sample metadata | Downloaded automatically from NCBI GEO | Source name and characteristics columns used |
| Illumina probe annotation | `illuminaHumanv4.db` R package | Maps probe IDs → gene symbols |

No local input files are required. All data is fetched at runtime via `GEOquery::getGEO()`.

---

## How to Run

### Option 1 — RStudio

Open the `.Rmd` file in RStudio and click **Knit**, or run from the R console:

```r
rmarkdown::render(
  input      = "GSE70311_GEO_original_genes_masked_signatures.Rmd",
  output_file = "GSE70311_GEO_original_genes_masked_signatures.html",
  envir      = new.env(parent = globalenv())
)
```

### Option 2 — Command line

```bash
Rscript -e "rmarkdown::render('GSE70311_GEO_original_genes_masked_signatures.Rmd')"
```

### Parameters (optional overrides)

The script accepts YAML params — override at render time if needed:

```r
rmarkdown::render(
  "GSE70311_GEO_original_genes_masked_signatures.Rmd",
  params = list(
    gse_accession = "GSE70311",   # GEO accession
    gpl_platform  = "GPL10558",   # Platform filter
    output_dir    = "./results"   # Where to save outputs
  )
)
```

---

## Expected Outputs

All files are written to `output_dir` (default: current directory).

```
output_dir/
├── GSE70311_all_metadata_from_GEO.csv                          # Raw GEO sample metadata
├── GSE70311_analysis_matrix_original_gene_names_masked_signatures.csv  # Full analysis matrix
├── GSE70311_signature_scores_long_masked_labels.csv            # Signature scores (long format)
├── GSE70311_signature_overall_auc_masked.csv                   # AUC per signature (overall)
├── GSE70311_signature_timepoint_auc_masked.csv                 # AUC per signature × time point
├── GSE70311_signature_plot_objects_masked_labels.rds           # R plot objects (list)
├── plots/
│   ├── GSE70311_Score 1_case_control.png                       # Boxplot: CASE vs Control
│   ├── GSE70311_Score 1_timepoint.png                          # Boxplot by time point
│   ├── GSE70311_Score 1_trajectory.png                         # Patient trajectory lines
│   └── ... (3 plots × 5 signatures = 15 PNGs total)
└── deg_outputs/
    ├── GSE70311_DEG_CASE_vs_Control_group_original_gene_names.csv  # Full DEG table
    ├── GSE70311_DEG_summary_original_gene_names.png
    ├── GSE70311_volcano_CASE_vs_Control_group_original_gene_names.png
    ├── GSE70311_MA_CASE_vs_Control_group_original_gene_names.png
    ├── GSE70311_heatmap_top_DEG_CASE_vs_Control_group_original_gene_names.png
    └── GSE70311_PCA_top_DEG_CASE_vs_Control_group_original_gene_names.png
```

The rendered HTML report is self-contained with embedded figures and masked signature code.

---

## Validation & Testing

### Built-in checks (run automatically)

| Check | Location | What it verifies |
|---|---|---|
| Missing package guard | `setup` chunk | Stops with a clear error if any required package is absent |
| GPL platform match | `download-geo` chunk | Stops if `GPL10558` is not found in the GEO object |
| Log2 transform detection | `download-geo` chunk | Auto-detects and applies log2 if data is on a linear scale |
| Probe filter | `gene-symbols` chunk | Removes probes with no gene symbol and `LOC*` loci; keeps highest-expressed probe per gene |
| DEG finite/variance filter | `deg-original-gene-names` chunk | Drops genes with non-finite values or zero variance before limma |
| `safe_auc()` function | `setup` chunk | Returns `NA` (no crash) if a group is missing or sample size is too small |

### Manual sanity checks after running

1. **Sample counts** — the `metadata` chunk prints a table of `case_status × time_point`. Expect ~15 Sepsis and ~15 SIRS patients across up to 9 time points.

2. **Expression matrix dimensions** — the `gene-symbols` chunk prints `n_genes` and `n_samples`. Expect ~14,000–18,000 genes and ~280–320 samples after probe filtering.

3. **DEG summary table** — the `deg-original-gene-names` chunk prints counts of "Higher in CASE", "Higher in Control group", and "Not significant". A reasonable result has hundreds to low thousands of significant genes at adj. p < 0.05 and |logFC| > 1.

4. **AUC table** — `overall_auc` should show AUC values between 0.5 and 1.0 for all five scores. AUC < 0.5 or `NA` for a signature indicates the required genes were not found in this dataset.

5. **Output files** — verify all 6 files in `deg_outputs/` and 15 PNGs in `plots/` are created and non-empty.

6. **HTML report** — open in a browser and confirm signature calculation code is hidden (only boxplots and AUC tables visible under "Hidden signature calculations").

---

## Key Design Notes

- **Signature masking:** The five signature formulas are computed in an `echo=FALSE` chunk, so they are invisible in the rendered HTML. Scores appear only as "Score 1" through "Score 5".
- **Gene symbol resolution:** When multiple probes map to the same gene, the probe with the highest mean expression is retained.
- **Time point harmonization:** Day 6 → Day 7; Days 26/28 → Day 21 (consolidates sparse late time points).
- **Diagnosis status:** Four patients (P09, P16, P36, P42, P47) have specific time points reclassified as `CASE` vs `Pre-CASE` based on clinical progression logic hard-coded in the metadata chunk.
