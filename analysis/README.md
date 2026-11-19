# Post-counts analysis (R / Quarto)

`phallusia_expression_analysis.qmd` takes the merged Salmon gene tables (`results/star_salmon/`, written by `03d_merge_quant.sbatch` of the pipeline; nf-core/rnaseq output works too) and the app database and produces:

| Output (in `analysis_output/`) | Content |
|---|---|
| `gene_summary.csv` | every gene: LRT statistic and adjusted p-value, temporal yes/no, expression-pattern cluster, peak stage, tau |
| `cluster_go_enrichment.csv` | GO enrichment of each cluster (propagated GO terms, hypergeometric test, BH) |
| `stage_markers.csv` | up to 25 stage-specific genes per stage (tau and TPM thresholds) |
| `stage_mean_vst.csv`, `cluster_centroids.csv`, `qc_samples.csv` | stage means, cluster profiles, library size and genes detected |

## Model

DESeq2 likelihood-ratio test, full model `~ replicate + stage` against `~ replicate`. `replicate` is the time course a library
belongs to: the two replicates of PRJEB25149 are two separate time courses (batches of embryos), and every library is more similar to
the same series at the neighbouring stage than to the other library of its own stage. PCA and sample correlations use VST values with
this series effect removed.

## Run it

R packages: DESeq2, SummarizedExperiment, ggplot2, pheatmap, cluster, DBI, RSQLite; GO.db and patchwork are optional (GO.db: parent terms and term names for the enrichment).

```bash
cd analysis
# paths are at the top of the notebook; or set them for one run:
PHM_RESULTS=/path/to/results/star_salmon PHM_DB=../app/data/phallusia.sqlite quarto render phallusia_expression_analysis.qmd
```
In RStudio: open the `.qmd` and press Render. It needs `salmon.merged.gene_counts.tsv` and `salmon.merged.gene_tpm.tsv`; with `salmon.merged.gene_lengths.tsv` (or the nf-core `.rds`) DESeq2 also uses transcript-length offsets.

## Tested how

`bash test/run_test.sh` simulates merged Salmon tables with planted patterns (6 temporal patterns of 150 genes, 2,400 flat genes, 150 genes expressed in a single stage, real GO terms planted per pattern, and a per-gene shift between the two replicate series as in the real data), runs all notebook code with knitr, with each of the three input types (SummarizedExperiment; counts + TPM + gene lengths; counts + TPM only), and checks what is recovered. Results of the last run (all three inputs gave the same picture):

- 100% of planted genes with base mean >= 20 called temporal; 99% of all planted genes
- **false discovery proportion 3.1-3.3% at a nominal 1%** (only 2 replicates per stage; read adjusted p-values as a ranking); the replicate-series shift is not called temporal (Nov 2026 run, model `~ replicate + stage`)
- cluster purity 99.9% at k = 10; 12 of 12 planted GO terms enriched in the right cluster; 150 of 150 single-stage genes found as markers
- not tested: Quarto itself (not installable in the test sandbox; Quarto uses the same knitr engine), and real data
