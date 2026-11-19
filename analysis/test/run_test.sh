#!/bin/bash
# Self-test of the analysis notebook on SIMULATED data with planted patterns.
#   bash run_test.sh [work_dir]
# Needs R with DESeq2, SummarizedExperiment, GO.db, RSQLite, ggplot2, pheatmap, cluster, knitr.
# It runs the notebook code with knitr (the engine Quarto uses) with the three input types it understands
# (SummarizedExperiment; count + TPM + gene-length tables; count + TPM tables), and checks what was recovered
# (check_recovery.R).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ANALYSIS="$(dirname "${HERE}")"
WORK="${1:-${HERE}/_work}"
META="$(dirname "${ANALYSIS}")/sample_metadata.csv"
rm -rf "${WORK}"; mkdir -p "${WORK}"

Rscript "${HERE}/make_test_data.R" "${WORK}/data" "${META}"

run_nb() {   # run_nb <name> <expected text in the "Input:" line> : knits a copy of the notebook in ${WORK}/<name>
    local name="$1" expect="$2"
    mkdir -p "${WORK}/${name}"; cp "${ANALYSIS}/phallusia_expression_analysis.qmd" "${WORK}/${name}/"
    ( cd "${WORK}/${name}"
      PHM_RESULTS="${WORK}/data/results/star_salmon" PHM_METADATA="${META}" PHM_DB="${WORK}/data/app.sqlite" \
      PHM_OUT="${WORK}/${name}/analysis_output" \
      Rscript -e 'knitr::opts_chunk$set(dev = "png", dpi = 100); knitr::knit("phallusia_expression_analysis.qmd", output = "nb.md", quiet = TRUE)' \
        > knit.log 2>&1 ) || { echo "notebook FAILED (${name}); see ${WORK}/${name}/knit.log"; tail -20 "${WORK}/${name}/knit.log"; exit 1; }
    grep -q "Input: ${expect}" "${WORK}/${name}/nb.md" || { echo "notebook used another input than expected (${expect})"; grep "Input:" "${WORK}/${name}/nb.md"; exit 1; }
    echo "== notebook ran (${name}); checks:"
    Rscript "${HERE}/check_recovery.R" "${WORK}/data" "${WORK}/${name}/analysis_output"
}

RS="${WORK}/data/results/star_salmon"

echo "##### input = SummarizedExperiment (length offsets)"
run_nb with_se "SummarizedExperiment"

echo; echo "##### input = count + TPM + gene-length tables (what 03d_merge_quant.sbatch writes)"
mv "${RS}/salmon.merged.gene.SummarizedExperiment.rds" "${WORK}/se.rds.bak"
run_nb tables_with_lengths "count, TPM and gene-length tables"

echo; echo "##### input = count + TPM tables only"
mv "${RS}/salmon.merged.gene_lengths.tsv" "${WORK}/len.tsv.bak"
run_nb tables_only "count table, rounded"
mv "${WORK}/se.rds.bak" "${RS}/salmon.merged.gene.SummarizedExperiment.rds"
mv "${WORK}/len.tsv.bak" "${RS}/salmon.merged.gene_lengths.tsv"
echo; echo "Test finished. Figures and tables are in ${WORK}/with_se"
