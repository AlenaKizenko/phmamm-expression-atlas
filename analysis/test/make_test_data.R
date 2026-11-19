# Simulated STAR + Salmon (merged gene tables) output with PLANTED temporal patterns and GO terms, to test the analysis notebook.
# All values are FAKE. Usage:  Rscript make_test_data.R <output_dir> [path/to/sample_metadata.csv]
#
# Creates in <output_dir>:
#   results/star_salmon/salmon.merged.gene_counts.tsv, salmon.merged.gene_tpm.tsv, salmon.merged.gene_lengths.tsv
#     (the tables written by 03d_merge_quant.sbatch / merge_salmon.py; same layout as nf-core/rnaseq)
#   results/star_salmon/salmon.merged.gene.SummarizedExperiment.rds   (nf-core only; assays counts, tpm, length, ...)
#   app.sqlite   (tables genes + eggnog, the only ones the notebook reads)
#   truth.csv    (gene_id, planted cluster 0 = flat, 1-6 = planted pattern, base mean)
#   truth_go.csv (planted cluster, GO id)
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args) >= 1) args[1] else "test_data"
meta_file <- if (length(args) >= 2) args[2] else "../../sample_metadata.csv"
suppressPackageStartupMessages({ library(SummarizedExperiment); library(GO.db) })
set.seed(42)

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[order(meta$dev_order, meta$replicate), ]
S <- nrow(meta)
stage_i <- as.integer(factor(meta$dev_order))          # 1..15
n_stage <- max(stage_i)
s <- seq_len(n_stage)

# six planted patterns (log2 scale), 150 genes each; 2400 flat genes; 100 very low genes
patterns <- list(
  maternal  = 3 * exp(-(s - 1) / 3),
  cleavage  = 3 * exp(-((s - 6.5)^2) / (2 * 1.5^2)),
  neurula   = 3 * exp(-((s - 9.5)^2) / (2 * 1.3^2)),
  tailbud   = 3 * exp(-((s - 12.5)^2) / (2 * 1.3^2)),
  late_rise = 3 * plogis((s - 11) / 1.2),
  larva     = 3 * exp(-((s - 15)^2) / (2 * 1.2^2))
)
n_per <- 150; n_flat <- 2400; n_low <- 100; n_spec <- 10      # n_spec: stage-specific genes per stage
G <- length(patterns) * n_per + n_flat + n_low + n_spec * n_stage
gene_id <- sprintf("LOC%06d", seq_len(G))
# plant: 1-6 = pattern, 0 = flat, -1 = very low, 101-115 = expressed in one stage only (100 + stage number)
plant <- c(rep(seq_along(patterns), each = n_per), rep(0L, n_flat), rep(-1L, n_low), 100L + rep(s, each = n_spec))

base <- exp(rnorm(G, 4.3, 1.5)); base[plant == -1] <- runif(n_low, 0.05, 0.6)
base[plant > 100]  <- runif(sum(plant > 100), 0.05, 0.3)      # near-zero background
disp <- 0.05 + 0.4 / sqrt(pmax(base, 1))
sf   <- runif(S, 0.7, 1.3)
lfc  <- matrix(0, G, S)
for (p in seq_along(patterns)) lfc[plant == p, ] <- rep(patterns[[p]][stage_i], each = sum(plant == p))
for (g in which(plant > 100)) lfc[g, ] <- 9 * (stage_i == plant[g] - 100)       # ~500-fold in one stage
# replicate-series (batch) effect, as in the real data: each replicate is its own time course.
# Every gene gets its own shift between the two series (sd 0.5 log2); it must NOT be called temporal.
batch_lfc <- rnorm(G, 0, 0.5)
lfc[, meta$replicate == 2] <- lfc[, meta$replicate == 2] + batch_lfc
mu <- base * outer(rep(1, G), sf) * 2^lfc
counts <- matrix(rnbinom(G * S, mu = mu, size = 1 / disp), G, S,
                 dimnames = list(gene_id, meta$sample))

glen <- exp(rnorm(G, 7.6, 0.5))
len  <- glen * matrix(runif(G * S, 0.95, 1.05), G, S, dimnames = list(gene_id, meta$sample))
rpk  <- counts / len
tpm  <- t(t(rpk) / colSums(rpk)) * 1e6

dir.create(file.path(out, "results", "star_salmon"), recursive = TRUE, showWarnings = FALSE)
rs <- file.path(out, "results", "star_salmon")
frac <- counts + matrix(runif(G * S, -0.4, 0.4), G, S); frac[frac < 0] <- 0     # Salmon estimates are fractional
write_tab <- function(m, f) {
  d <- data.frame(gene_id = rownames(m), gene_name = rownames(m), round(m, 3), check.names = FALSE)
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_tab(frac, file.path(rs, "salmon.merged.gene_counts.tsv"))
write_tab(tpm,  file.path(rs, "salmon.merged.gene_tpm.tsv"))
write_tab(len,  file.path(rs, "salmon.merged.gene_lengths.tsv"))
se <- SummarizedExperiment(assays = list(counts = counts + 0, tpm = tpm, length = len,
                                         counts_scaled = counts + 0, counts_length_scaled = counts + 0),
                           colData = DataFrame(sample = meta$sample, row.names = meta$sample))
saveRDS(se, file.path(rs, "salmon.merged.gene.SummarizedExperiment.rds"))

# --- GO terms: real GO ids. Planted terms for each pattern (65% of the pattern's genes) + random background terms
planted_go <- list(
  c("GO:0006412", "GO:0042254"),   # translation, ribosome biogenesis
  c("GO:0006260", "GO:0000278"),   # DNA replication, mitotic cell cycle
  c("GO:0007399", "GO:0007411"),   # nervous system development, axon guidance
  c("GO:0001525", "GO:0007517"),   # angiogenesis, muscle organ development
  c("GO:0007155", "GO:0030198"),   # cell adhesion, extracellular matrix organization
  c("GO:0007600", "GO:0006936")    # sensory perception, muscle contraction
)
stopifnot(all(unlist(planted_go) %in% keys(GO.db)))
bp <- keys(GO.db, keytype = "GOID")
bp <- bp[Ontology(bp) %in% "BP"]
pool <- setdiff(sample(bp, 400), unlist(planted_go))
go_of <- vector("list", G); names(go_of) <- gene_id
for (g in seq_len(G)) {
  terms <- character()
  if (runif(1) < 0.6) terms <- sample(pool, sample(1:4, 1))
  if (plant[g] %in% 1:6) for (t in planted_go[[plant[g]]]) if (runif(1) < 0.65) terms <- c(terms, t)
  go_of[[g]] <- unique(terms)
}
has_go <- lengths(go_of) > 0

library(DBI); library(RSQLite)
db <- file.path(out, "app.sqlite"); unlink(db)
con <- dbConnect(SQLite(), db)
dbWriteTable(con, "genes", data.frame(gene_id = gene_id, symbol = gene_id,
             description = paste("simulated gene", plant), biotype = "protein_coding"))
dbWriteTable(con, "eggnog", data.frame(gene_id = gene_id[has_go],
             gos = vapply(go_of[has_go], paste, "", collapse = ",")))
dbDisconnect(con)

stage_label <- ifelse(meta$hotta_stage == "" | is.na(meta$hotta_stage), meta$stage, paste0(meta$stage, " (", meta$hotta_stage, ")"))
stage_names <- unique(stage_label)       # in developmental order
write.csv(data.frame(gene_id = gene_id, planted = plant, base = round(base, 2),
                     stage_label = ifelse(plant > 100, stage_names[pmax(plant - 100, 1)], NA)),
          file.path(out, "truth.csv"), row.names = FALSE)
write.csv(data.frame(planted = rep(seq_along(planted_go), lengths(planted_go)), go_id = unlist(planted_go)),
          file.path(out, "truth_go.csv"), row.names = FALSE)
cat("Wrote", G, "genes x", S, "samples to", out, "\n")
