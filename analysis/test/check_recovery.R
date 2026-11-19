# Checks that the notebook recovers what make_test_data.R planted.
# Usage: Rscript check_recovery.R <test_data_dir> <analysis_output_dir>
args <- commandArgs(trailingOnly = TRUE)
td <- args[1]; ao <- args[2]
truth <- read.csv(file.path(td, "truth.csv"))
tgo   <- read.csv(file.path(td, "truth_go.csv"))
gs    <- read.csv(file.path(ao, "gene_summary.csv"))
enr   <- read.csv(file.path(ao, "cluster_go_enrichment.csv"))
d <- merge(truth, gs, by = "gene_id")

ok <- TRUE
check <- function(name, pass, detail) {
  cat(sprintf("%-4s %s: %s\n", if (pass) "PASS" else "FAIL", name, detail)); if (!pass) ok <<- FALSE
}
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) x * (x - 1) / 2
  sij <- sum(c2(tab)); si <- sum(c2(rowSums(tab))); sj <- sum(c2(colSums(tab))); ex <- si * sj / c2(n)
  (sij - ex) / (0.5 * (si + sj) - ex)
}

pl <- d[d$planted %in% 1:6, ]
r_hi <- mean(pl$temporal[pl$base >= 20], na.rm = TRUE)
check("recall of planted genes (base mean >= 20)", r_hi >= 0.90, sprintf("%.1f%% called temporal", 100 * r_hi))
check("recall of all planted genes", mean(pl$temporal, na.rm = TRUE) >= 0.75,
      sprintf("%.1f%% (low-count genes are expected to be missed)", 100 * mean(pl$temporal, na.rm = TRUE)))
fl <- d[d$planted == 0 & d$tested, ]
spec <- d[d$planted > 100, ]
called <- d[!is.na(d$temporal) & d$temporal, ]
fdp <- mean(called$planted %in% c(0, -1))
check("false discovery proportion among temporal genes", fdp <= 0.05,
      sprintf("%.1f%% of %d called genes are flat genes (nominal FDR 1%%); %.2f%% of the %d tested flat genes were called",
              100 * fdp, nrow(called), 100 * mean(fl$temporal, na.rm = TRUE), nrow(fl)))
lo <- d[d$planted == -1, ]
check("very-low genes are filtered out", mean(lo$tested) < 0.5, sprintf("%.0f%% of them tested", 100 * mean(lo$tested)))

pt <- pl[!is.na(pl$cluster), ]
k  <- length(unique(gs$cluster[!is.na(gs$cluster)]))
# purity: share of planted genes whose cluster is dominated by their own pattern (splitting a pattern over two clusters is
# harmless, merging different patterns is not). Adjusted Rand index is printed for information only.
maj <- tapply(pt$planted, pt$cluster, function(x) as.integer(names(which.max(table(x)))))
purity <- mean(maj[as.character(pt$cluster)] == pt$planted)
check("cluster purity", purity >= 0.9, sprintf("%.1f%% of planted genes sit in a cluster dominated by their own pattern (k = %d; ARI = %.2f)",
      100 * purity, k, ari(pt$planted, pt$cluster)))
check("no pattern is merged away", all(1:6 %in% maj), sprintf("patterns that dominate at least one cluster: %s", paste(sort(unique(maj)), collapse = ",")))

# stage markers: genes expressed in one stage only must come out with that peak stage and a high tau
mk <- read.csv(file.path(ao, "stage_markers.csv"))
check("stage-specific genes: peak stage correct", mean(spec$peak_stage == spec$stage_label) >= 0.95,
      sprintf("%d of %d planted genes", sum(spec$peak_stage == spec$stage_label), nrow(spec)))
good <- !is.na(spec$tau) & spec$tau >= 0.8 & spec$max_stage_tpm >= 5 & spec$peak_stage == spec$stage_label
check("stage-specific genes: right peak stage, tau >= 0.8", mean(good) >= 0.9,
      sprintf("%d of %d planted genes", sum(good), nrow(spec)))
inmk <- spec$gene_id %in% mk$gene_id
check("stage-specific genes listed in stage_markers.csv", mean(inmk) >= 0.9, sprintf("%d of %d", sum(inmk), nrow(spec)))
check("no flat gene is listed as marker", !any(mk$gene_id %in% d$gene_id[d$planted == 0]), sprintf("%d markers listed in total", nrow(mk)))

# enrichment: planted GO terms should be significant in the cluster that holds most of that pattern's genes
hit <- 0
for (p in 1:6) {
  best <- as.integer(names(which.max(table(pt$cluster[pt$planted == p]))))
  for (go in tgo$go_id[tgo$planted == p]) {
    s <- enr[enr$cluster == best & enr$go_id == go, ]
    pass <- nrow(s) == 1 && s$padj < 0.05
    hit <- hit + pass
    cat(sprintf("       pattern %d -> cluster %d, %s: %s\n", p, best, go,
                if (nrow(s)) sprintf("padj = %.2g, fold = %.1f", s$padj, s$fold) else "not tested"))
  }
}
check("planted GO terms found", hit / nrow(tgo) >= 0.8, sprintf("%d of %d", hit, nrow(tgo)))
cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
