#!/usr/bin/env Rscript
# =============================================================================
# benchmark-6k-TWAS.R
#
# Evaluate TWAS hits against the tested-cell test set.
#
# Design principle
# ----------------
# The exam (test/test-6k-TWAS.tsv) is the set of (gene, disease) cells that
# TWAS actually tested -- one row per tested cell, with LABEL = 1 iff the cell
# is in universe-6k. There are no Cartesian / unmeasurable rows; every
# LABEL = 0 row is a genuine "tested but not in universe" negative.
#
# Pipeline:
#   1. Load exam and TWAS results (TWAS.RData); apply the same QC filter
#      (pred_perf_r2 >= 0.01) used by create-test-6k-TWAS.R.
#   2. Shrink dat to (gene_name, TRAIT) pairs in the exam.
#   3. q = p.adjust(pvalue, "fdr") on the shrunken object; hit iff q < 0.05.
#   4. Each exam cell is "predicted positive" if ANY of its constituent
#      (gene, TRAIT) rows is a hit. Build 2x2 vs LABEL and report enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-6k-TWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-6k-TWAS.tsv", sep = "\t", quote = "")

load("result-compiled/TWAS.RData")               # provides `dat`
setDT(dat)
dat <- dat[pred_perf_r2 >= 0.01]                 # mirror create-test QC

# -----------------------------------------------------------------------------
# Explode exam TRAIT into (gene, disease, TRAIT) long form.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, `Disease indication`, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene_name, TRAIT) pairs in the exam.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(gene_name = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs, by = c("gene_name", "TRAIT"))

# -----------------------------------------------------------------------------
# FDR within the tested scope; q < 0.05 significant.
# -----------------------------------------------------------------------------
dat.sub[, q := p.adjust(pvalue, method = "fdr")]
dat.sub[, sig := !is.na(q) & q < 0.05]

# -----------------------------------------------------------------------------
# Propagate significance to exam cells.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[sig == TRUE,
                            .(`Target gene name` = gene_name, TRAIT)])
sig.pairs[, sig := TRUE]

exam.scored <- merge(
  exam.long, sig.pairs,
  by = c("Target gene name", "TRAIT"),
  all.x = TRUE
)
exam.scored[is.na(sig), sig := FALSE]

pred.pos.cells <- unique(
  exam.scored[sig == TRUE, .(`Target gene name`, `Disease indication`)]
)
pred.pos.cells[, predicted := 1L]

exam.pred <- merge(
  exam, pred.pos.cells,
  by = c("Target gene name", "Disease indication"),
  all.x = TRUE
)
exam.pred[is.na(predicted), predicted := 0L]

# -----------------------------------------------------------------------------
# Contingency table and enrichment (over the tested scope only -- no Cartesian)
# -----------------------------------------------------------------------------
TP <- exam.pred[LABEL == 1 & predicted == 1, .N]
FP <- exam.pred[LABEL == 0 & predicted == 1, .N]
FN <- exam.pred[LABEL == 1 & predicted == 0, .N]
TN <- exam.pred[LABEL == 0 & predicted == 0, .N]

base.rate     <- (TP + FN) / (TP + FP + FN + TN)
observed.rate <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
fold.enrich   <- observed.rate / base.rate

ft <- fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE))

cat("--- TWAS Enrichment ---\n")
cat(sprintf("Contingency table:  TP=%d  FP=%d  FN=%d  TN=%d\n",
            TP, FP, FN, TN))
cat(sprintf("Base rate:          %.4f\n", base.rate))
cat(sprintf("Observed rate:      %.4f\n", observed.rate))
cat(sprintf("Fold enrichment:    %.2f\n", fold.enrich))
cat(sprintf("Odds ratio:         %.2f\n", ft$estimate))
cat(sprintf("95%% CI:             [%.2f, %.2f]\n",
            ft$conf.int[1], ft$conf.int[2]))
cat(sprintf("Fisher p-value:     %.2e\n", ft$p.value))

# -----------------------------------------------------------------------------
# True positives: tested + hit + in universe
# -----------------------------------------------------------------------------
tp.cells <- exam.pred[LABEL == 1 & predicted == 1,
                      .(`Target gene name`, `Disease indication`)]
setorder(tp.cells, `Target gene name`, `Disease indication`)
cat("\n--- True positives ---\n")
print(tp.cells, nrows = Inf)
