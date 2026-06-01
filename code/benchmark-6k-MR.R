#!/usr/bin/env Rscript
# =============================================================================
# benchmark-6k-MR.R
#
# Evaluate MR hits against the tested-cell test set.
#
# Design principle
# ----------------
# The exam (test/test-6k-MR.tsv) is the set of (gene, disease) cells that MR
# actually tested -- one row per tested cell, with LABEL = 1 iff the cell is
# in universe-6k. There are no Cartesian / unmeasurable rows; every LABEL = 0
# row is a genuine "tested but not in universe" negative.
#
# Pipeline:
#   1. Load exam and MR results (MR.RData); apply the same QC filter
#      (drop exception in {PQTL_EMPTY, PQTL_NO_SIGNAL}) used by create-test.
#   2. Translate dat$protein (OlinkID) -> gene symbol via lookup-Olink.tsv.
#   3. Shrink dat to (gene, TRAIT) pairs in the exam, then FDR within scope.
#   4. q = p.adjust(p, "fdr"); hit iff q < 0.05. Non-numeric p strings -> NA -> 1.
#   5. Each exam cell is "predicted positive" if ANY of its constituent
#      (gene, TRAIT) rows is a hit. Build 2x2 vs LABEL and report enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-6k-MR.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-6k-MR.tsv", sep = "\t", quote = "")

load("result-compiled/MR.RData")                 # provides `dat`
setDT(dat)
# Mirror create-test QC: keep only rows that yielded a usable MR estimate.
dat <- dat[exception %in% c("blank", "ONLY_ONE_IV", "SUMSTAT_NO_SIGNAL")]

olink <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Translate protein (OlinkID) -> gene symbol
# -----------------------------------------------------------------------------
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

dat <- merge(dat, prot2gene, by = "protein")

# -----------------------------------------------------------------------------
# Explode exam TRAIT into (gene, disease, TRAIT) long form.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, `Disease indication`, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene, TRAIT) pairs in the exam.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(`Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs,
                 by = c("Target gene name", "TRAIT"))

# -----------------------------------------------------------------------------
# p -> numeric (NAs from coercion treated as 1), then FDR.
# -----------------------------------------------------------------------------
dat.sub[, p.num := suppressWarnings(as.numeric(p))]
dat.sub[is.na(p.num), p.num := 1]
dat.sub[, q := p.adjust(p.num, method = "fdr")]
dat.sub[, sig := q < 0.05]

# -----------------------------------------------------------------------------
# Propagate significance to exam cells.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[sig == TRUE, .(`Target gene name`, TRAIT)])
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

cat("--- MR Enrichment ---\n")
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
