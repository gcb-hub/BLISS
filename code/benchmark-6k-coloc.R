#!/usr/bin/env Rscript
# =============================================================================
# benchmark-6k-coloc.R
#
# Evaluate coloc hits against the tested-cell test set.
#
# Design principle
# ----------------
# The exam (test/test-6k-coloc.tsv) is the set of (gene, disease) cells that
# coloc actually tested -- one row per tested cell, with LABEL = 1 iff the
# cell is in universe-6k. There are no Cartesian / unmeasurable rows; every
# LABEL = 0 row is a genuine "tested but not in universe" negative.
#
# Pipeline:
#   1. Load exam and coloc results (coloc.RData). No row-level QC for coloc.
#   2. Translate dat$protein (OlinkID) -> gene symbol via lookup-Olink.tsv.
#   3. Shrink dat to (gene, TRAIT) pairs in the exam.
#   4. Hit iff H4 >= 0.8 (non-numeric H4 -> 0).
#   5. Each exam cell is "predicted positive" if ANY of its constituent
#      (gene, TRAIT) rows is a hit. Build 2x2 vs LABEL and report enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-6k-coloc.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-6k-coloc.tsv", sep = "\t", quote = "")

load("result-compiled/coloc.RData")              # provides `dat`
setDT(dat)

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
# H4 -> numeric (NAs from coercion treated as 0); threshold at 0.8.
# -----------------------------------------------------------------------------
dat.sub[, H4.num := suppressWarnings(as.numeric(H4))]
dat.sub[is.na(H4.num), H4.num := 0]
dat.sub[, sig := H4.num >= 0.8]

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

cat("--- coloc Enrichment ---\n")
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
