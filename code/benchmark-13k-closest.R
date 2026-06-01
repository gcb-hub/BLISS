#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-closest.R
#
# Evaluate the "closest gene" baseline against the 13k gene x MeSH-disease
# "exam". For each MVP lead variant, the prediction is the
# `Lead Variant Gene (Nearest)` column from
# data/loci-from-MVP-corrected.csv (EUR rows only) -- there
# is no significance threshold, so every (gene, TRAIT) pair from MVP that
# falls inside the exam's tested scope counts as a positive prediction.
#
#   1. Load exam (test/test-13k-closest.tsv) and MVP loci.
#   2. Synthesize `dat` from MVP: (Target gene name, TRAIT) where the gene is
#      the nearest gene to a lead variant. No protein->gene translation.
#   3. Shrink `dat` to rows whose (gene, TRAIT) pair appears in the exam
#      (TRAIT column exploded on ";").
#   4. No threshold: every row in the shrunken dat is "significant".
#   5. For each exam (gene, MESH) cell, predicted positive if ANY matching
#      shrunken dat row is significant. Report 2x2 + Fisher enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-13k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-13k-closest.tsv", sep = "\t", quote = "")

loci <- fread("data/loci-from-MVP-corrected.csv")
if (colnames(loci)[1] != "Population") setnames(loci, 1, "Population")
loci <- loci[Population == "EUR"]

# -----------------------------------------------------------------------------
# Synthesize `dat` from MVP loci.
# -----------------------------------------------------------------------------
loci[, p.num := suppressWarnings(as.numeric(`P-Value`))]
dat <- unique(loci[
  nzchar(`Lead Variant Gene (Nearest)`) &
  !is.na(`Lead Variant Gene (Nearest)`) &
  nzchar(Trait) & !is.na(Trait) &
  !is.na(p.num) & p.num < 5e-8,
  .(`Target gene name` = `Lead Variant Gene (Nearest)`,
    TRAIT              = Trait)
])

# -----------------------------------------------------------------------------
# Explode exam TRAIT (semicolon-separated) into (gene, MESH, TRAIT) long form.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, MESH, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene, TRAIT) pairs present in the exam.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(`Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs,
                 by = c("Target gene name", "TRAIT"))

# -----------------------------------------------------------------------------
# No threshold: every surviving (gene, TRAIT) pair is a positive prediction.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[, .(`Target gene name`, TRAIT)])
sig.pairs[, sig := TRUE]

# -----------------------------------------------------------------------------
# Propagate significance to exam (gene, MESH) cells.
# -----------------------------------------------------------------------------
exam.scored <- merge(
  exam.long, sig.pairs,
  by = c("Target gene name", "TRAIT"),
  all.x = TRUE
)
exam.scored[is.na(sig), sig := FALSE]

pred.pos.cells <- unique(
  exam.scored[sig == TRUE, .(`Target gene name`, MESH)]
)
pred.pos.cells[, predicted := 1L]

exam.pred <- merge(
  exam, pred.pos.cells,
  by = c("Target gene name", "MESH"),
  all.x = TRUE
)
exam.pred[is.na(predicted), predicted := 0L]

# -----------------------------------------------------------------------------
# Contingency table and enrichment
# -----------------------------------------------------------------------------
TP <- exam.pred[LABEL == 1 & predicted == 1, .N]
FP <- exam.pred[LABEL == 0 & predicted == 1, .N]
FN <- exam.pred[LABEL == 1 & predicted == 0, .N]
TN <- exam.pred[LABEL == 0 & predicted == 0, .N]

base.rate     <- (TP + FN) / (TP + FP + FN + TN)
observed.rate <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
fold.enrich   <- observed.rate / base.rate

ft <- fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE))

cat("--- closest Enrichment ---\n")
cat(sprintf("Contingency table:  TP=%d  FP=%d  FN=%d  TN=%d\n",
            TP, FP, FN, TN))
cat(sprintf("Base rate:          %.4f\n", base.rate))
cat(sprintf("Observed rate:      %.4f\n", observed.rate))
cat(sprintf("Fold enrichment:    %.2f\n", fold.enrich))
cat(sprintf("Odds ratio:         %.2f\n", ft$estimate))
cat(sprintf("95%% CI:             [%.2f, %.2f]\n",
            ft$conf.int[1], ft$conf.int[2]))
cat(sprintf("Fisher p-value:     %.2e\n", ft$p.value))
