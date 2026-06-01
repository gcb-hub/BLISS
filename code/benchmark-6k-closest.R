#!/usr/bin/env Rscript
# =============================================================================
# benchmark-6k-closest.R
#
# Evaluate the "closest gene" baseline against the 6k tested-cell test set.
# For each MVP lead variant, the prediction is the
# `Lead Variant Gene (Nearest)` column from data/loci-from-MVP.csv -- there
# is no significance threshold, so every (gene, TRAIT) pair from MVP that
# falls inside the exam's tested scope counts as a positive prediction.
#
# Pipeline:
#   1. Load exam (test/test-6k-closest.tsv) and MVP loci.
#   2. Synthesize `dat` from MVP: (Target gene name, TRAIT) where the gene is
#      the nearest gene to a lead variant. No protein->gene translation,
#      no QC filter beyond dropping empty/NA gene/Trait.
#   3. Shrink `dat` to rows whose (gene, TRAIT) pair appears in the exam
#      (TRAIT column exploded on ";").
#   4. No threshold: every row in the shrunken dat is "significant".
#   5. For each exam (gene, Disease indication) cell, predicted positive if
#      ANY matching shrunken dat row is significant. Build 2x2 + Fisher.
#
# Known degeneracy
# ----------------
# Because the exam was built FROM MVP (gene, TRAIT) pairs and every such
# pair is "sig=TRUE", every exam cell becomes predicted=1 -- so FN=TN=0,
# Fisher's OR is Inf, and fold enrichment = 1.0. This is an inherent
# property of closest under the tested-cells framework, not a bug. The
# meaningful comparison for closest is across methods (its base rate vs
# coloc/MR/PWAS/TWAS), or via validate-6k-closest.R.
#
# Run from the project root:
#   Rscript code/benchmark-6k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-6k-closest.tsv", sep = "\t", quote = "")

loci <- fread("data/loci-from-MVP.csv")
if (colnames(loci)[1] != "Trait") setnames(loci, 1, "Trait")

# -----------------------------------------------------------------------------
# Synthesize `dat` from MVP loci.
# -----------------------------------------------------------------------------
loci[, p.num := suppressWarnings(as.numeric(`P-value`))]
dat <- unique(loci[
  nzchar(`Lead Variant Gene (Nearest)`) &
  !is.na(`Lead Variant Gene (Nearest)`) &
  nzchar(Trait) & !is.na(Trait) &
  !is.na(p.num) & p.num < 5e-8,
  .(`Target gene name` = `Lead Variant Gene (Nearest)`,
    TRAIT              = Trait)
])

# -----------------------------------------------------------------------------
# Explode exam TRAIT (semicolon-separated) into long form.
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
# No threshold: every surviving (gene, TRAIT) pair is a positive prediction.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[, .(`Target gene name`, TRAIT)])
sig.pairs[, sig := TRUE]

# -----------------------------------------------------------------------------
# Propagate significance to exam cells.
# -----------------------------------------------------------------------------
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

ft <- tryCatch(
  fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE)),
  error = function(e) NULL
)

cat("--- closest Enrichment ---\n")
cat(sprintf("Contingency table:  TP=%d  FP=%d  FN=%d  TN=%d\n",
            TP, FP, FN, TN))
cat(sprintf("Base rate:          %.4f\n", base.rate))
cat(sprintf("Observed rate:      %.4f\n", observed.rate))
cat(sprintf("Fold enrichment:    %.2f\n", fold.enrich))
if (is.null(ft) || (FN == 0 && TN == 0)) {
  cat("Odds ratio:         N/A (degenerate 2x2)\n")
  cat("95% CI:             N/A (degenerate 2x2)\n")
  cat("Fisher p-value:     N/A (degenerate 2x2)\n")
} else {
  cat(sprintf("Odds ratio:         %.2f\n", ft$estimate))
  cat(sprintf("95%% CI:             [%.2f, %.2f]\n",
              ft$conf.int[1], ft$conf.int[2]))
  cat(sprintf("Fisher p-value:     %.2e\n", ft$p.value))
}

# -----------------------------------------------------------------------------
# True positives: tested + hit + in universe
# -----------------------------------------------------------------------------
tp.cells <- exam.pred[LABEL == 1 & predicted == 1,
                      .(`Target gene name`, `Disease indication`)]
setorder(tp.cells, `Target gene name`, `Disease indication`)
cat("\n--- True positives ---\n")
print(tp.cells, nrows = Inf)
