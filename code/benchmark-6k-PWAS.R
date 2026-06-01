#!/usr/bin/env Rscript
# =============================================================================
# benchmark-6k-PWAS.R
#
# Evaluate PWAS hits against the tested-cell test set.
#
# Design principle
# ----------------
# The exam (test/test-6k-PWAS.tsv) is the set of (gene, disease) cells that
# PWAS actually tested -- one row per tested cell, with LABEL = 1 iff the cell
# appears in universe-6k. There are no Cartesian / unmeasurable rows; every
# LABEL = 0 row is a genuine "tested but not in universe" negative.
#
# Pipeline:
#   1. Load exam and PWAS results (PWAS.RData); apply the same QC filter
#      (r2.c >= 0.01) used by create-test-6k-PWAS.R.
#   2. Shrink `dat` down to rows whose (ID.inner, TRAIT) pair appears in the
#      exam's TRAIT column (exploded on ";"). This restricts the FDR universe
#      to the tested scope used in the exam.
#   3. q = p.adjust(p.a.c, method = "fdr") on the shrunken object;
#      hit iff q < 0.05.
#   4. Each exam cell is "predicted positive" if ANY of its constituent
#      (gene, TRAIT) rows is a hit.
#   5. Build the 2x2 against exam LABEL and report enrichment +
#      Fisher exact test.
#
# Run from the project root:
#   Rscript code/benchmark-6k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-6k-PWAS.tsv", sep = "\t", quote = "")

load("result-compiled/PWAS.RData")   # provides `dat`
setDT(dat)
dat <- dat[r2.c >= 0.01]             # mirror create-test QC

# -----------------------------------------------------------------------------
# Explode exam TRAIT into (gene, disease, TRAIT) long form.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, `Disease indication`, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene, TRAIT) pairs in the exam, then FDR within that scope.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(ID.inner = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs, by = c("ID.inner", "TRAIT"))

# -- significance threshold: uncomment ONE line below --
 sig.method <- "fdr"; cutoff <- 0.05
# sig.method <- "fdr"; cutoff <- 0.10
# sig.method <- "fdr"; cutoff <- 0.01
# sig.method <- "fdr"; cutoff <- 0.001
# sig.method <- "bonferroni"; cutoff <- 0.05
# sig.method <- "bonferroni"; cutoff <- 0.01

dat.sub[, q := p.adjust(p.a.c, method = sig.method)]
dat.sub[, sig := !is.na(q) & q < cutoff]

# -----------------------------------------------------------------------------
# Propagate significance to exam cells: a cell is predicted positive if ANY
# of its constituent (gene, TRAIT) rows is significant.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[sig == TRUE,
                            .(`Target gene name` = ID.inner, TRAIT)])
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
# Contingency and enrichment
#
# All four cells of the 2x2 sit inside the tested scope:
#   TP = tested, hit,     in universe
#   FP = tested, hit,     not in universe
#   FN = tested, not hit, in universe
#   TN = tested, not hit, not in universe
# Base rate = (TP + FN) / N_tested  -- the "background" universe-density among
# every cell PWAS actually tested. Observed rate = TP / (TP + FP). Their ratio
# is the fold enrichment among hits.
# -----------------------------------------------------------------------------
TP <- exam.pred[LABEL == 1 & predicted == 1, .N]
FP <- exam.pred[LABEL == 0 & predicted == 1, .N]
FN <- exam.pred[LABEL == 1 & predicted == 0, .N]
TN <- exam.pred[LABEL == 0 & predicted == 0, .N]

base.rate     <- (TP + FN) / (TP + FP + FN + TN)
observed.rate <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
fold.enrich   <- observed.rate / base.rate

ft <- fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE))

cat(sprintf("--- PWAS Enrichment (%s < %g) ---\n", sig.method, cutoff))
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
