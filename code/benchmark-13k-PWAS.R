#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-PWAS.R
#
# Evaluate PWAS results against the 13k gene x MeSH-disease "exam".
#
#   1. Load exam (test/test-13k-PWAS.tsv) and PWAS results (PWAS.RData).
#   2. Shrink `dat` down to rows whose (ID.inner, TRAIT) pair appears in the
#      exam's TRAIT column (exploded on ";").
#   3. Recompute FDR via p.adjust(dat$p.a.c, method = "fdr") on the shrunken
#      object; call q < 0.05 significant.
#   4. For each exam row, mark it "predicted positive" if ANY shrunken dat row
#      matching its (gene, TRAIT-set) is significant.
#   5. Build the 2x2 contingency table against exam LABEL and report
#      enrichment metrics + Fisher exact test.
#
# Run from the project root:
#   Rscript code/benchmark-13k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-13k-PWAS.tsv", sep = "\t", quote = "")

load("result-compiled/PWAS.RData")   # provides `dat`
setDT(dat)

# -----------------------------------------------------------------------------
# Explode exam TRAIT (semicolon-separated) into long form: one row per
# (gene, MESH, TRAIT). Rows with empty TRAIT cannot match any dat row and are
# dropped here -- they will fall through as predicted=0 in the join below.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, MESH, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat: keep only rows whose (ID.inner, TRAIT) pair appears in exam.long.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(ID.inner = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs, by = c("ID.inner", "TRAIT"))

# -----------------------------------------------------------------------------
# FDR on the shrunken object; significance threshold q < 0.05.
# -----------------------------------------------------------------------------
dat.sub[, q := p.adjust(p.a.c, method = "fdr")]
dat.sub[, sig := !is.na(q) & q < 0.05]

# -----------------------------------------------------------------------------
# Join significance back onto the long exam, then collapse to (gene, MESH):
# a cell is predicted positive if ANY of its constituent (gene, TRAIT) rows
# is significant.
# -----------------------------------------------------------------------------
sig.pairs <- unique(dat.sub[sig == TRUE, .(`Target gene name` = ID.inner, TRAIT)])
sig.pairs[, sig := TRUE]

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
# Contingency and enrichment
# -----------------------------------------------------------------------------
TP <- exam.pred[LABEL == 1 & predicted == 1, .N]
FP <- exam.pred[LABEL == 0 & predicted == 1, .N]
FN <- exam.pred[LABEL == 1 & predicted == 0, .N]
TN <- exam.pred[LABEL == 0 & predicted == 0, .N]

base.rate     <- (TP + FN) / (TP + FP + FN + TN)
observed.rate <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
fold.enrich   <- observed.rate / base.rate

ft <- fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE))

cat("--- PWAS Enrichment ---\n")
cat(sprintf("Contingency table:  TP=%d  FP=%d  FN=%d  TN=%d\n",
            TP, FP, FN, TN))
cat(sprintf("Base rate:          %.4f\n", base.rate))
cat(sprintf("Observed rate:      %.4f\n", observed.rate))
cat(sprintf("Fold enrichment:    %.2f\n", fold.enrich))
cat(sprintf("Odds ratio:         %.2f\n", ft$estimate))
cat(sprintf("95%% CI:             [%.2f, %.2f]\n",
            ft$conf.int[1], ft$conf.int[2]))
cat(sprintf("Fisher p-value:     %.2e\n", ft$p.value))
