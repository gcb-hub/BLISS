#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-TWAS.R
#
# Evaluate TWAS results against the 13k gene x MeSH-disease "exam".
#
#   1. Load exam (test/test-13k-TWAS.tsv) and TWAS results (TWAS.RData).
#   2. Shrink `dat` to rows whose (gene_name, TRAIT) pair appears in the exam
#      (TRAIT column exploded on ";").
#   3. q = p.adjust(dat$pvalue, method = "fdr") on the shrunken object;
#      q < 0.05 is significant.
#   4. For each exam (gene, MESH) cell, predicted positive if ANY matching
#      shrunken dat row is significant. Report 2x2 + Fisher enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-13k-TWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-13k-TWAS.tsv", sep = "\t", quote = "")

load("result-compiled/TWAS.RData")               # provides `dat`
setDT(dat)

# -----------------------------------------------------------------------------
# Explode exam TRAIT (semicolon-separated) into (gene, MESH, TRAIT) long form.
# -----------------------------------------------------------------------------
exam.long <- exam[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(`Target gene name`, MESH, LABEL)]
exam.long[, TRAIT := trimws(TRAIT)]
exam.long <- exam.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene_name, TRAIT) pairs present in the exam.
# -----------------------------------------------------------------------------
exam.pairs <- unique(exam.long[, .(gene_name = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, exam.pairs, by = c("gene_name", "TRAIT"))

# -----------------------------------------------------------------------------
# FDR on the shrunken object; q < 0.05 significant.
# -----------------------------------------------------------------------------
dat.sub[, q := p.adjust(pvalue, method = "fdr")]
dat.sub[, sig := !is.na(q) & q < 0.05]

# -----------------------------------------------------------------------------
# Propagate significance to exam (gene, MESH) cells.
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
