#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-coloc.R
#
# Evaluate coloc results against the 13k gene x MeSH-disease "exam".
#
#   1. Load exam (test/test-13k-coloc.tsv) and coloc results (coloc.RData).
#   2. Translate dat$protein (OlinkID) -> gene symbol via lookup-Olink.tsv.
#   3. Shrink `dat` to rows whose (gene, TRAIT) pair appears in the exam
#      (TRAIT column exploded on ";").
#   4. Coerce dat$H4 to numeric (NAs from coercion treated as 0); significance
#      threshold is H4 >= 0.8.
#   5. For each exam (gene, MESH) cell, predicted positive if ANY matching
#      shrunken dat row is significant. Report 2x2 + Fisher enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-13k-coloc.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-13k-coloc.tsv", sep = "\t", quote = "")

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
# H4 -> numeric (NAs from coercion treated as 0); threshold at 0.8.
# -----------------------------------------------------------------------------
dat.sub[, H4.num := suppressWarnings(as.numeric(H4))]
dat.sub[is.na(H4.num), H4.num := 0]
dat.sub[, sig := H4.num >= 0.8]

# -----------------------------------------------------------------------------
# Propagate significance to exam (gene, MESH) cells.
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
