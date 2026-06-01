#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-MR.R
#
# Evaluate MR results against the 13k gene x MeSH-disease "exam".
#
#   1. Load exam (test/test-13k-MR.tsv) and MR results (MR.RData).
#   2. Translate dat$protein (OlinkID) -> gene symbol via lookup-Olink.tsv.
#   3. Shrink `dat` to rows whose (gene, TRAIT) pair appears in the exam
#      (TRAIT column exploded on ";").
#   4. Coerce dat$p to numeric (NAs -> 1), then q = p.adjust(p, method = "fdr");
#      q < 0.05 is significant.
#   5. For each exam (gene, MESH) cell, predicted positive if ANY matching
#      shrunken dat row is significant. Report 2x2 + Fisher enrichment.
#
# Run from the project root:
#   Rscript code/benchmark-13k-MR.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
exam <- fread("test/test-13k-MR.tsv", sep = "\t", quote = "")

load("result-compiled/MR.RData")                  # provides `dat`
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
# p -> numeric (coercion NAs treated as 1), then FDR; sig := q < 0.05.
# -----------------------------------------------------------------------------
dat.sub[, p.num := suppressWarnings(as.numeric(p))]
dat.sub[is.na(p.num), p.num := 1]
dat.sub[, q := p.adjust(p.num, method = "fdr")]
dat.sub[, sig := q < 0.05]

# -----------------------------------------------------------------------------
# Propagate significance to exam (gene, MESH) cells: predicted positive if ANY
# constituent (gene, TRAIT) in dat.sub is significant.
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
