#!/usr/bin/env Rscript
# =============================================================================
# validate-6k-TWAS.R
#
# Report two numbers for TWAS:
#   Number 1 (denominator): rows in data/universe-6k.tsv that TWAS can
#                            reasonably answer.
#   Number 2 (numerator):   rows in the denominator validated by the
#                            significant portion of TWAS results.
#
# Pipeline:
#   1. Load test/denom-6k-TWAS.tsv (from create-denom-6k-TWAS.R).
#   2. Load TWAS.RData, apply the same QC filter (pred_perf_r2 >= 0.01).
#   3. Shrink dat to (gene_name, TRAIT) pairs in the denominator, FDR within
#      that scope; hit iff q < 0.05.
#   4. Mark each denominator row as validated iff any of its candidate
#      (gene, TRAIT) pairs is a hit.
#
# Outputs:
#   stdout: N1, N2, proportion
#   test/validated-6k-TWAS.tsv -- denominator file plus a `validated` flag
#
# Run from the project root:
#   Rscript code/validate-6k-TWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
denom <- fread("test/denom-6k-TWAS.tsv", sep = "\t", quote = "")

load("result-compiled/TWAS.RData")               # provides `dat`
setDT(dat)
dat <- dat[pred_perf_r2 >= 0.01]                 # mirror create-denom QC

# -----------------------------------------------------------------------------
# Explode denom TRAIT into long (row, gene, TRAIT) form.
# -----------------------------------------------------------------------------
denom[, denom_row := .I]
denom.long <- denom[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(denom_row, `Target gene name`)]
denom.long[, TRAIT := trimws(TRAIT)]
denom.long <- denom.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene_name, TRAIT) in the denominator; FDR within scope.
# -----------------------------------------------------------------------------
scope.pairs <- unique(denom.long[, .(gene_name = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, scope.pairs, by = c("gene_name", "TRAIT"))

dat.sub[, q := p.adjust(pvalue, method = "fdr")]
dat.sub[, sig := !is.na(q) & q < 0.05]

sig.pairs <- unique(dat.sub[sig == TRUE,
                            .(`Target gene name` = gene_name, TRAIT)])
sig.pairs[, sig := TRUE]

# -----------------------------------------------------------------------------
# Mark validated rows
# -----------------------------------------------------------------------------
denom.scored <- merge(denom.long, sig.pairs,
                      by = c("Target gene name", "TRAIT"),
                      all.x = TRUE)
denom.scored[is.na(sig), sig := FALSE]

validated_rows <- unique(denom.scored[sig == TRUE, denom_row])
denom[, validated := as.integer(denom_row %in% validated_rows)]

N1 <- nrow(denom)
N2 <- sum(denom$validated)

# -----------------------------------------------------------------------------
# Save and report
# -----------------------------------------------------------------------------
out <- copy(denom)
out[, denom_row := NULL]
fwrite(out, "test/validated-6k-TWAS.tsv", sep = "\t", quote = FALSE)

cat("--- TWAS validation ---\n")
cat(sprintf("Denominator (N1): %d\n", N1))
cat(sprintf("Numerator   (N2): %d\n", N2))
cat(sprintf("Proportion N2/N1: %.4f\n",
            if (N1 > 0) N2 / N1 else NA_real_))
