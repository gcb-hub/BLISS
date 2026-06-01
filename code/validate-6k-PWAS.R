#!/usr/bin/env Rscript
# =============================================================================
# validate-6k-PWAS.R
#
# Report two numbers for PWAS:
#   Number 1 (denominator): rows in data/universe-6k.tsv that PWAS can
#                            reasonably answer.
#   Number 2 (numerator):   rows in the denominator that are actually
#                            validated by the significant portion of PWAS
#                            results (i.e. at least one of the row's candidate
#                            TRAITs is significant for its gene).
#
# Pipeline:
#   1. Load test/denom-6k-PWAS.tsv (produced by create-denom-6k-PWAS.R).
#   2. Load PWAS.RData and apply the same QC filter used to build the
#      denominator (r2.c >= 0.01).
#   3. Shrink dat to (gene, TRAIT) pairs present in the denominator, then
#      compute FDR within that scope. Hit iff q < 0.05.
#   4. Mark each denominator row as validated iff any of its candidate
#      (gene, TRAIT) pairs is a hit.
#
# Outputs:
#   stdout: N1, N2, proportion
#   test/validated-6k-PWAS.tsv -- denominator file plus a `validated` flag
#
# Run from the project root:
#   Rscript code/validate-6k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
denom <- fread("test/denom-6k-PWAS.tsv", sep = "\t", quote = "")

load("result-compiled/PWAS.RData")   # provides `dat`
setDT(dat)
dat <- dat[r2.c >= 0.01]             # mirror create-denom QC

# -----------------------------------------------------------------------------
# Explode denom TRAIT into long (row, gene, TRAIT) form for the scope / join.
# -----------------------------------------------------------------------------
denom[, denom_row := .I]
denom.long <- denom[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(denom_row, `Target gene name`)]
denom.long[, TRAIT := trimws(TRAIT)]
denom.long <- denom.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene, TRAIT) pairs in the denominator and apply FDR.
# -----------------------------------------------------------------------------
scope.pairs <- unique(denom.long[, .(ID.inner = `Target gene name`, TRAIT)])
dat.sub <- merge(dat, scope.pairs, by = c("ID.inner", "TRAIT"))

dat.sub[, q := p.adjust(p.a.c, method = "fdr")]
dat.sub[, sig := !is.na(q) & q < 0.05]

sig.pairs <- unique(dat.sub[sig == TRUE,
                            .(`Target gene name` = ID.inner, TRAIT)])
sig.pairs[, sig := TRUE]

# -----------------------------------------------------------------------------
# A denom row is validated iff any of its candidate (gene, TRAIT) pairs is sig.
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
fwrite(out, "test/validated-6k-PWAS.tsv", sep = "\t", quote = FALSE)

cat("--- PWAS validation ---\n")
cat(sprintf("Denominator (N1): %d\n", N1))
cat(sprintf("Numerator   (N2): %d\n", N2))
cat(sprintf("Proportion N2/N1: %.4f\n",
            if (N1 > 0) N2 / N1 else NA_real_))
