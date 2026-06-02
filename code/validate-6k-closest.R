#!/usr/bin/env Rscript
# =============================================================================
# validate-6k-closest.R
#
# Report two numbers for closest:
#   Number 1 (denominator): rows in data/universe-6k.tsv that closest can
#                            reasonably answer.
#   Number 2 (numerator):   rows in the denominator validated by the MVP
#                            (gene, TRAIT) predictions.
#
# Pipeline:
#   1. Load test/denom-6k-closest.tsv (from create-denom-6k-closest.R).
#   2. Load data/loci-from-MVP-corrected.csv and synthesize `dat` of MVP (gene, TRAIT)
#      pairs (drop empty/NA Lead Variant Gene (Nearest) and Trait).
#   3. Shrink dat to (gene, TRAIT) pairs in the denominator. No threshold:
#      every surviving pair is a hit (closest has no p-value).
#   4. Mark each denominator row as validated iff any of its candidate
#      (gene, TRAIT) pairs is a hit.
#
# Outputs:
#   stdout: N1, N2, proportion
#   test/validated-6k-closest.tsv -- denominator file plus a `validated` flag
#
# Run from the project root:
#   Rscript code/validate-6k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
denom <- fread("test/denom-6k-closest.tsv", sep = "\t", quote = "")

loci <- fread("data/loci-from-MVP-corrected.csv")

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
# Explode denom TRAIT into long (row, gene, TRAIT) form.
# -----------------------------------------------------------------------------
denom[, denom_row := .I]
denom.long <- denom[nzchar(TRAIT), .(
  TRAIT = unlist(strsplit(TRAIT, ";"))
), by = .(denom_row, `Target gene name`)]
denom.long[, TRAIT := trimws(TRAIT)]
denom.long <- denom.long[nzchar(TRAIT)]

# -----------------------------------------------------------------------------
# Shrink dat to (gene, TRAIT) in the denominator. No threshold.
# -----------------------------------------------------------------------------
scope.pairs <- unique(denom.long[, .(`Target gene name`, TRAIT)])
dat.sub <- merge(dat, scope.pairs, by = c("Target gene name", "TRAIT"))

sig.pairs <- unique(dat.sub[, .(`Target gene name`, TRAIT)])
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
fwrite(out, "test/validated-6k-closest.tsv", sep = "\t", quote = FALSE)

cat("--- closest validation ---\n")
cat(sprintf("Denominator (N1): %d\n", N1))
cat(sprintf("Numerator   (N2): %d\n", N2))
cat(sprintf("Proportion N2/N1: %.4f\n",
            if (N1 > 0) N2 / N1 else NA_real_))
