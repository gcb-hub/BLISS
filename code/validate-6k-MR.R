#!/usr/bin/env Rscript
# =============================================================================
# validate-6k-MR.R
#
# Report two numbers for MR:
#   Number 1 (denominator): rows in data/universe-6k.tsv that MR can
#                            reasonably answer.
#   Number 2 (numerator):   rows in the denominator validated by the
#                            significant portion of MR results.
#
# Pipeline:
#   1. Load test/denom-6k-MR.tsv (from create-denom-6k-MR.R).
#   2. Load MR.RData, apply the same QC filter (exception in blank /
#      ONLY_ONE_IV / SUMSTAT_NO_SIGNAL) and translate protein -> gene via
#      lookup-Olink.
#   3. Shrink dat to (gene, TRAIT) pairs in the denominator, FDR within that
#      scope, hit iff q < 0.05. Non-numeric p strings -> NA -> 1.
#   4. Mark each denominator row as validated iff any of its candidate
#      (gene, TRAIT) pairs is a hit.
#
# Outputs:
#   stdout: N1, N2, proportion
#   test/validated-6k-MR.tsv -- denominator file plus a `validated` flag
#
# Run from the project root:
#   Rscript code/validate-6k-MR.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
denom <- fread("test/denom-6k-MR.tsv", sep = "\t", quote = "")

load("result-compiled/MR.RData")                 # provides `dat`
setDT(dat)
dat <- dat[exception %in% c("blank", "ONLY_ONE_IV", "SUMSTAT_NO_SIGNAL")]

olink <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# protein (OlinkID) -> gene symbol
# -----------------------------------------------------------------------------
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

dat <- merge(dat, prot2gene, by = "protein")

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
# Shrink dat to (gene, TRAIT) in the denominator; FDR within scope.
# -----------------------------------------------------------------------------
scope.pairs <- unique(denom.long[, .(`Target gene name`, TRAIT)])
dat.sub <- merge(dat, scope.pairs, by = c("Target gene name", "TRAIT"))

dat.sub[, p.num := suppressWarnings(as.numeric(p))]
dat.sub[is.na(p.num), p.num := 1]
dat.sub[, q := p.adjust(p.num, method = "fdr")]
dat.sub[, sig := q < 0.05]

sig.pairs <- unique(dat.sub[sig == TRUE, .(`Target gene name`, TRAIT)])
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
fwrite(out, "test/validated-6k-MR.tsv", sep = "\t", quote = FALSE)

cat("--- MR validation ---\n")
cat(sprintf("Denominator (N1): %d\n", N1))
cat(sprintf("Numerator   (N2): %d\n", N2))
cat(sprintf("Proportion N2/N1: %.4f\n",
            if (N1 > 0) N2 / N1 else NA_real_))
