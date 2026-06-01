#!/usr/bin/env Rscript
# =============================================================================
# validate-6k-cohorts.R
#
# Run validate-6k-PWAS.R-style analysis for additional cohorts:
#   deCODE       -- SomaScan EUR, probe-level h2 filter (h2-SomaScan-EUR.csv)
#   ARIC         -- SomaScan EUR, gene-level  h2 filter (h2-SomaScan-EUR.csv)
#   ARIC-AA      -- SomaScan AFR, gene-level  h2 filter (h2-SomaScan-AFR.csv)
#   UKB-AFR+EUR  -- Olink, r2.c >= 0.01 (mirrors PWAS/UKB-PPP)
#
# For each cohort:
#   * Build a cohort-specific denominator from data/universe-6k.tsv: a row is
#     "answerable" if the cohort has a QC-passing (gene, TRAIT) result whose
#     mapped ontology IDs overlap the row's Disease indication IDs.
#   * Validate: FDR within the denominator scope on p.a.c, q < 0.05; row is
#     validated if any of its candidate (gene, TRAIT) pairs is significant.
#   * Report N2 (validated) / N1 (denominator) / N2/N1.
#
# Run from the project root:
#   Rscript code/validate-6k-cohorts.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

universe <- fread("data/universe-6k.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

lookup.long <- lookup[, .(
  id = unlist(strsplit(traitFromSourceMappedIds, ";"))
), by = TRAIT]
lookup.long[, id := trimws(id)]
lookup.long <- lookup.long[nzchar(id)]

universe[, row_id := .I]
universe.long <- universe[, .(
  id = trimws(unlist(strsplit(`Disease indication ID`, ",")))
), by = .(row_id, `Target gene name`, `Disease indication`)]
universe.long <- universe.long[nzchar(id)]

# -----------------------------------------------------------------------------
# Per-cohort loader: returns a data.table with columns (gene, TRAIT, p.a.c)
# after applying the cohort-appropriate QC filter.
# -----------------------------------------------------------------------------
load.cohort <- function(name) {
  if (name == "deCODE") {
    load("result-compiled/deCODE.RData"); setDT(dat)
    h2 <- fread("data/h2-SomaScan-EUR.csv")
    h2.seqids <- h2[`cis-h2` > 0.01, `SOMAmer ID`]
    id.parts <- tstrsplit(dat$ID.inner, "_", fixed = TRUE)
    dat[, SeqId := paste0("SeqId_", id.parts[[4]], "_", id.parts[[5]])]
    dat[, gene  := id.parts[[6]]]
    dat <- dat[SeqId %in% h2.seqids]
    return(dat[, .(gene, TRAIT, p.a.c)])
  }
  if (name == "ARIC") {
    load("result-compiled/ARIC.RData"); setDT(dat)
    h2 <- fread("data/h2-SomaScan-EUR.csv")
    h2.genes <- unique(h2[`cis-h2` > 0.01, `Entrez Gene Symbol`])
    h2.genes <- h2.genes[nzchar(h2.genes)]
    dat <- dat[ID.inner %in% h2.genes]
    return(dat[, .(gene = ID.inner, TRAIT, p.a.c)])
  }
  if (name == "ARIC-AA") {
    load("result-compiled/ARIC-AA.RData"); setDT(dat)
    h2 <- fread("data/h2-SomaScan-AFR.csv")
    h2.genes <- unique(h2[`cis-h2` > 0.01, `Entrez Gene Symbol`])
    h2.genes <- h2.genes[nzchar(h2.genes)]
    dat <- dat[ID.inner %in% h2.genes]
    return(dat[, .(gene = ID.inner, TRAIT, p.a.c)])
  }
  if (name == "UKB-AFR+EUR") {
    load("result-compiled/UKB-AFR+EUR.RData"); setDT(dat)
    dat <- dat[r2.c >= 0.01]
    return(dat[, .(gene = ID.inner, TRAIT, p.a.c)])
  }
  stop("unknown cohort: ", name)
}

# -----------------------------------------------------------------------------
# Validation pipeline given a (gene, TRAIT, p.a.c) data.table for a cohort.
# Returns named list with N1, N2, ratio.
# -----------------------------------------------------------------------------
validate.one <- function(dat) {
  dat.gt <- unique(dat[, .(`Target gene name` = gene, TRAIT)])

  gene.trait.id <- merge(dat.gt, lookup.long,
                         by = "TRAIT", allow.cartesian = TRUE)

  hits <- merge(universe.long, gene.trait.id,
                by = c("Target gene name", "id"),
                allow.cartesian = TRUE)

  row.cand <- hits[, .(
    TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
  ), by = row_id]

  denom <- merge(universe, row.cand, by = "row_id")

  denom[, denom_row := .I]
  denom.long <- denom[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(denom_row, `Target gene name`)]
  denom.long[, TRAIT := trimws(TRAIT)]
  denom.long <- denom.long[nzchar(TRAIT)]

  scope.pairs <- unique(denom.long[, .(gene = `Target gene name`, TRAIT)])
  dat.sub <- merge(dat, scope.pairs, by = c("gene", "TRAIT"))

  dat.sub[, q := p.adjust(p.a.c, method = "fdr")]
  dat.sub[, sig := !is.na(q) & q < 0.05]

  sig.pairs <- unique(dat.sub[sig == TRUE,
                              .(`Target gene name` = gene, TRAIT, sig = TRUE)])

  denom.scored <- merge(denom.long, sig.pairs,
                        by = c("Target gene name", "TRAIT"),
                        all.x = TRUE)
  denom.scored[is.na(sig), sig := FALSE]

  validated_rows <- unique(denom.scored[sig == TRUE, denom_row])

  N1 <- nrow(denom)
  N2 <- length(validated_rows)
  list(N1 = N1, N2 = N2)
}

cohorts <- c("deCODE", "ARIC", "ARIC-AA", "UKB-AFR+EUR")
results <- list()
for (c in cohorts) {
  cat("=== ", c, " ===\n", sep = "")
  dat.c <- load.cohort(c)
  r <- validate.one(dat.c)
  results[[c]] <- r
  cat(sprintf("Answerable (N1): %d\n", r$N1))
  cat(sprintf("Validated  (N2): %d\n", r$N2))
  cat(sprintf("Recovery rate:   %.1f%%\n\n",
              100 * (if (r$N1 > 0) r$N2 / r$N1 else NA_real_)))
}

cat("=== summary ===\n")
cat("Method\tValidated trios recovered\tAnswerable trios\tRecovery rate\n")
for (c in cohorts) {
  r <- results[[c]]
  rate <- if (r$N1 > 0) sprintf("%.1f%%", 100 * r$N2 / r$N1) else "NA"
  cat(sprintf("%s\t%d\t%d\t%s\n", c, r$N2, r$N1, rate))
}
