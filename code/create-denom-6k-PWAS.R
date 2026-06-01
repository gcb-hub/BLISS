#!/usr/bin/env Rscript
# =============================================================================
# create-denom-6k-PWAS.R
#
# Build the PWAS "denominator" file: the subset of rows in
# data/universe-6k.tsv that PWAS can reasonably answer, i.e. cells where PWAS
# has a QC-passing result for the row's `Target gene name` with a TRAIT whose
# ontology ID overlaps the row's `Disease indication ID`. Raw row granularity
# is preserved (no dedup on gene x disease), so repeat drugs against the same
# (gene, disease) cell each contribute their own row.
#
# "Reasonably answered", per universe row (g, d, disease ids D):
#   PWAS has at least one row in `dat` (post-QC: r2.c >= 0.01) with
#   ID.inner = g and TRAIT whose mapped ontology IDs (via lookup-MVP's
#   traitFromSourceMappedIds) intersect D.
#
# Output:
#   test/denom-6k-PWAS.tsv -- all original columns of universe-6k.tsv plus a
#                              TRAIT column of `;`-joined candidate TRAITs
#                              (the TRAITs that could validate this row in the
#                              accompanying validate-6k-PWAS.R script).
#
# Run from the project root:
#   Rscript code/create-denom-6k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/PWAS.RData")              # provides `dat`
setDT(dat)

universe <- fread("data/universe-6k.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: QC-filtered PWAS results in (gene, TRAIT) form
# -----------------------------------------------------------------------------
dat.f  <- dat[r2.c >= 0.01]
dat.gt <- unique(dat.f[, .(`Target gene name` = ID.inner, TRAIT)])

# -----------------------------------------------------------------------------
# Step 2: TRAIT -> ontology ID (long)
# -----------------------------------------------------------------------------
lookup.long <- lookup[, .(
  id = unlist(strsplit(traitFromSourceMappedIds, ";"))
), by = TRAIT]
lookup.long[, id := trimws(id)]
lookup.long <- lookup.long[nzchar(id)]

# Method's (gene, TRAIT, id) triples: every way a PWAS (gene, TRAIT) result
# can reach some ontology ID.
gene.trait.id <- merge(dat.gt, lookup.long,
                       by = "TRAIT", allow.cartesian = TRUE)

# -----------------------------------------------------------------------------
# Step 3: Universe -> per-row (gene, id) long form, preserving row_id
# -----------------------------------------------------------------------------
universe[, row_id := .I]
universe.long <- universe[, .(
  id = trimws(unlist(strsplit(`Disease indication ID`, ",")))
), by = .(row_id, `Target gene name`, `Disease indication`)]
universe.long <- universe.long[nzchar(id)]

# -----------------------------------------------------------------------------
# Step 4: Join -- rows that PWAS can answer, with their candidate TRAITs
# -----------------------------------------------------------------------------
hits <- merge(universe.long, gene.trait.id,
              by = c("Target gene name", "id"),
              allow.cartesian = TRUE)

row.cand <- hits[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = row_id]

# -----------------------------------------------------------------------------
# Step 5: Materialise the denominator file at raw universe-row granularity
# -----------------------------------------------------------------------------
denom <- merge(universe, row.cand, by = "row_id")
setorder(denom, row_id)
denom[, row_id := NULL]

fwrite(denom, "test/denom-6k-PWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("Universe total:  %d\n", nrow(universe)))
cat(sprintf("Denominator:     %d\n", nrow(denom)))
