#!/usr/bin/env Rscript
# =============================================================================
# create-denom-6k-TWAS.R
#
# Build the TWAS "denominator" file: the subset of rows in
# data/universe-6k.tsv that TWAS can reasonably answer, at raw universe-row
# granularity (no gene x disease dedup).
#
# "Reasonably answered", per universe row (g, d, disease ids D):
#   TWAS has at least one row in `dat` (post-QC: pred_perf_r2 >= 0.01) with
#   gene_name = g and TRAIT whose mapped ontology IDs (via lookup-MVP's
#   traitFromSourceMappedIds) intersect D.
#
# Output:
#   test/denom-6k-TWAS.tsv -- original universe-6k.tsv columns plus a TRAIT
#                              column listing the `;`-joined candidate TRAITs
#                              that could validate this row.
#
# Run from the project root:
#   Rscript code/create-denom-6k-TWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/TWAS.RData")               # provides `dat`
setDT(dat)

universe <- fread("data/universe-6k.tsv",     sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: QC-filtered TWAS results in (gene, TRAIT) form
# -----------------------------------------------------------------------------
dat.f  <- dat[pred_perf_r2 >= 0.01]
dat.gt <- unique(dat.f[, .(`Target gene name` = gene_name, TRAIT)])

# -----------------------------------------------------------------------------
# Step 2: TRAIT -> ontology ID (long)
# -----------------------------------------------------------------------------
lookup.long <- lookup[, .(
  id = unlist(strsplit(traitFromSourceMappedIds, ";"))
), by = TRAIT]
lookup.long[, id := trimws(id)]
lookup.long <- lookup.long[nzchar(id)]

gene.trait.id <- merge(dat.gt, lookup.long,
                       by = "TRAIT", allow.cartesian = TRUE)

# -----------------------------------------------------------------------------
# Step 3: Universe -> per-row (gene, id) long, preserving row_id
# -----------------------------------------------------------------------------
universe[, row_id := .I]
universe.long <- universe[, .(
  id = trimws(unlist(strsplit(`Disease indication ID`, ",")))
), by = .(row_id, `Target gene name`, `Disease indication`)]
universe.long <- universe.long[nzchar(id)]

# -----------------------------------------------------------------------------
# Step 4: Join -- rows that TWAS can answer, with their candidate TRAITs
# -----------------------------------------------------------------------------
hits <- merge(universe.long, gene.trait.id,
              by = c("Target gene name", "id"),
              allow.cartesian = TRUE)

row.cand <- hits[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = row_id]

# -----------------------------------------------------------------------------
# Step 5: Materialise the denominator at raw universe-row granularity
# -----------------------------------------------------------------------------
denom <- merge(universe, row.cand, by = "row_id")
setorder(denom, row_id)
denom[, row_id := NULL]

fwrite(denom, "test/denom-6k-TWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("Universe total:  %d\n", nrow(universe)))
cat(sprintf("Denominator:     %d\n", nrow(denom)))
