#!/usr/bin/env Rscript
# =============================================================================
# create-denom-6k-closest.R
#
# Build the closest "denominator" file: the subset of rows in
# data/universe-6k.tsv that closest can reasonably answer, at raw
# universe-row granularity (no gene x disease dedup).
#
# "Reasonably answered", per universe row (g, d, disease ids D):
#   MVP has at least one row with `Lead Variant Gene (Nearest)` = g and
#   `Trait` = t such that t maps via lookup-MVP$traitFromSourceMappedIds
#   to an ontology id in D.
#
# Output:
#   test/denom-6k-closest.tsv -- original universe-6k.tsv columns plus a
#                                 TRAIT column listing the `;`-joined
#                                 candidate TRAITs that could validate this
#                                 row.
#
# Run from the project root:
#   Rscript code/create-denom-6k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
loci <- fread("data/loci-from-MVP-corrected.csv")
loci <- loci[nzchar(`Lead Variant Gene (Nearest)`) &
             !is.na(`Lead Variant Gene (Nearest)`)]
loci <- loci[nzchar(Trait) & !is.na(Trait)]
loci[, p.num := suppressWarnings(as.numeric(`P-Value`))]
loci <- loci[!is.na(p.num) & p.num < 5e-8]

universe <- fread("data/universe-6k.tsv",     sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: MVP loci in (gene, TRAIT) form.
# -----------------------------------------------------------------------------
dat.gt <- unique(loci[, .(
  `Target gene name` = `Lead Variant Gene (Nearest)`,
  TRAIT              = Trait
)])

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
# Step 4: Join -- rows that closest can answer, with their candidate TRAITs
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

fwrite(denom, "test/denom-6k-closest.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("Universe total:  %d\n", nrow(universe)))
cat(sprintf("Denominator:     %d\n", nrow(denom)))
