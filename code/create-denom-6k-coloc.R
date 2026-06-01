#!/usr/bin/env Rscript
# =============================================================================
# create-denom-6k-coloc.R
#
# Build the coloc "denominator" file: the subset of rows in
# data/universe-6k.tsv that coloc can reasonably answer, at raw universe-row
# granularity (no gene x disease dedup).
#
# "Reasonably answered", per universe row (g, d, disease ids D):
#   coloc has at least one row in `dat` (post-QC: exception not in
#   NO_OVERLAP / PQTL_EMPTY / PQTL_NO_SIGNAL) whose `protein` maps via
#   lookup-Olink to gene g and whose TRAIT maps via lookup-MVP to an ontology
#   ID that intersects D.
#
# Output:
#   test/denom-6k-coloc.tsv -- original universe-6k.tsv columns plus a TRAIT
#                               column listing the `;`-joined candidate TRAITs
#                               that could validate this row.
#
# Run from the project root:
#   Rscript code/create-denom-6k-coloc.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/coloc.RData")              # provides `dat`
setDT(dat)
dat <- dat[!exception %in% c("NO_OVERLAP", "PQTL_EMPTY", "PQTL_NO_SIGNAL")]

universe <- fread("data/universe-6k.tsv",      sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv",  sep = "\t", quote = "")
olink    <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: coloc results in (gene, TRAIT) form
# -----------------------------------------------------------------------------
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

dat.gt <- unique(merge(dat[, .(protein, TRAIT)], prot2gene,
                       by = "protein")[, .(`Target gene name`, TRAIT)])

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
# Step 4: Join -- rows that coloc can answer, with their candidate TRAITs
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

fwrite(denom, "test/denom-6k-coloc.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("Universe total:  %d\n", nrow(universe)))
cat(sprintf("Denominator:     %d\n", nrow(denom)))
