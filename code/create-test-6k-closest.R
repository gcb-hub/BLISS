#!/usr/bin/env Rscript
# =============================================================================
# create-test-6k-closest.R
#
# Parallel to create-test-6k-coloc.R, but driven by data/loci-from-MVP.csv.
# The "closest" method's prediction for each MVP lead variant is the
# `Lead Variant Gene (Nearest)` column. There is no row-level QC beyond
# dropping rows where that column (or `Trait`) is empty/NA -- those loci
# have no usable prediction.
#
# Design principle (same as coloc / PWAS / MR / TWAS)
# ---------------------------------------------------
# A (gene, disease) cell counts as a NEGATIVE only if closest actually
# tested it. Cells the method cannot reach are UNMEASURABLE and excluded.
#
# Tested scope, per cell (g, d):
#   MVP has at least one row with `Lead Variant Gene (Nearest)` = g and
#   `Trait` = t such that t maps via lookup-MVP$traitFromSourceMappedIds
#   to one of d's `Disease indication ID` ontology IDs.
#
# LABEL:
#   1 if (g, d) appears in universe-6k, else 0.
#
# Outputs:
#   data/gene-6k-closest.txt
#   data/disease-6k-closest.txt
#   data/universe-6k-closest.tsv
#   test/test-6k-closest.tsv
#
# Run from the project root:
#   Rscript code/create-test-6k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
loci <- fread("data/loci-from-MVP.csv")
if (colnames(loci)[1] != "Trait") setnames(loci, 1, "Trait")
loci <- loci[nzchar(`Lead Variant Gene (Nearest)`) &
             !is.na(`Lead Variant Gene (Nearest)`)]
loci <- loci[nzchar(Trait) & !is.na(Trait)]
loci[, p.num := suppressWarnings(as.numeric(`P-value`))]
loci <- loci[!is.na(p.num) & p.num < 5e-8]

universe <- fread("data/universe-6k.tsv",     sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: MVP loci in (gene, TRAIT) form -- closest's analogue of `dat.gt`.
# No protein->gene translation; the gene is already an HGNC symbol.
# -----------------------------------------------------------------------------
dat.gt <- unique(loci[, .(
  `Target gene name` = `Lead Variant Gene (Nearest)`,
  TRAIT              = Trait
)])

# -----------------------------------------------------------------------------
# Step 2: TRAIT -> ontology ID (long form)
# -----------------------------------------------------------------------------
lookup.long <- lookup[, .(
  id = unlist(strsplit(traitFromSourceMappedIds, ";"))
), by = TRAIT]
lookup.long[, id := trimws(id)]
lookup.long <- lookup.long[nzchar(id)]

# -----------------------------------------------------------------------------
# Step 3: Universe disease (name) -> ontology IDs (long form)
# -----------------------------------------------------------------------------
disease.id.long <- universe[, .(
  id = trimws(unlist(strsplit(`Disease indication ID`, ",")))
), by = `Disease indication`]
disease.id.long <- unique(disease.id.long[nzchar(id)])

# -----------------------------------------------------------------------------
# Step 4: Tested (gene, disease) cells -- the natural background
# -----------------------------------------------------------------------------
gene.trait.id <- merge(dat.gt, lookup.long,
                       by = "TRAIT", allow.cartesian = TRUE)
tested.triples <- merge(gene.trait.id, disease.id.long,
                        by = "id", allow.cartesian = TRUE)

test <- tested.triples[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, `Disease indication`)]

# -----------------------------------------------------------------------------
# Step 5: LABEL against the full universe-6k
# -----------------------------------------------------------------------------
universe.dedup <- unique(universe, by = c("Target gene name", "Disease indication"))
pos.key <- universe.dedup[, paste(`Target gene name`, `Disease indication`,
                                  sep = "||")]
test[, LABEL := as.integer(
  paste(`Target gene name`, `Disease indication`, sep = "||") %in% pos.key
)]

setcolorder(test, c("Target gene name", "Disease indication", "TRAIT", "LABEL"))
fwrite(test, "test/test-6k-closest.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Helper outputs
# -----------------------------------------------------------------------------
genes.tested    <- sort(unique(test[["Target gene name"]]))
diseases.tested <- sort(unique(test[["Disease indication"]]))
writeLines(genes.tested,    "data/gene-6k-closest.txt")
writeLines(diseases.tested, "data/disease-6k-closest.txt")

pos.pairs <- test[LABEL == 1, .(`Target gene name`, `Disease indication`)]
univ.sub  <- merge(universe.dedup, pos.pairs,
                   by = c("Target gene name", "Disease indication"))
fwrite(univ.sub, "data/universe-6k-closest.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes tested:     %d\n", length(genes.tested)))
cat(sprintf("diseases tested:  %d\n", length(diseases.tested)))
cat(sprintf("test cells:       %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
cat(sprintf("universe rows in scope: %d\n", nrow(univ.sub)))
