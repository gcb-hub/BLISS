#!/usr/bin/env Rscript
# =============================================================================
# create-test-6k-coloc.R
#
# Parallel to create-test-6k-PWAS.R, driven by coloc.RData. Rows whose
# `exception` is NO_OVERLAP, PQTL_EMPTY, or PQTL_NO_SIGNAL are dropped first --
# those cells were not actually analyzed and should not count as tested.
#
# Design principle
# ----------------
# A (gene, disease) cell counts as a NEGATIVE only if coloc actually tested
# it. Cells the method could not analyze are UNMEASURABLE and excluded -- they
# are NOT treated as negatives.
#
# Tested scope, per cell (g, d):
#   coloc has at least one row in `dat` whose `protein` maps via lookup-Olink
#   to gene g and whose TRAIT maps via lookup-MVP to one of d's ontology IDs.
#
# LABEL:
#   1 if (g, d) appears in universe-6k, else 0.
#
# Outputs:
#   data/gene-6k-coloc.txt
#   data/disease-6k-coloc.txt
#   data/universe-6k-coloc.tsv
#   test/test-6k-coloc.tsv
#
# Run from the project root:
#   Rscript code/create-test-6k-coloc.R
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
fwrite(test, "test/test-6k-coloc.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Helper outputs
# -----------------------------------------------------------------------------
genes.tested    <- sort(unique(test[["Target gene name"]]))
diseases.tested <- sort(unique(test[["Disease indication"]]))
writeLines(genes.tested,    "data/gene-6k-coloc.txt")
writeLines(diseases.tested, "data/disease-6k-coloc.txt")

pos.pairs <- test[LABEL == 1, .(`Target gene name`, `Disease indication`)]
univ.sub  <- merge(universe.dedup, pos.pairs,
                   by = c("Target gene name", "Disease indication"))
fwrite(univ.sub, "data/universe-6k-coloc.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes tested:     %d\n", length(genes.tested)))
cat(sprintf("diseases tested:  %d\n", length(diseases.tested)))
cat(sprintf("test cells:       %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
cat(sprintf("universe rows in scope: %d\n", nrow(univ.sub)))
