#!/usr/bin/env Rscript
# =============================================================================
# create-test-6k-MR.R
#
# Parallel to create-test-6k-PWAS.R, driven by MR.RData.
#
# Design principle
# ----------------
# A (gene, disease) cell counts as a NEGATIVE only if MR actually tested it.
# Cells the method could not analyze (no usable pQTL signal, or no gene-trait
# result reachable through the lookup) are UNMEASURABLE and excluded -- they
# are NOT treated as negatives.
#
# Tested scope, per cell (g, d):
#   MR has at least one row in `dat` whose `exception` is neither PQTL_EMPTY
#   nor PQTL_NO_SIGNAL (i.e. the relevance assumption is met), whose `protein`
#   maps via lookup-Olink to gene g, and whose TRAIT maps via lookup-MVP to
#   one of d's ontology IDs from universe-6k.
#
# LABEL:
#   1 if (g, d) appears in universe-6k, else 0.
#
# Outputs:
#   data/gene-6k-MR.txt
#   data/disease-6k-MR.txt
#   data/universe-6k-MR.tsv
#   test/test-6k-MR.tsv
#
# Run from the project root:
#   Rscript code/create-test-6k-MR.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/MR.RData")                # provides `dat`
setDT(dat)

universe <- fread("data/universe-6k.tsv",      sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv",  sep = "\t", quote = "")
olink    <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: QC-filtered MR results in (gene, TRAIT) form
# -----------------------------------------------------------------------------
# Keep only MR rows that produced a usable estimate. Categories kept:
#   blank             -- clean MR result
#   ONLY_ONE_IV       -- single-IV point estimate (still usable)
#   SUMSTAT_NO_SIGNAL -- weak outcome but a numeric p was emitted
# Excluded categories (no usable result -> unmeasurable, not negative):
#   PQTL_EMPTY, PQTL_NO_SIGNAL  -- relevance assumption violated
#   NO_IV, NO_OVERLAP           -- no instrument / no SNP overlap
dat.f <- dat[exception %in% c("blank", "ONLY_ONE_IV", "SUMSTAT_NO_SIGNAL")]

# Translate protein (OlinkID) -> gene symbol (HGNC.symbol).
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

dat.gt <- unique(merge(dat.f[, .(protein, TRAIT)], prot2gene,
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
fwrite(test, "test/test-6k-MR.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Helper outputs
# -----------------------------------------------------------------------------
genes.tested    <- sort(unique(test[["Target gene name"]]))
diseases.tested <- sort(unique(test[["Disease indication"]]))
writeLines(genes.tested,    "data/gene-6k-MR.txt")
writeLines(diseases.tested, "data/disease-6k-MR.txt")

pos.pairs <- test[LABEL == 1, .(`Target gene name`, `Disease indication`)]
univ.sub  <- merge(universe.dedup, pos.pairs,
                   by = c("Target gene name", "Disease indication"))
fwrite(univ.sub, "data/universe-6k-MR.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes tested:     %d\n", length(genes.tested)))
cat(sprintf("diseases tested:  %d\n", length(diseases.tested)))
cat(sprintf("test cells:       %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
cat(sprintf("universe rows in scope: %d\n", nrow(univ.sub)))
