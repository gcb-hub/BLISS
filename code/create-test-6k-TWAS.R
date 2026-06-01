#!/usr/bin/env Rscript
# =============================================================================
# create-test-6k-TWAS.R
#
# Parallel to create-test-6k-PWAS.R, driven by TWAS.RData. Gene column is
# dat$gene_name (already a symbol). Row-level QC: pred_perf_r2 >= 0.01.
#
# Design principle
# ----------------
# A (gene, disease) cell counts as a NEGATIVE only if TWAS actually tested it.
# Cells the method could not analyze are UNMEASURABLE and excluded -- they are
# NOT treated as negatives.
#
# Tested scope, per cell (g, d):
#   TWAS has at least one row in `dat` (post-QC: pred_perf_r2 >= 0.01) with
#   gene_name = g and TRAIT mapping via lookup-MVP to one of d's ontology IDs.
#
# LABEL:
#   1 if (g, d) appears in universe-6k, else 0.
#
# Outputs:
#   data/gene-6k-TWAS.txt
#   data/disease-6k-TWAS.txt
#   data/universe-6k-TWAS.tsv
#   test/test-6k-TWAS.tsv
#
# Run from the project root:
#   Rscript code/create-test-6k-TWAS.R
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
fwrite(test, "test/test-6k-TWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Helper outputs
# -----------------------------------------------------------------------------
genes.tested    <- sort(unique(test[["Target gene name"]]))
diseases.tested <- sort(unique(test[["Disease indication"]]))
writeLines(genes.tested,    "data/gene-6k-TWAS.txt")
writeLines(diseases.tested, "data/disease-6k-TWAS.txt")

pos.pairs <- test[LABEL == 1, .(`Target gene name`, `Disease indication`)]
univ.sub  <- merge(universe.dedup, pos.pairs,
                   by = c("Target gene name", "Disease indication"))
fwrite(univ.sub, "data/universe-6k-TWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes tested:     %d\n", length(genes.tested)))
cat(sprintf("diseases tested:  %d\n", length(diseases.tested)))
cat(sprintf("test cells:       %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
cat(sprintf("universe rows in scope: %d\n", nrow(univ.sub)))
