#!/usr/bin/env Rscript
# =============================================================================
# create-test-13k-coloc.R
#
# Parallel to create-test-13k-MR.R, but driven by coloc.RData. Rows whose
# `exception` is NO_OVERLAP, PQTL_EMPTY, or PQTL_NO_SIGNAL are dropped before
# anything else -- those cells were not actually analyzed.
#
# As with the 13k PWAS / MR pipelines:
#   - universe-13k.tsv carries one MeSH ID per row in `indication_mesh_id`,
#     so TRAIT -> disease is a direct lookup via lookup-MVP.tsv's `MESH`.
#   - LABEL in the test file comes from succ_3_a (TRUE -> 1, else 0).
#   - dat$protein (OlinkID) is translated to gene symbol via lookup-Olink.tsv.
#
# Outputs:
#   data/gene-13k-coloc.txt
#   data/disease-13k-coloc.txt
#   data/universe-13k-coloc.tsv
#   test/test-13k-coloc.tsv     (Target gene name | MESH | TRAIT | LABEL)
#
# Run from the project root:
#   Rscript code/create-test-13k-coloc.R
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

universe <- fread("data/universe-13k-old.tsv",  sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv",  sep = "\t", quote = "")
olink    <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")

# universe filter disabled: using full Minikel universe

# -----------------------------------------------------------------------------
# Part 1: gene list
# -----------------------------------------------------------------------------
# Translate protein (OlinkID) -> gene symbol via the Olink lookup.
prot2gene <- unique(olink[, .(protein = OlinkID, gene = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(gene)]

genes.init <- unique(merge(
  data.table(protein = unique(dat$protein)),
  prot2gene, by = "protein"
)$gene)

genes.univ <- unique(universe$gene)
genes.keep <- intersect(genes.init, genes.univ)

writeLines(genes.keep, "data/gene-13k-coloc.txt")

# -----------------------------------------------------------------------------
# Part 2: disease list
# -----------------------------------------------------------------------------
# Translate dat$TRAIT into MeSH IDs via lookup's `MESH` column (blanks dropped).
trait2mesh <- lookup[nzchar(MESH), .(TRAIT, MESH)]

mesh.init <- unique(trait2mesh[TRAIT %in% unique(dat$TRAIT)]$MESH)
mesh.univ <- unique(universe$indication_mesh_id)
diseases.keep <- intersect(mesh.init, mesh.univ)

writeLines(diseases.keep, "data/disease-13k-coloc.txt")

# -----------------------------------------------------------------------------
# Part 3: shrink + dedup universe
# -----------------------------------------------------------------------------
# Dedup on (gene, indication_mesh_id) preferring succ_3_a == TRUE so positive
# labels survive ties.
univ.sub <- universe[gene %in% genes.keep &
                     indication_mesh_id %in% diseases.keep]

setorder(univ.sub, gene, indication_mesh_id, -succ_3_a)
univ.sub <- unique(univ.sub, by = c("gene", "indication_mesh_id"))

fwrite(univ.sub, "data/universe-13k-coloc.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Part 4: test file
# -----------------------------------------------------------------------------
# Scaffold: one row per (gene, MESH); LABEL = 1 iff succ_3_a == TRUE, else 0.
test <- univ.sub[, .(
  `Target gene name` = gene,
  MESH               = indication_mesh_id,
  LABEL              = as.integer(succ_3_a)
)]
test[is.na(LABEL), LABEL := 0L]

# Aggregate supporting coloc TRAIT codes per (gene, MESH).
dat.long <- merge(
  dat[, .(protein, TRAIT)],
  prot2gene, by = "protein"
)[, .(`Target gene name` = gene, TRAIT)]

dat.long <- merge(dat.long, trait2mesh,
                  by = "TRAIT", allow.cartesian = TRUE)

trait.agg <- dat.long[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, MESH)]

test <- merge(test, trait.agg,
              by = c("Target gene name", "MESH"),
              all.x = TRUE)
test[is.na(TRAIT), TRAIT := ""]

setcolorder(test, c("Target gene name", "MESH", "TRAIT", "LABEL"))

fwrite(test, "test/test-13k-coloc.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes kept:       %d\n", length(genes.keep)))
cat(sprintf("diseases kept:    %d\n", length(diseases.keep)))
cat(sprintf("universe rows:    %d\n", nrow(univ.sub)))
cat(sprintf("test rows:        %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
