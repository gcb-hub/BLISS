#!/usr/bin/env Rscript
# =============================================================================
# create-test-13k-PWAS.R
#
# Build a benchmarking test set that pairs PWAS gene-trait associations with
# the 13k drug-indication universe (universe-13k.tsv).
#
# Unlike the 6k pipeline (which used semicolon-separated EFO/MONDO IDs), the
# 13k universe carries a single MeSH ID per row in `indication_mesh_id`, so
# the TRAIT -> disease mapping is a direct lookup via the `MESH` column of
# dictionary/lookup-MVP.tsv.
#
# Outputs:
#   data/gene-13k-PWAS.txt      -- PWAS-supported genes that appear in universe
#   data/disease-13k-PWAS.txt   -- MeSH IDs shared by PWAS (via lookup) & universe
#   data/universe-13k-PWAS.tsv  -- universe shrunken + deduped to those
#   data/test-13k-PWAS.tsv      -- one row per (gene, MESH) with LABEL = succ_3_a
#                                  (Target gene name | MESH | TRAIT | LABEL)
#
# Run from the project root:
#   Rscript code/create-test-13k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/PWAS.RData")              # provides `dat`
setDT(dat)

universe <- fread("data/universe-13k-old.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# universe filter disabled: using full Minikel universe

# -----------------------------------------------------------------------------
# Part 1: gene list
# -----------------------------------------------------------------------------
# Filter PWAS by heritability threshold r2.c >= 0.01, take unique ID.inner as
# the initial gene set, then intersect with universe's `gene` column.
dat.f <- dat[r2.c >= 0.01]

genes.init <- unique(dat.f$ID.inner)
genes.univ <- unique(universe$gene)
genes.keep <- intersect(genes.init, genes.univ)

writeLines(genes.keep, "data/gene-13k-PWAS.txt")

# -----------------------------------------------------------------------------
# Part 2: disease list
# -----------------------------------------------------------------------------
# Translate dat$TRAIT into MeSH IDs via lookup's `MESH` column (only ~2,810 of
# 4,173 lookup rows have a MeSH mapping; blanks are dropped).
trait2mesh <- lookup[nzchar(MESH), .(TRAIT, MESH)]

mesh.init <- unique(trait2mesh[TRAIT %in% unique(dat$TRAIT)]$MESH)
mesh.univ <- unique(universe$indication_mesh_id)
diseases.keep <- intersect(mesh.init, mesh.univ)

writeLines(diseases.keep, "data/disease-13k-PWAS.txt")

# -----------------------------------------------------------------------------
# Part 3: shrink + dedup universe
# -----------------------------------------------------------------------------
# Keep rows where gene is in the gene list AND indication_mesh_id is in the
# disease list. Dedup on (gene, indication_mesh_id), preferring rows with
# succ_3_a == TRUE so positive labels survive ties.
univ.sub <- universe[gene %in% genes.keep &
                     indication_mesh_id %in% diseases.keep]

setorder(univ.sub, gene, indication_mesh_id, -succ_3_a)
univ.sub <- unique(univ.sub, by = c("gene", "indication_mesh_id"))

fwrite(univ.sub, "data/universe-13k-PWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Part 4: test file
# -----------------------------------------------------------------------------
# Scaffold: one row per (gene, MESH) from the deduped shrunken universe.
# LABEL = 1 iff succ_3_a == TRUE for that row, else 0.
test <- univ.sub[, .(
  `Target gene name` = gene,
  MESH               = indication_mesh_id,
  LABEL              = as.integer(succ_3_a)
)]
test[is.na(LABEL), LABEL := 0L]

# Aggregate supporting PWAS TRAIT codes per (gene, MESH): join filtered PWAS
# with the TRAIT->MESH lookup, then collapse.
dat.long <- merge(
  dat.f[, .(`Target gene name` = ID.inner, TRAIT)],
  trait2mesh,
  by = "TRAIT", allow.cartesian = TRUE
)

trait.agg <- dat.long[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, MESH)]

test <- merge(test, trait.agg,
              by = c("Target gene name", "MESH"),
              all.x = TRUE)
test[is.na(TRAIT), TRAIT := ""]

setcolorder(test, c("Target gene name", "MESH", "TRAIT", "LABEL"))

fwrite(test, "test/test-13k-PWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes kept:       %d\n", length(genes.keep)))
cat(sprintf("diseases kept:    %d\n", length(diseases.keep)))
cat(sprintf("universe rows:    %d\n", nrow(univ.sub)))
cat(sprintf("test rows:        %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
