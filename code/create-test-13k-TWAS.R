#!/usr/bin/env Rscript
# =============================================================================
# create-test-13k-TWAS.R
#
# Parallel to create-test-13k-PWAS.R, but driven by TWAS.RData. Gene column
# is dat$gene_name (already a symbol). Filter rows by pred_perf_r2 >= 0.01
# before taking the initial gene list.
#
# As with the other 13k pipelines:
#   - universe-13k.tsv carries one MeSH ID per row in `indication_mesh_id`,
#     so TRAIT -> disease is a direct lookup via lookup-MVP.tsv's `MESH`.
#   - LABEL in the test file comes from succ_3_a (TRUE -> 1, else 0).
#
# Outputs:
#   data/gene-13k-TWAS.txt
#   data/disease-13k-TWAS.txt
#   data/universe-13k-TWAS.tsv
#   test/test-13k-TWAS.tsv     (Target gene name | MESH | TRAIT | LABEL)
#
# Run from the project root:
#   Rscript code/create-test-13k-TWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/TWAS.RData")               # provides `dat`
setDT(dat)

universe <- fread("data/universe-13k-old.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# universe filter disabled: using full Minikel universe

# -----------------------------------------------------------------------------
# Part 1: gene list
# -----------------------------------------------------------------------------
# Filter TWAS by prediction-model R^2 >= 0.01; gene_name is already a symbol.
dat.f <- dat[pred_perf_r2 >= 0.01]

genes.init <- unique(dat.f$gene_name)
genes.univ <- unique(universe$gene)
genes.keep <- intersect(genes.init, genes.univ)

writeLines(genes.keep, "data/gene-13k-TWAS.txt")

# -----------------------------------------------------------------------------
# Part 2: disease list
# -----------------------------------------------------------------------------
# Translate dat$TRAIT into MeSH IDs via lookup's `MESH` column (blanks dropped).
trait2mesh <- lookup[nzchar(MESH), .(TRAIT, MESH)]

mesh.init <- unique(trait2mesh[TRAIT %in% unique(dat$TRAIT)]$MESH)
mesh.univ <- unique(universe$indication_mesh_id)
diseases.keep <- intersect(mesh.init, mesh.univ)

writeLines(diseases.keep, "data/disease-13k-TWAS.txt")

# -----------------------------------------------------------------------------
# Part 3: shrink + dedup universe
# -----------------------------------------------------------------------------
# Dedup on (gene, indication_mesh_id) preferring succ_3_a == TRUE so positive
# labels survive ties.
univ.sub <- universe[gene %in% genes.keep &
                     indication_mesh_id %in% diseases.keep]

setorder(univ.sub, gene, indication_mesh_id, -succ_3_a)
univ.sub <- unique(univ.sub, by = c("gene", "indication_mesh_id"))

fwrite(univ.sub, "data/universe-13k-TWAS.tsv", sep = "\t", quote = FALSE)

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

# Aggregate supporting TWAS TRAIT codes per (gene, MESH).
dat.long <- merge(
  dat.f[, .(`Target gene name` = gene_name, TRAIT)],
  trait2mesh, by = "TRAIT", allow.cartesian = TRUE
)

trait.agg <- dat.long[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, MESH)]

test <- merge(test, trait.agg,
              by = c("Target gene name", "MESH"),
              all.x = TRUE)
test[is.na(TRAIT), TRAIT := ""]

setcolorder(test, c("Target gene name", "MESH", "TRAIT", "LABEL"))

fwrite(test, "test/test-13k-TWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes kept:       %d\n", length(genes.keep)))
cat(sprintf("diseases kept:    %d\n", length(diseases.keep)))
cat(sprintf("universe rows:    %d\n", nrow(univ.sub)))
cat(sprintf("test rows:        %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
