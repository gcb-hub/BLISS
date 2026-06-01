#!/usr/bin/env Rscript
# =============================================================================
# create-test-13k-closest.R
#
# Parallel to create-test-13k-coloc.R, but driven by
# data/loci-from-MVP-corrected.csv (EUR rows only).
# The "closest" method's prediction for each MVP lead variant is the
# `Lead Variant Gene (Nearest)` column. Rows with an empty nearest-gene are
# dropped before anything else -- those loci have no prediction.
#
# As with the 13k coloc / MR / PWAS / TWAS pipelines:
#   - universe-13k carries one MeSH ID per row in `indication_mesh_id`,
#     so TRAIT -> disease is a direct lookup via lookup-MVP.tsv's `MESH`.
#   - LABEL in the test file comes from succ_3_a (TRUE -> 1, else 0).
#
# Per user direction, the truth source is data/universe-13k-old.tsv (13k rows),
# not data/universe-13k.tsv.
#
# Outputs:
#   data/gene-13k-closest.txt
#   data/disease-13k-closest.txt
#   data/universe-13k-closest.tsv
#   test/test-13k-closest.tsv     (Target gene name | MESH | TRAIT | LABEL)
#
# Run from the project root:
#   Rscript code/create-test-13k-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
loci <- fread("data/loci-from-MVP-corrected.csv")
if (colnames(loci)[1] != "Population") setnames(loci, 1, "Population")
loci <- loci[Population == "EUR"]
loci <- loci[nzchar(`Lead Variant Gene (Nearest)`) &
             !is.na(`Lead Variant Gene (Nearest)`)]
loci <- loci[nzchar(Trait) & !is.na(Trait)]
loci[, p.num := suppressWarnings(as.numeric(`P-Value`))]
loci <- loci[!is.na(p.num) & p.num < 5e-8]

universe <- fread("data/universe-13k-old.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# universe filter disabled: using full Minikel universe

# Synthesize the equivalent of `dat`: one row per (gene, TRAIT) prediction.
dat <- unique(loci[, .(
  `Target gene name` = `Lead Variant Gene (Nearest)`,
  TRAIT              = Trait
)])

# -----------------------------------------------------------------------------
# Part 1: gene list
# -----------------------------------------------------------------------------
genes.init <- unique(dat$`Target gene name`)
genes.univ <- unique(universe$gene)
genes.keep <- intersect(genes.init, genes.univ)

writeLines(genes.keep, "data/gene-13k-closest.txt")

# -----------------------------------------------------------------------------
# Part 2: disease list
# -----------------------------------------------------------------------------
trait2mesh <- lookup[nzchar(MESH), .(TRAIT, MESH)]

mesh.init <- unique(trait2mesh[TRAIT %in% unique(dat$TRAIT)]$MESH)
mesh.univ <- unique(universe$indication_mesh_id)
diseases.keep <- intersect(mesh.init, mesh.univ)

writeLines(diseases.keep, "data/disease-13k-closest.txt")

# -----------------------------------------------------------------------------
# Part 3: shrink + dedup universe
# -----------------------------------------------------------------------------
univ.sub <- universe[gene %in% genes.keep &
                     indication_mesh_id %in% diseases.keep]

setorder(univ.sub, gene, indication_mesh_id, -succ_3_a)
univ.sub <- unique(univ.sub, by = c("gene", "indication_mesh_id"))

fwrite(univ.sub, "data/universe-13k-closest.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Part 4: test file
# -----------------------------------------------------------------------------
test <- univ.sub[, .(
  `Target gene name` = gene,
  MESH               = indication_mesh_id,
  LABEL              = as.integer(succ_3_a)
)]
test[is.na(LABEL), LABEL := 0L]

dat.long <- merge(dat, trait2mesh,
                  by = "TRAIT", allow.cartesian = TRUE)

trait.agg <- dat.long[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, MESH)]

test <- merge(test, trait.agg,
              by = c("Target gene name", "MESH"),
              all.x = TRUE)
test[is.na(TRAIT), TRAIT := ""]

setcolorder(test, c("Target gene name", "MESH", "TRAIT", "LABEL"))

fwrite(test, "test/test-13k-closest.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes kept:       %d\n", length(genes.keep)))
cat(sprintf("diseases kept:    %d\n", length(diseases.keep)))
cat(sprintf("universe rows:    %d\n", nrow(univ.sub)))
cat(sprintf("test rows:        %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
