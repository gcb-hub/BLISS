#!/usr/bin/env Rscript
# =============================================================================
# create-test-6k-PWAS.R
#
# Build the PWAS evaluation set following the "natural background" enrichment
# design.
#
# Design principle
# ----------------
# A (gene, disease) cell counts as a NEGATIVE only if PWAS actually tested it.
# Cells the method could not analyze (no gene-trait result reachable through
# the lookup) are UNMEASURABLE and are excluded -- treating them as negatives
# would conflate "tested and not significant" with "never tested", which is the
# central flaw the supervisor flagged.
#
# Tested scope, per cell (g, d):
#   PWAS has at least one row in `dat` (post-QC: r2.c >= 0.01) whose
#   (ID.inner = g, TRAIT) maps -- via lookup-MVP traitFromSourceMappedIds --
#   to one of d's ontology IDs from universe-6k's `Disease indication ID`.
#
# LABEL:
#   1 if (g, d) appears in universe-6k (drug-target-indication trio), else 0.
#   The universe is NOT pre-restricted to the PWAS-covered gene/disease axes:
#   pairs that the universe contains but PWAS cannot test simply do not appear
#   in the tested scope (they are unmeasurable, not FN).
#
# Outputs:
#   data/gene-6k-PWAS.txt      -- genes appearing in the tested scope
#   data/disease-6k-PWAS.txt   -- diseases appearing in the tested scope
#   data/universe-6k-PWAS.tsv  -- universe rows whose (gene, disease) pair is
#                                 in the tested scope (i.e. the LABEL = 1 set)
#   test/test-6k-PWAS.tsv      -- one row per tested (gene, disease) cell with
#                                 columns: Target gene name | Disease indication
#                                 | TRAIT (semicolon-joined supporting traits)
#                                 | LABEL
#
# Run from the project root:
#   Rscript code/create-test-6k-PWAS.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------
load("result-compiled/PWAS.RData")              # provides `dat`
setDT(dat)

universe <- fread("data/universe-6k.tsv", sep = "\t", quote = "")
lookup   <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")

# -----------------------------------------------------------------------------
# Step 1: QC-filtered method results in (gene, TRAIT) form
# -----------------------------------------------------------------------------
# r2.c >= 0.01 is the heritability cutoff that defines an "analyzable" gene
# for PWAS. Note: we do NOT restrict to universe genes here -- the tested
# scope must include genes the method analyzed but that are not in universe;
# those cells become genuine LABEL = 0 negatives.
dat.f  <- dat[r2.c >= 0.01]
dat.gt <- unique(dat.f[, .(`Target gene name` = ID.inner, TRAIT)])

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
# Join PWAS (gene, TRAIT) with TRAIT->id, then with universe disease->id, so a
# tested triple (gene, TRAIT, Disease indication) exists iff PWAS produced a
# result for `gene` against a `TRAIT` whose mapped ontology ID matches one of
# `Disease indication`'s IDs.
gene.trait.id <- merge(dat.gt, lookup.long,
                       by = "TRAIT", allow.cartesian = TRUE)
tested.triples <- merge(gene.trait.id, disease.id.long,
                        by = "id", allow.cartesian = TRUE)

# Aggregate to one row per tested (gene, disease) cell; keep the TRAIT(s) that
# made the cell testable (used downstream by benchmark for FDR scope).
test <- tested.triples[, .(
  TRAIT = paste(sort(unique(TRAIT)), collapse = ";")
), by = .(`Target gene name`, `Disease indication`)]

# -----------------------------------------------------------------------------
# Step 5: LABEL against the full universe-6k (deduped on gene x disease name)
# -----------------------------------------------------------------------------
universe.dedup <- unique(universe, by = c("Target gene name", "Disease indication"))
pos.key <- universe.dedup[, paste(`Target gene name`, `Disease indication`,
                                  sep = "||")]
test[, LABEL := as.integer(
  paste(`Target gene name`, `Disease indication`, sep = "||") %in% pos.key
)]

setcolorder(test, c("Target gene name", "Disease indication", "TRAIT", "LABEL"))
fwrite(test, "test/test-6k-PWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Helper outputs
# -----------------------------------------------------------------------------
# Genes / diseases that survived into the tested scope (semantics differ from
# the prior design: these are no longer "method ∩ universe" but "method covers
# a testable cell involving this gene/disease").
genes.tested    <- sort(unique(test[["Target gene name"]]))
diseases.tested <- sort(unique(test[["Disease indication"]]))
writeLines(genes.tested,    "data/gene-6k-PWAS.txt")
writeLines(diseases.tested, "data/disease-6k-PWAS.txt")

# Shrunken universe = positive (LABEL = 1) cells with their original metadata.
pos.pairs <- test[LABEL == 1, .(`Target gene name`, `Disease indication`)]
univ.sub  <- merge(universe.dedup, pos.pairs,
                   by = c("Target gene name", "Disease indication"))
fwrite(univ.sub, "data/universe-6k-PWAS.tsv", sep = "\t", quote = FALSE)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat(sprintf("genes tested:     %d\n", length(genes.tested)))
cat(sprintf("diseases tested:  %d\n", length(diseases.tested)))
cat(sprintf("test cells:       %d  (positives: %d)\n",
            nrow(test), sum(test$LABEL)))
cat(sprintf("universe rows in scope: %d\n", nrow(univ.sub)))
