#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-vs-closest.R
#
# Fixed-FDR set-difference benchmark for each omics method (PWAS/BLISS, MR,
# coloc, TWAS) vs. closest-gene. For each method:
#
#   1. Compute predictions on the headline (full Minikel) universe with FDR
#      adjusted over the method's own full prediction set.
#   2. Restrict the EVALUATION set (not the FDR set) to cells outside the
#      closest-gene prediction set.
#   3. Report contingency, PPV, fold enrichment, and Fisher exact OR + p.
#
# The closest-gene row is reported on its own headline universe for reference
# (it cannot be "vs itself").
#
# Run from the project root:
#   Rscript code/benchmark-13k-vs-closest.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Closest-fire (gene, MESH) cells, from raw loci + TRAIT->MESH lookup.
# These are the cells closest-gene predicts positive on.
# -----------------------------------------------------------------------------
lookup <- fread("dictionary/lookup-MVP.tsv", sep = "\t", quote = "")
trait2mesh <- lookup[nzchar(MESH), .(TRAIT, MESH)]

loci <- fread("data/loci-from-MVP-corrected.csv")
if (colnames(loci)[1] != "Population") setnames(loci, 1, "Population")
loci[, p.num := suppressWarnings(as.numeric(`P-Value`))]
loci.eur <- loci[Population == "EUR" &
                 nzchar(`Lead Variant Gene (Nearest)`) &
                 nzchar(Trait) & !is.na(p.num) & p.num < 5e-8]
closest.pairs <- unique(loci.eur[, .(
  gene = `Lead Variant Gene (Nearest)`, TRAIT = Trait
)])
closest.cells <- unique(merge(
  closest.pairs, trait2mesh, by = "TRAIT", allow.cartesian = TRUE
)[, .(`Target gene name` = gene, MESH)])
closest.cells[, in.closest := TRUE]

# -----------------------------------------------------------------------------
# Helper: take a test set + predicted-positive (gene, MESH) cells and compute
# the contingency for two views:
#   (a) headline (full test set)
#   (b) vs-closest (test set restricted to cells outside closest.cells)
# Returns a data.table with both rows.
# -----------------------------------------------------------------------------
score_views <- function(method.name, test, pred.cells) {
  test <- merge(test, pred.cells, by = c("Target gene name", "MESH"),
                all.x = TRUE)
  test[is.na(predicted), predicted := 0L]

  test <- merge(test, closest.cells,
                by = c("Target gene name", "MESH"), all.x = TRUE)
  test[is.na(in.closest), in.closest := FALSE]

  summarise <- function(d, label) {
    TP <- d[LABEL == 1 & predicted == 1, .N]
    FP <- d[LABEL == 0 & predicted == 1, .N]
    FN <- d[LABEL == 1 & predicted == 0, .N]
    TN <- d[LABEL == 0 & predicted == 0, .N]
    base.rate     <- (TP + FN) / (TP + FP + FN + TN)
    observed.rate <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    fold          <- observed.rate / base.rate
    ft <- tryCatch(
      fisher.test(matrix(c(TP, FP, FN, TN), 2, byrow = TRUE)),
      error = function(e) NULL
    )
    data.table(
      method = method.name, view = label,
      TP = TP, FP = FP, FN = FN, TN = TN,
      total.pred = TP + FP,
      PPV = round(observed.rate, 3),
      base.rate = round(base.rate, 3),
      fold = round(fold, 2),
      OR = if (!is.null(ft)) round(ft$estimate, 2) else NA_real_,
      CI.lo = if (!is.null(ft)) round(ft$conf.int[1], 2) else NA_real_,
      CI.hi = if (!is.null(ft)) round(ft$conf.int[2], 2) else NA_real_,
      p.value = if (!is.null(ft)) ft$p.value else NA_real_
    )
  }

  rbind(
    summarise(test,                                         "headline"),
    summarise(test[in.closest == FALSE],                    "vs-closest")
  )
}

# Score a method: load its dat, build (gene, TRAIT) significance with FDR
# over the test-restricted dat (mirroring the per-method benchmark scripts).
score_method <- function(method.name, test_path, dat_file, dat_var,
                         scoring_fn) {
  test <- fread(test_path, sep = "\t", quote = "")
  e <- new.env()
  load(dat_file, envir = e)
  dat <- get(dat_var, envir = e)
  setDT(dat)
  pred.cells <- scoring_fn(dat, test)
  score_views(method.name, test, pred.cells)
}

# Generic predict_cells from significant (gene, TRAIT) pairs to (gene, MESH).
cells_from_sig <- function(test, sig.pairs) {
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  scored <- merge(exam.long, sig.pairs,
                  by = c("Target gene name", "TRAIT"), all.x = TRUE)
  scored[is.na(sig), sig := FALSE]
  unique(scored[sig == TRUE, .(`Target gene name`, MESH, predicted = 1L)])
}

# -----------------------------------------------------------------------------
# PWAS / BLISS
# -----------------------------------------------------------------------------
score_pwas <- function(dat, test) {
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  exam.pairs <- unique(exam.long[, .(ID.inner = `Target gene name`, TRAIT)])
  dat.sub <- merge(dat, exam.pairs, by = c("ID.inner", "TRAIT"))
  dat.sub[, q := p.adjust(p.a.c, method = "fdr")]
  dat.sub[, sig := !is.na(q) & q < 0.05]
  sig.pairs <- unique(dat.sub[sig == TRUE,
                              .(`Target gene name` = ID.inner, TRAIT,
                                sig = TRUE)])
  cells_from_sig(test, sig.pairs)
}

# -----------------------------------------------------------------------------
# MR
# -----------------------------------------------------------------------------
olink <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

score_mr <- function(dat, test) {
  dat <- dat[exception %in% c("blank", "ONLY_ONE_IV", "SUMSTAT_NO_SIGNAL")]
  dat <- merge(dat, prot2gene, by = "protein")
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  exam.pairs <- unique(exam.long[, .(`Target gene name`, TRAIT)])
  dat.sub <- merge(dat, exam.pairs, by = c("Target gene name", "TRAIT"))
  dat.sub[, p.num := suppressWarnings(as.numeric(p))]
  dat.sub[is.na(p.num), p.num := 1]
  dat.sub[, q := p.adjust(p.num, method = "fdr")]
  dat.sub[, sig := q < 0.05]
  sig.pairs <- unique(dat.sub[sig == TRUE,
                              .(`Target gene name`, TRAIT, sig = TRUE)])
  cells_from_sig(test, sig.pairs)
}

# -----------------------------------------------------------------------------
# coloc (H4 >= 0.8 as significance; no FDR)
# -----------------------------------------------------------------------------
score_coloc <- function(dat, test) {
  dat <- merge(dat, prot2gene, by = "protein")
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  exam.pairs <- unique(exam.long[, .(`Target gene name`, TRAIT)])
  dat.sub <- merge(dat, exam.pairs, by = c("Target gene name", "TRAIT"))
  dat.sub[, H4.num := suppressWarnings(as.numeric(H4))]
  dat.sub[is.na(H4.num), H4.num := 0]
  dat.sub[, sig := H4.num >= 0.8]
  sig.pairs <- unique(dat.sub[sig == TRUE,
                              .(`Target gene name`, TRAIT, sig = TRUE)])
  cells_from_sig(test, sig.pairs)
}

# -----------------------------------------------------------------------------
# TWAS
# -----------------------------------------------------------------------------
score_twas <- function(dat, test) {
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  exam.pairs <- unique(exam.long[, .(gene_name = `Target gene name`, TRAIT)])
  dat.sub <- merge(dat, exam.pairs, by = c("gene_name", "TRAIT"))
  dat.sub[, q := p.adjust(pvalue, method = "fdr")]
  dat.sub[, sig := !is.na(q) & q < 0.05]
  sig.pairs <- unique(dat.sub[sig == TRUE,
                              .(`Target gene name` = gene_name, TRAIT,
                                sig = TRUE)])
  cells_from_sig(test, sig.pairs)
}

# -----------------------------------------------------------------------------
# closest (predicted from raw loci; no FDR). Reported on its own test only.
# -----------------------------------------------------------------------------
score_closest <- function(dat = NULL, test) {
  loci.tab <- fread("data/loci-from-MVP-corrected.csv")
  if (colnames(loci.tab)[1] != "Population") setnames(loci.tab, 1, "Population")
  loci.tab[, p.num := suppressWarnings(as.numeric(`P-Value`))]
  dat.C <- unique(loci.tab[
    Population == "EUR" &
      nzchar(`Lead Variant Gene (Nearest)`) &
      nzchar(Trait) & !is.na(p.num) & p.num < 5e-8,
    .(`Target gene name` = `Lead Variant Gene (Nearest)`,
      TRAIT              = Trait)
  ])
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  exam.pairs <- unique(exam.long[, .(`Target gene name`, TRAIT)])
  dat.sub <- merge(dat.C, exam.pairs, by = c("Target gene name", "TRAIT"))
  sig.pairs <- unique(dat.sub[, .(`Target gene name`, TRAIT, sig = TRUE)])
  cells_from_sig(test, sig.pairs)
}

# -----------------------------------------------------------------------------
# Run everything
# -----------------------------------------------------------------------------
results <- rbindlist(list(
  score_method("BLISS",         "test/test-13k-PWAS.tsv",
               "result-compiled/PWAS.RData",  "dat", score_pwas),
  score_method("MR",            "test/test-13k-MR.tsv",
               "result-compiled/MR.RData",    "dat", score_mr),
  score_method("Colocalization", "test/test-13k-coloc.tsv",
               "result-compiled/coloc.RData", "dat", score_coloc),
  score_method("TWAS",          "test/test-13k-TWAS.tsv",
               "result-compiled/TWAS.RData",  "dat", score_twas)
))

# closest on its own test, headline view only
test.C <- fread("test/test-13k-closest.tsv", sep = "\t", quote = "")
pred.C <- score_closest(NULL, test.C)
test.C2 <- merge(test.C, pred.C, by = c("Target gene name", "MESH"),
                 all.x = TRUE)
test.C2[is.na(predicted), predicted := 0L]
TP <- test.C2[LABEL == 1 & predicted == 1, .N]
FP <- test.C2[LABEL == 0 & predicted == 1, .N]
FN <- test.C2[LABEL == 1 & predicted == 0, .N]
TN <- test.C2[LABEL == 0 & predicted == 0, .N]
ft <- fisher.test(matrix(c(TP, FP, FN, TN), 2, byrow = TRUE))
closest.row <- data.table(
  method = "closest gene", view = "headline",
  TP = TP, FP = FP, FN = FN, TN = TN,
  total.pred = TP + FP,
  PPV = round(TP / (TP + FP), 3),
  base.rate = round((TP + FN) / (TP + FP + FN + TN), 3),
  fold = round((TP / (TP + FP)) / ((TP + FN) / (TP + FP + FN + TN)), 2),
  OR = round(ft$estimate, 2),
  CI.lo = round(ft$conf.int[1], 2),
  CI.hi = round(ft$conf.int[2], 2),
  p.value = ft$p.value
)
results <- rbind(results, closest.row)

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
cat("\n=== Benchmark: omics methods evaluated on headline universe and on\n",
    "    closest-gene's complement (fixed-FDR set difference). ===\n\n",
    sep = "")
setorder(results, method, view)
fmt <- function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                          formatC(p, format = "g", digits = 2))
out <- results[, .(
  method, view,
  `Val/Total` = sprintf("%d/%d", TP, total.pred),
  PPV  = sprintf("%.0f%%", 100 * PPV),
  Fold = sprintf("%.2fx", fold),
  `OR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", OR, CI.lo, CI.hi),
  P = sapply(p.value, fmt)
)]
print(out, nrows = nrow(out))

fwrite(results, "result-compiled/benchmark-vs-closest.tsv",
       sep = "\t", quote = FALSE)
cat("\nWrote: result-compiled/benchmark-vs-closest.tsv\n")
