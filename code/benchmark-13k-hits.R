#!/usr/bin/env Rscript
# =============================================================================
# benchmark-13k-hits.R
#
# Dump the predicted-positive (gene, MESH) cells for each of the 4 omics
# methods (BLISS, MR, Coloc, TWAS) plus the closest-gene reference, using the
# same scoring as benchmark-13k-vs-closest.R. Output: a long table suitable
# for UpSet plotting.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source_env <- new.env()
# Reuse scoring fns from the vs-closest script by sourcing it in an env, but
# stop before it runs the big rbindlist at the bottom. Simpler to redefine
# the small bits we need here.

olink <- fread("dictionary/lookup-Olink.tsv", sep = "\t", quote = "")
prot2gene <- unique(olink[, .(protein = OlinkID,
                              `Target gene name` = HGNC.symbol)])
prot2gene <- prot2gene[nzchar(`Target gene name`)]

cells_from_sig <- function(test, sig.pairs) {
  exam.long <- test[nzchar(TRAIT), .(
    TRAIT = unlist(strsplit(TRAIT, ";"))
  ), by = .(`Target gene name`, MESH, LABEL)]
  exam.long[, TRAIT := trimws(TRAIT)]
  exam.long <- exam.long[nzchar(TRAIT)]
  scored <- merge(exam.long, sig.pairs,
                  by = c("Target gene name", "TRAIT"), all.x = TRUE)
  scored[is.na(sig), sig := FALSE]
  unique(scored[sig == TRUE, .(`Target gene name`, MESH)])
}

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

score_closest_cells <- function(test) {
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
# Run each method: load its test set, score, attach LABEL.
# -----------------------------------------------------------------------------
run <- function(method.name, test_path, dat_file, dat_var, scoring_fn) {
  test <- fread(test_path, sep = "\t", quote = "")
  e <- new.env()
  load(dat_file, envir = e)
  dat <- get(dat_var, envir = e)
  setDT(dat)
  cells <- scoring_fn(dat, test)
  cells <- merge(cells, unique(test[, .(`Target gene name`, MESH, LABEL)]),
                 by = c("Target gene name", "MESH"), all.x = TRUE)
  cells[, method := method.name]
  cells[, .(method, `Target gene name`, MESH, LABEL)]
}

hits <- rbindlist(list(
  run("BLISS",          "test/test-13k-PWAS.tsv",
      "result-compiled/PWAS.RData",  "dat", score_pwas),
  run("MR",             "test/test-13k-MR.tsv",
      "result-compiled/MR.RData",    "dat", score_mr),
  run("Colocalization", "test/test-13k-coloc.tsv",
      "result-compiled/coloc.RData", "dat", score_coloc),
  run("TWAS",           "test/test-13k-TWAS.tsv",
      "result-compiled/TWAS.RData",  "dat", score_twas)
))

# closest-gene reference on its own test
test.C <- fread("test/test-13k-closest.tsv", sep = "\t", quote = "")
cells.C <- score_closest_cells(test.C)
cells.C <- merge(cells.C, unique(test.C[, .(`Target gene name`, MESH, LABEL)]),
                 by = c("Target gene name", "MESH"), all.x = TRUE)
cells.C[, method := "closest"]
hits <- rbind(hits, cells.C[, .(method, `Target gene name`, MESH, LABEL)])

setorder(hits, method, `Target gene name`, MESH)

fwrite(hits, "result-compiled/benchmark-13k-hits-long.tsv",
       sep = "\t", quote = FALSE)

# Also a wide indicator matrix on the union of all (gene, MESH) cells called
# by any method — convenient for UpSet packages that want a 0/1 matrix.
all.cells <- unique(hits[, .(`Target gene name`, MESH)])
for (m in unique(hits$method)) {
  hit.m <- hits[method == m, .(`Target gene name`, MESH)]
  hit.m[, (m) := 1L]
  all.cells <- merge(all.cells, hit.m,
                     by = c("Target gene name", "MESH"), all.x = TRUE)
}
# LABEL from the union of test sets (any test marking the cell LABEL==1)
label.union <- rbindlist(lapply(
  list("test/test-13k-PWAS.tsv","test/test-13k-MR.tsv",
       "test/test-13k-coloc.tsv","test/test-13k-TWAS.tsv",
       "test/test-13k-closest.tsv"),
  function(p) fread(p, sep="\t", quote="")[, .(`Target gene name`, MESH, LABEL)]
))
label.max <- label.union[, .(LABEL = max(LABEL, na.rm = TRUE)),
                         by = .(`Target gene name`, MESH)]
all.cells <- merge(all.cells, label.max,
                   by = c("Target gene name", "MESH"), all.x = TRUE)
for (col in setdiff(names(all.cells), c("Target gene name","MESH","LABEL"))) {
  all.cells[is.na(get(col)), (col) := 0L]
}
setcolorder(all.cells, c("Target gene name","MESH","LABEL",
                         "BLISS","MR","Colocalization","TWAS","closest"))
setorder(all.cells, -LABEL, `Target gene name`, MESH)

fwrite(all.cells, "result-compiled/benchmark-13k-hits-wide.tsv",
       sep = "\t", quote = FALSE)

# Quick sanity printout: per-method hit counts and TP counts
summary <- hits[, .(total.hits = .N, validated = sum(LABEL == 1, na.rm=TRUE)),
                by = method]
cat("--- Per-method hit counts (should match headline column) ---\n")
print(summary)
cat("\nWrote:\n",
    "  result-compiled/benchmark-13k-hits-long.tsv  (one row per method-hit)\n",
    "  result-compiled/benchmark-13k-hits-wide.tsv  (indicator matrix)\n",
    sep = "")
