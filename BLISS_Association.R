# Install required packages and load libraries
cat("Starting PWAS analysis...\n")
options(repos = c(CRAN = "https://cloud.r-project.org/"))

package <- c("data.table", "dplyr", "optparse")
package <- package[!(package %in% installed.packages()[, "Package"])]

if (length(package)) {
    cat("Installing packages:", paste(package, collapse = ", "), "\n")
    suppressMessages(install.packages(package, quiet = TRUE))
}

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

#################
# Parse options #
#################

option_list <- list(
    make_option("--path_sumstats", type = "character", default = NA, action = "store", help = "Path to the processed GWAS summary dataset file."
    ),
    make_option("--n", type = "numeric", default = NA, action = "store", help = "Sample size of GWAS summary dataset."
    ),
    make_option("--model", type = "character", default = "UKBPPP_EUR", action = "store", help = "Protein prediction models to be used."
    ),
    make_option("--chr", type = "numeric", default = NA, action = "store", help = "Chromosome to work on. When unspecified, will process chromosome 1-22."
    ),
    make_option("--output_dir", type = "character", default = getwd(), action = "store", help = "Directory to stor output file. Default is current directory."
    ),
    make_option("--output_name", type = "character", default = "PWAS_results", action = "store", help = "Name of the output file"
    ),
    make_option("--output_augmented", type = "logical", default = FALSE, action = "store", help = "Whether to output the augmented results. Default is FALSE."
    ),
    make_option("--output_twas_fusion", type = "logical", default = FALSE, action = "store", help = "Format the final output like TWAS-fusion. Default is FALSE."
    ),
    make_option("--clean_slate", type = "logical", default = FALSE, action = "store", help = "Whether to clean the slate and remove existing results before running the analysis. Default is FALSE."
    )
)

# Parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

ss.path     <- opt$path_sumstats
n.sumstats  <- opt$n
model       <- opt$model
CHR         <- opt$chr
output.dir  <- opt$output_dir
output.name <- opt$output_name
output.aug  <- opt$output_augmented
output.tf   <- opt$output_twas_fusion
clean.slate <- opt$clean_slate

cat("Arguments: model =", model, "| chr =", ifelse(is.na(CHR), "1-22", paste(CHR, collapse = ",")), 
    "| output =", output.dir, output.name, "\n")

# Fix CHR
if (is.na(CHR)) {
    CHR <- 1:22
}

# Ensure that the summary statistics file exists
if (!file.exists(ss.path)) {
    stop(paste0("Summary statistics file ", ss.path, " not found!\n"))
}

# Ensure that the model files exist
if (paste0("model/", model) %>% dir.exists()) {
    if (paste0("model/", model, "/.manifest") %>% file.exists()) {
        paste0("Model ", model, " found.\nCorresponding manifest file found.\n") %>% cat()
    } else {
        stop(
            paste0("Model ", model, " is found, but with no manifest file!\n", "Please make sure the proper model files are in ", getwd(), ".")
        )
    }
} else {
    stop(
        paste0("Model ", model, " not found!\n", "Please make sure the proper model files are in ", getwd(), ".")
    )
}

########################################
# Old functions from TestAssociation.r #
########################################

# ComputeAlpha
ComputeAlpha <- function(w, Z, n, n_0, matrix.LD) {
    # BIG SIGMA
    SIGMA <- t(w) %*% matrix.LD %*% w * (n_0 - 1) / n_0

    # SIGMA <- (t(ref %*% w) %*% (ref %*% w)) / n_0
    SIGMA <- SIGMA[1, 1]

    # Alpha hat
    alpha <- ((w %*% Z) / sqrt(n)) / SIGMA
    alpha <- alpha[1, 1]

    # small sigma
    sigma <- sqrt(1 - 2 * ((w %*% Z) / sqrt(n)) * alpha + alpha ^ 2 * SIGMA)
    sigma <- sigma[1, 1]

    # se of alpha
    se <- sqrt((1 / n + 1 / n_0) * (SIGMA * alpha ^ 2) / SIGMA + sigma ^ 2 / (n * SIGMA))

    return(c(alpha, se))
}

# allele.qc
allele.qc <- function(a1, a2, ref1, ref2) {
    ref <- ref1
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip1 <- flip
    ref <- ref2
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip2 <- flip
    snp <- list()
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}

# TestAssociation
TestAssociation <- function(sumstats, weight, SS.original, matrix.LD, n.sumstats, n.ref) {

    # Find common SNPs
    snps.common <- intersect(sumstats$SNP, SS.original$SNP)

    # Keep only common SNPs
    sumstats.temp <- sumstats[sumstats$SNP %in% snps.common, ]

    rownames(sumstats.temp) <- sumstats.temp$SNP
    names(weight) <- SS.original$SNP

    sumstats.temp <- sumstats.temp[snps.common, ]
    SS.original <- SS.original[snps.common, ]
    weight <- weight[snps.common]
    matrix.LD <- matrix.LD[snps.common, snps.common]

    # Allele-flip the phenotype ss and weight w.r.t. reference panel
    qc <- allele.qc(
        sumstats.temp$A1,
        sumstats.temp$A2,
        SS.original$A1,
        SS.original$A2
    )

    sumstats.temp$Z[qc$flip] <- -1 * sumstats.temp$Z[qc$flip]

    if (!("N" %in% colnames(sumstats.temp))) {
        sumstats.temp["N"] <- n.sumstats
    } else {
        n.sumstats <- floor(mean(sumstats.temp$N))
    }

    # Remove strand ambiguous SNPs (if any)
    if (sum(!qc$keep) > 0) {
        sumstats.temp <- sumstats.temp[keep, ]
        SS.original   <- SS.original[keep, ]
        weight        <- weight[keep]
        matrix.LD     <- matrix.LD[keep, keep]
    }

    ################
    # Classic PWAS #
    ################

    # Compute TWAS z-score, r2, and p-value
    z.twas  <- as.numeric(weight %*% sumstats.temp$Z)
    r2.twas <- as.numeric(weight %*% matrix.LD %*% weight)

    p.classic <- 2 * (pnorm(abs(z.twas / sqrt(r2.twas)), lower.tail = FALSE))

    ####################
    # Alternative PWAS #
    ####################
    
    if (length(weight) > 1) {
        out <- ComputeAlpha(
            w   = weight,
            Z   = sumstats.temp$Z,
            n   = n.sumstats,
            n_0 = n.ref,
            matrix.LD = matrix.LD
        )

        out[3] <- 2 * (pnorm(abs(out[1] / out[2]), lower.tail = FALSE))
    } else {
        out <- c(NA, NA, NA)
    }

    BEST.ID <- sumstats.temp$SNP[which.max(abs(sumstats.temp$Z))[1]]
    BEST.Z  <- sumstats.temp$Z[which.max(abs(sumstats.temp$Z))[1]]

    output <- c(z.twas, sqrt(r2.twas), p.classic, out, BEST.ID, BEST.Z)

    names(output) <- c(
        "beta.classic","se.classic", "p.classic",
        "beta_alt", "se_alt", "p_alt",
        "BEST.GWAS.ID", "BEST.GWAS.Z"
    )

    return(output)
}

recover_corr_matrix <- function(stored_data) {
    # Input validation
    if (!is.list(stored_data)) {
        stop("Input must be a list (output from store_upper_triangle)")
    }

    required_names <- c("upper_values", "n_vars", "var_names")
    if (!all(required_names %in% names(stored_data))) {
        stop("Input list must contain: upper_values, n_vars, var_names")
    }

    # Extract components
    upper_values <- stored_data$upper_values
    n_vars <- stored_data$n_vars
    var_names <- stored_data$var_names

    # Validate dimensions
    expected_length <- n_vars * (n_vars - 1) / 2
    if (length(upper_values) != expected_length) {
        stop(paste("Expected", expected_length, "upper triangle values for", 
                  n_vars, "variables, but got", length(upper_values)))
    }

    # Initialize matrix with zeros
    corr_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)

    # Fill upper triangle
    corr_matrix[upper.tri(corr_matrix)] <- upper_values

    # Make symmetric by copying upper triangle to lower triangle
    corr_matrix[lower.tri(corr_matrix)] <- t(corr_matrix)[lower.tri(corr_matrix)]

    # Set diagonal to 1 (assuming correlation matrix)
    diag(corr_matrix) <- 1

    # Add variable names if available
    if (!is.null(var_names) && length(var_names) == n_vars) {
        rownames(corr_matrix) <- var_names
        colnames(corr_matrix) <- var_names
    }

    return(corr_matrix)
}

########################################
# Functions used to handle file states #
########################################

# GetState
GetState <- function(destination) {
    if (destination %>% paste0(., ".finished") %>% file.exists()) {
        state <- 2
    } else if (destination %>% file.exists()) {
        state <- 1
    } else {
        state <- 0
    }

    return(state)
}

# InitializeFile
InitializeFile <- function(state.current, destination, header) {
    if (state.current == 0) {
        file.current <- matrix(integer(0), nrow = 0, ncol = length(header)) %>% as.data.frame()

        names(file.current) <- header

        write.table(
            file.current,
            file = destination,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    } else {
        file.current <- fread(destination, data.table = FALSE)
    }

    return(file.current)
}

# MarkFinished
MarkFinished <- function(destination) {
    paste0("mv ", destination, " ", destination, ".finished") %>% system()

    invisible(NULL)
}

############
# Manifest #
############

# Read the manifest file
manifest <- paste0("model/", model, "/.manifest") %>%
    fread(., data.table = FALSE) %>%
    filter(chromosome %in% CHR)

cat("Found", nrow(manifest), "protein models for chromosome(s):", paste(CHR, collapse = ", "), "\n")

if (nrow(manifest) == 0) {
    stop(paste0("No protein models found for chromosome(s): ", paste(CHR, collapse = ", "), "."))
}

######################
# Summary statistics #
######################

# Trial read the summary statistics file
ss <- fread(ss.path, data.table = FALSE, header = TRUE, showProgress = FALSE, nrows = 10)

if (is.na(n.sumstats)) {
    if (!("N" %in% colnames(ss))) {
        stop("The summary statistics file does not contain a column named 'N'. Please provide the sample size with --n.")
    }

    if (all(is.na(ss$N))) {
        stop("The summary statistics file contains no valid sample size information. Please provide the sample size with --n.")
    }
}

# Read the summary statistics file
# Be cautious with na.omit() as it may remove too many rows unnecessarily
cat("Loading summary statistics...\n")
ss.pheno <- fread(ss.path, data.table = FALSE, header = TRUE, showProgress = FALSE)
ss.pheno <- ss.pheno[ss.pheno$CHR %in% CHR, ] %>%
    na.omit() %>%
    filter(!duplicated(SNP))

cat("Loaded", nrow(ss.pheno), "SNPs after filtering\n")

if (is.na(n.sumstats)) {
    n.sumstats <- median(ss.pheno$N, na.rm = TRUE)
    cat("Sample size:", n.sumstats, "\n")
}

##################
# Main iteration #
##################

# Create the output directory if it does not exist
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

destination <- file.path(output.dir, output.name)

# Header
if (!output.tf) {
    header <- c(
        "protein", "beta", "se", "p", "model_r2", "n_snps", "MHC"
    )
} else {
    header <- c(
        "FILE", "ID", "CHR", "P0", "P1", "HSQ", "BEST.GWAS.ID", "BEST.GWAS.Z",
        "EQTL.ID", "EQTL.R2", "EQTL.Z", "EQTL.GWAS.Z", "NSNP", "MODEL", "MODELCV.R2",
        "MODELCV.PV", "TWAS.Z", "TWAS.P"
    )

    output.aug <- FALSE
}


# Get state
state.current <- GetState(destination)

if (state.current != 2) {
    file.current <- InitializeFile(state.current, destination, header)
} else {
    if (!clean.slate) {
        stop(paste0("The file ", destination, " is already finished. Please remove the finished file to re-run the analysis or set clean_slate to TRUE.\n"))
    } else {
        paste0(destination, ".finished") %>% file.remove()
        file.current <- InitializeFile(0, destination, header)
    }
}

# Wipe the file if clean slate is requested
if (clean.slate) {
    if (file.exists(destination)) {
        file.remove(destination)
    }

    file.current <- InitializeFile(0, destination, header)
}

# Get index
if (nrow(file.current) == 0) {
    index.current <- 1
} else {
    index.current <- which(manifest$protein == file.current$protein[nrow(file.current)]) + 1
    cat("Resuming from protein", index.current, "\n")
}

# Setup progress bar
total <- nrow(manifest)
to_do <- total - (index.current - 1)
cat("Processing", to_do, "proteins...\n")
pb <- txtProgressBar(min = 0, max = to_do, style = 3)

# Main iteration
for (j in index.current:nrow(manifest)) {

    current_progress <- j - index.current + 1
    setTxtProgressBar(pb, current_progress)

    # Initialize the update vector
    update <- length(header) %>% numeric()

    # 1: protein name
    update[1] <- manifest$protein[j]

    file <- paste0("model/", model, "/", manifest$filename[j])
    if (file %>% file.exists()) {
        load(file)

        matrix.LD <- recover_corr_matrix(matrix.LD)
        if (all(CHECKED$alpha == 0)) {
            CHECKED$alpha <- rep(1 / length(CHECKED$alpha), length(CHECKED$alpha))
        }

        weight <- as.matrix(weight) %*% CHECKED$alpha

        tryCatch(
            expr = {
                temp <- TestAssociation(
                    sumstats    = ss.pheno,
                    weight      = weight,
                    SS.original = ss,
                    matrix.LD   = matrix.LD,
                    n.sumstats  = n.sumstats,
                    n.ref       = manifest$n_ref[j]
                )
            },
            error = function(e) {
                cat("\nError processing", manifest$protein[j], ":", conditionMessage(e), "\n")
                temp <<- rep(NA, 6)
            }
        )

        update[2:4] <- temp[4:6]
        update[5] <- mean(R2.PUMAS, na.rm = TRUE)
        update[6] <- length(weight)
        update[7] <- manifest$MHC[j]

        if (output.tf) {
            # Re-do update
            update[1] <- file
            update[2] <- manifest$protein[j]
            update[3] <- manifest$chromosome[j]
            update[4] <- manifest$start[j]
            update[5] <- manifest$end[j]
            update[6] <- mean(R2.PUMAS, na.rm = TRUE)
            update[7] <- temp[7]
            update[8] <- temp[8]
            update[9] <- NA
            update[10] <- NA
            update[11] <- NA
            update[12] <- NA
            update[13] <- length(weight)
            update[14] <- "BLISS"
            update[15] <- NA
            update[16] <- NA
            update[17] <- as.numeric(temp[4]) / as.numeric(temp[5])
            update[18] <- temp[6]
        }
    }

    # Update
    cat(
        update %>% paste(., collapse = "\t") %>% paste0(., "\n"),
        file = destination,
        append = TRUE
    )
}

# Close progress bar
close(pb)

cat(sprintf("\nCompleted processing %d proteins!\n", to_do))

# Augment results if requested
if (output.aug) {
    cat("Augmenting results...\n")
    # Read the final results and augment with additional information from the manifest
    result <- fread(destination, data.table = FALSE, header = TRUE)

    manifest <- manifest[, c("protein", "filename", "chromosome", "start", "end")]

    result <- left_join(result, manifest, by = "protein")
    result["q"] <- p.adjust(result$p, method = "fdr")

    # Save the augmented results
    write.table(
        result,
        file = destination,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
}

# Rename the intermediate file to "finished"
MarkFinished(destination)
cat("Analysis complete! Results saved to:", paste0(destination, ".finished"), "\n")
