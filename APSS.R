APSS <- function(directory.working, filename, auto = FALSE, do.return = FALSE, BIG = 2) {
    # Packages
    suppressMessages(library(data.table))
    suppressMessages(library(magrittr))
    suppressMessages(library(R.utils))
    suppressMessages(library(stringr))

    # Length of separator
    length.separator <- 120
    separator <- rep("=", length.separator) %>% paste0(., collapse = "")

    # Set working directory
    setwd(directory.working)

    # Read the input file
    check.BIG <- FALSE
    if (file.size(filename) >= BIG * 1024 * 1024 * 1024) {
        check.BIG <- TRUE
        cat(paste0("Input dataset is VERY VERY VERY BIG. We have to be selective :-<<<"), sep = "\n")
        cat(separator, sep = "\n")
        raw.try <- fread(filename, showProgress = FALSE, nrows = 100) %>% as.data.frame()

        # Print the head of the input dataset
        cat("Following are the first few lines of the dataset:", sep = "\n")
        cat(separator, sep = "\n")
        print(head(raw.try))
        cat(separator, sep = "\n")

        # Ask user to pick some columns 
        cat("Press ENTER to skip a column! Write something if you want to keep a column. ", sep = "\n")
        cat(separator, sep = "\n")
        keep <- character(ncol(raw.try))
        for (index in 1:ncol(raw.try)) {
            keep[index] <- readline(prompt = paste0("Keep ", colnames(raw.try)[index], "? "))
        }
        cat(separator, sep = "\n")
        keep <- which(keep != "")

        rm(raw.try)

        raw <- fread(filename, showProgress = FALSE, select = keep, fill = TRUE) %>% as.data.frame()
    } else {
        raw <- fread(filename, showProgress = FALSE, fill = TRUE) %>% as.data.frame()
    }

    if (check.BIG) {
        cat(paste0(nrow(raw), " lines read from the selected raw summary statistics file!"), sep = "\n")
    } else {
        cat(paste0(nrow(raw), " lines read from the raw summary statistics file!"), sep = "\n")
    }
    cat(separator, sep = "\n")

    # Get the header
    header.inner <- colnames(raw)

    # Auto mode or not?
    if ("V1" %in% header.inner | "V2" %in% header.inner | "V3" %in% header.inner & auto) {
        cat("I don't see any header! Switch to interactive mode!", sep = "\n")
        cat(separator, sep = "\n")
    }

    # Start interacting
    if (auto == FALSE) {
        # Print the head of the input data set
        cat("Following are the first few lines of the dataset:", sep = "\n")
        cat(separator, sep = "\n")
        print(head(raw))
        cat(separator, sep = "\n")

        # Ask user for reasonable column names
        cat("Feel free to input whatever column names you want as long as they are unique! Press ENTER to skip a column!", sep = "\n")
        for (index in 1:ncol(raw)) {
            header.inner[index] <- readline(prompt = paste0("What would be a proper column name for ", colnames(raw)[index], "? "))
        }
        cat(separator, sep = "\n")
    }

    # Make a copy of the desired header
    header.outer <- header.inner

    ######################
    # Process the header #
    ######################

    # Initialize the header
    header.inner <- tolower(header.inner)

    # SNP
    try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps","rsids")
    header.inner[header.inner %in% try.snp] <- "SNP"

    # A1
    try.a1 <- c("a1", "allele1", "allele_1", "nea", "non_effect_allele", "reference_allele", "inc_allele", "ref", "a1lele1", "al1ele1")
    header.inner[header.inner %in% try.a1] <- "A1"

    # A2
    try.a2 <- c("a2", "allele2", "allele_2", "other_allele",  "dec_allele", "alt", "a0","effect_allele","ea")
    header.inner[header.inner %in% try.a2] <- "A2"

    # Z-score
    try.z <- c("zscore", "z-score", "gc_zscore", "z")
    header.inner[header.inner %in% try.z] <- "Z"

    # P
    try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p")
    header.inner[header.inner %in% try.p] <- "P"

    # Beta
    try.beta <- c("b", "beta", "effects", "effect")
    header.inner[header.inner %in% try.beta] <- "BETA"

    # Odds ratio
    try.or <- c("or")
    header.inner[header.inner %in% try.or] <- "ODDS_RATIO"

    # Log odds
    try.logodds <- c("log_odds", "logor", "log_or")
    header.inner[header.inner %in% try.logodds] <- "LOG_ODDS"

    # MAF
    try.maf <- c("eaf", "frq", "maf", "frq_u", "f_u", "freq","af_alt")
    header.inner[header.inner %in% try.maf] <- "MAF"

    # INFO
    try.info <- c("info", "info_score")
    header.inner[header.inner %in% try.info] <- "INFO"

    # Chromosome
    try.chromosome <- c("chrom", "ch", "chr", "chromosome","#chrom")
    header.inner[header.inner %in% try.chromosome] <- "CHROMOSOME"

    # Position
    try.position <- c("pos", "posit", "position", "bp", "bpos")
    header.inner[header.inner %in% try.position] <- "POSITION"

    # Standard error
    try.se <- c("se", "sebeta", "beta_se")
    header.inner[header.inner %in% try.se] <- "SE"

    # CHR plus POSITION
    try.chrpos <- c("chr:pos", "chr_pos", "chr-pos", "chrpos")
    header.inner[header.inner %in% try.chrpos] <- "CHR_POS"

    # Samplesize
    try.samplesize <- c("n", "samplesize", "num_samples", "sample")
    header.inner[header.inner %in% try.samplesize] <- "SAMPLESIZE"

    # Update the header
    colnames(raw) <- header.inner

    # Drop some rows
    n.start <- ncol(raw)

    raw <- raw[, which(colnames(raw) != "")]
    header.outer <- header.outer[header.outer != ""]
    header.inner <- header.inner[header.inner != ""]

    n.end <- ncol(raw)
    cat(paste0(n.start - n.end, " columns were dropped from the input dataset!"), sep = "\n")
    cat(separator, sep = "\n")

    # Double-check the class of each column and coerce if needed
    list.coerce <- c("Z", "P", "BETA", "ODDS_RATIO", "LOG_ODDS", "MAF", "INFO", "SE")

    options(warn = -1)
    if ("POSITION" %in% header.inner) {
        if (class(raw$POSITION) == "character") {
            class(raw$POSITION) <- "integer"
            cat(paste0("Column POSITION has wrong class and has been coerced to integer."), sep = "\n")
            cat(separator, sep = "\n")
        }
    }

    for (i in 1:length(header.inner)) {
        if (header.inner[i] %in% list.coerce) {
            if (class(raw[, header.inner[i]]) != "numeric") {
                class(raw[, header.inner[i]]) <- "numeric"
                cat(paste0("Column ", header.inner[i], " has wrong class and has been coerced to numeric."), sep = "\n")
                cat(separator, sep = "\n")
            }
        }
    }
    options(warn = 0)

    # Drop rows with missing values
    do.missing <- readline(prompt = "Drop rows with missing values? Y/N ")

    if (do.missing == "Y") {
        n.start <- nrow(raw)

        raw <- raw[complete.cases(raw), ]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")

        if ((n.start - n.end) / n.start >= 0.8) {
            cat(paste0(round((n.start - n.end) / n.start * 100, digits = 2), "% of rows were removed!"), sep = "\n")
            cat("Please check the dataset manually because I am sensing something fishy.")
            Sys.sleep(5)
            .Internal(.invokeRestart(list(NULL, NULL), NULL))
        }
        cat(separator, sep = "\n")
    }

    # Split CHR_POS
    if ("CHR_POS" %in% header.inner) {
        do.split <- readline(prompt = paste0("Do you want me to split the ID column? (This could take a while.) Y/N ")) %>% toupper()
        cat(separator, sep = "\n")
        if (do.split == "Y") {
            split.temp <- str_split(raw$CHR_POS, ":") %>% unlist()
            if (sum(grepl("chr", raw$CHR_POS)) != 0) {
                split.temp[2 * (1:nrow(raw)) - 1] <- gsub("chr", "", split.temp[2 * (1:nrow(raw)) - 1])
            }

            split.temp <- as.numeric(split.temp)
            split.temp <- matrix(split.temp, nrow = nrow(raw), ncol = 2, byrow = TRUE) %>% as.data.frame()
            names(split.temp) <- c("V1", "V2")

            header.inner <- c(header.inner, "CHROMOSOME", "POSITION")
            raw["CHROMOSOME"] <- split.temp[, 1]
            raw["POSITION"]   <- split.temp[, 2]

            # Print the head of split.temp
            cat("Following are the first few lines of the splitted column:", sep = "\n")
            cat(separator, sep = "\n")
            print(head(split.temp))
            cat(separator, sep = "\n")

            # Ask user for reasonable column names
            cat("Feel free to input whatever column names you want as long as they are unique! Press ENTER to skip a column!", sep = "\n")
            for (index in 1:ncol(split.temp)) {
                header.outer[length(header.outer) + 1] <- readline(prompt = paste0("What would be a proper column name for ", colnames(split.temp)[index], "? "))
            }
            cat(separator, sep = "\n")
        }
    }

    # Problematic rsid
    test.1 <- grepl(":", raw$SNP)
    test.2 <- grepl("_", raw$SNP)
    test.3 <- grepl("chr", raw$SNP)
    if (sum(test.1 & test.3) == nrow(raw)) {
        cat("I noticed a pattern (chrX:position) in SNPs' IDs!", sep = "\n")
        do.split <- readline(prompt = paste0("Do you want me to split the ID column? (This could take a while.) Y/N "))
        cat(separator, sep = "\n")
        if (do.split == "Y") {
            split.temp <- str_split(raw$SNP, ":") %>% unlist()
            split.temp[2 * (1:nrow(raw)) - 1] <- gsub("chr", "", split.temp[2 * (1:nrow(raw)) - 1])

            split.temp <- as.numeric(split.temp)
            split.temp <- matrix(split.temp, nrow = nrow(raw), ncol = 2, byrow = TRUE) %>% as.data.frame()
            names(split.temp) <- c("V1", "V2")

            header.inner <- c(header.inner, "CHROMOSOME", "POSITION")
            raw["CHROMOSOME"] <- split.temp[, 1]
            raw["POSITION"]   <- split.temp[, 2]

            # Print the head of split.temp
            cat("Following are the first few lines of the splitted column:", sep = "\n")
            cat(separator, sep = "\n")
            print(head(split.temp))
            cat(separator, sep = "\n")

            # Ask user for reasonable column names
            cat("Feel free to input whatever column names you want as long as they are unique! Press ENTER to skip a column!", sep = "\n")
            for (index in 1:ncol(split.temp)) {
                header.outer[length(header.outer) + 1] <- readline(prompt = paste0("What would be a proper column name for ", colnames(split.temp)[index], "? "))
            }
            cat(separator, sep = "\n")
        }
    }

    if (sum(test.1) != 0) {
        cat(paste0(sum(test.1), " SNPs have ':' in their IDs. Please be careful!"), sep = "\n")
        cat(separator, sep = "\n")
    }

    if (sum(test.2) != 0) {
        cat(paste0(sum(test.2), " SNPs have '_' in their IDs. Please be careful!"), sep = "\n")
        cat(separator, sep = "\n")
    }

    # Compare raw sumstats with the "list" (if provided)
    do.hapmap3 <- readline(prompt = paste0("Subset the dataset by HapMap3? Y/N "))
    cat(separator, sep = "\n")

    if (do.hapmap3 == "Y") {
        if (!file.exists("w_hm3.snplist.bz2")) {
            cat("I can't find HapMap3 list file and will try to download one from the Internet.", sep = "\n")
            cat(separator, sep = "\n")
            download.file("https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2", destfile = "w_hm3.snplist.bz2", quiet = TRUE)
        }
        list.hapmap3 <- as.data.frame(fread("w_hm3.snplist.bz2"))

        n.start <- nrow(raw)

        raw <- raw[raw$SNP %in% list.hapmap3$SNP, ]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
    }

    # Drop rows with MAF<0.01
    if ("MAF" %in% header.inner) {
        do.maf <- readline(prompt = paste0("Drop rows with MAF<=0.01? Y/N "))
        if (do.maf == "Y") {
        n.start <- nrow(raw)

        raw <- raw[raw$MAF >= 0.01,]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
        } else {
        cat(separator, sep = "\n")
        }
    }

    # Drop rows with INFO<0.9
    if ("INFO" %in% header.inner) {
        do.info <- readline(prompt = paste0("Drop rows with INFO <= 0.9? Y/N "))
        if (do.info == "Y") {
        n.start <- nrow(raw)

        raw <- raw[raw$INFO >= 0.9,]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
        } else {
        cat(separator, sep = "\n")
        }
    }

    # Drop rows not on chromosome 1-22
    if ("CHROMOSOME" %in% header.inner) {
        raw$CHROMOSOME <- as.character(raw$CHROMOSOME)
        if (length(union(unique(raw$CHROMOSOME), as.character(1:22))) != 22) {
            do.dropXY <- readline(prompt = paste0("Drop SNPs on chromosome X and Y? Y/N "))
            if (do.dropXY == "Y") {
                n.start <- nrow(raw)
                raw <- raw[raw$CHROMOSOME %in% as.character(1:22), ]
                n.end <- nrow(raw)
                cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
                cat(separator, sep = "\n")
            }
        }
    }

    # Drop rows for being on the MHC region
    if ("CHROMOSOME" %in% header.inner & "POSITION" %in% header.inner) {
        do.MHC <- readline(prompt = paste0("Drop rows on MHC region? Y/N "))
        if (do.MHC == "Y") {
            n.start <- nrow(raw)
            raw$CHROMOSOME <- as.character(raw$CHROMOSOME)
            raw <- raw[!(raw$CHROMOSOME == "6" & raw$POSITION >= 26000000 & raw$POSITION <= 34000000),]
            n.end <- nrow(raw)
            cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
            cat(separator, sep = "\n")
        } else {
            cat(separator, sep = "\n")
        }
    }

    # Transform P into usual range (0,1)
    if ("P" %in% header.inner) {
        if (max(raw$P) > 1 & min(raw$P) < 0) {
            do.transform <- readline(prompt = paste0("Transform P into usual range (0,1)? Y/N "))
            if (do.transform == "Y") {
                raw <- transform(raw, P = 10 ^ (-P))
                cat(separator, sep = "\n")
            } else {
                cat(separator, sep = "\n")
            }
        }
    }

    # Process A1 and A2
    if ("A1" %in% header.inner & "A2" %in% header.inner) {
        # MAKE A1 AND A2 CAPS
        raw$A1 <- toupper(raw$A1)
        raw$A2 <- toupper(raw$A2)

        # Erroneous A1/A2
        do.error <- readline(prompt = paste0("Drop rows with erroneous A1/A2 ? Y/N "))
        if (do.error == "Y") {
        n.start <- nrow(raw)

        raw <- raw[raw$A1 == 'A' | raw$A1 == 'C' | raw$A1 == 'T' | raw$A1 == 'G', ]
        raw <- raw[raw$A2 == 'A' | raw$A2 == 'C' | raw$A2 == 'T' | raw$A2 == 'G', ]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
        } else {
        cat(separator, sep = "\n")
        }

        # Ambigious SNPs
        do.ambi <- readline(prompt = paste0("Drop ambigious SNPs? Y/N "))
        if (do.ambi == "Y") {
        n.start <- nrow(raw)

        raw <- raw[!((raw$A1 == 'A' & raw$A2 == 'T') | (raw$A1 == 'T' & raw$A2 == 'A') | (raw$A1 == 'C' & raw$A2 == 'G') | (raw$A1 == 'G' & raw$A2 == 'C')),]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
        } else {
        cat(separator, sep = "\n")
        }
    } else {
        cat("WARNING: There might be some issue with your settings or the data itself! Because I didn't find any A1/A2 columns!", sep = "\n")
    }

    # Missing z-score?
    calculate.z <- FALSE
    if (!("Z" %in% header.inner)) {
        do.z <- readline(prompt = paste0("No z-score column is found. Do you want me to create one? Y/N "))
        if (do.z == "Y") {
            if ("BETA" %in% header.inner & "SE" %in% header.inner) {
                raw["Z"] <- raw$BETA / raw$SE
                calculate.z <- TRUE
            } else if ("ODDS_RATIO" %in% header.inner & "SE" %in% header.inner) {
                raw["Z"] <- log(raw$ODDS_RATIO) / raw$SE
                calculate.z <- TRUE
            } else if ("LOG_ODDS" %in% header.inner & "SE" %in% header.inner) {
                raw["Z"] <- raw$LOG_ODDS / raw$SE
                calculate.z <- TRUE
            } else if ("BETA" %in% header.inner & "P" %in% header.inner) {
                raw["Z"] <- sign(raw$BETA) * abs(qnorm(raw$P / 2))
                calculate.z <- TRUE
            } else if ("ODDS_RATIO" %in% header.inner & "P" %in% header.inner) {
                raw["Z"] <- sign(log(raw$ODDS_RATIO)) * abs(qnorm(raw$P / 2))
                calculate.z <- TRUE
            } else if ("LOG_ODDS" %in% header.inner & "P" %in% header.inner) {
                raw["Z"] <- sign(raw$ODDS_RATIO) * abs(qnorm(raw$P / 2))
                calculate.z <- TRUE
            } else {
                cat("I can't calculate z-score based on the information I have. SAD FACE EMOJI.", sep = "\n")
            }

            if (sum(is.na(raw$Z)) != 0) {
                n.start <- nrow(raw)
                raw <- raw[!is.na(raw$Z),]
                n.end <- nrow(raw)
                cat(paste0(n.start - n.end, " rows removed for having invalid z-score!"), sep = "\n")
            }
        }
        cat(separator, sep = "\n")
    }

    # Calculate chi^2 and sieve out observations with chi^2 larger than 80
    if ("P" %in% header.inner) {
        do.chi2 <- readline(prompt = paste0("Sieve out observations with chi^2 larger than 80? Y/N "))
        if (do.chi2 == "Y") {
        n.start <- nrow(raw)

        raw <- raw[qchisq(raw$P, 1, lower.tail = FALSE) <= 80,]

        n.end <- nrow(raw)
        cat(paste0(n.start - n.end, " rows removed!"), sep = "\n")
        cat(separator, sep = "\n")
        } else {
        cat(separator, sep = "\n")
        }
    }

    # Update the header
    if (calculate.z) {
        header.inner <- c(header.inner, "Z")
        header.outer <- c(header.outer, "Z")
    }
    colnames(raw) <- header.outer

    # Display output
    cat("Now take a final look at the output!", sep = "\n")
    print(head(raw))

    cat(separator, sep = "\n")

    # Output
    do.output <- readline(prompt = "Do you want to write the processed summary statistics to your hard drive? Y/N ") %>% toupper()
    if (do.output == "Y") {
        char.write <- readline(prompt = "Please input your desired filename (with extension and if no input is given, a default name will be assigned):")
        if (char.write == "") {
            char.write <- paste0(filename, ".sumstats")
        }
        fwrite(raw, file = char.write, sep = " ", row.names = FALSE, quote = FALSE)

        CHROMOSOME <- raw[, header.inner == "CHROMOSOME"]
        for (i in unique(CHROMOSOME)) {
            fwrite(raw[CHROMOSOME == i, ], file = paste0(filename, "-", i, ".sumstats"), sep = " ", row.names = FALSE, quote = FALSE)
        }
    }

    # Return?
    if (do.return) {
        return(raw)
    }
}
