# List of packages to install
packages_to_install <- c("data.table", "BEDMatrix", "dplyr", "MASS","optparse")

# Loop through each package and install if not already installed
for (pkg in packages_to_install) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
}

# Get the full path of the R script
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- script_path[grep("--file=", script_path)]
script_path <- sub("--file=", "", script_path)

# Extract the directory containing the R script
script_dir <- dirname(script_path)

# Set the working directory to the directory containing the R script
setwd(script_dir)



source("BLISSAssociation_Support.R")

option_list <- list(
make_option("--sumstats", type = "character", default = NA, action = "store", help = "The processed GWAS summary dataset file."
),
make_option("--sumstats_dir", type = "character", default = NA, action = "store", help = "The directory of GWAS summary dataset."
),
make_option("--N", type = "numeric", default = NA, action = "store", help = "The sample size of GWAS summary dataset. A numeric value."
),
make_option("--weights_models", type = "character", default = "UKB_EUR", action = "store", help = "The protein prediction models to be used."
),
make_option("--CHR", type = "numeric", default = NA, action = "store", help = "The chrosome to be runned. When unspecified, all chromsomes will be ran."
),
make_option("--output_dir", type = "character", default = getwd(), action = "store", help = "The directory of output file. Default is current directory."
),
make_option("--output_name", type = "character", default = "PWAS_results", action = "store", help = "The output file name."
)
)

opt <- parse_args(OptionParser(option_list = option_list))

sumstatfile     <- opt$sumstats
gwassum.dir        <- opt$sumstats_dir
n.sumstats        <- opt$N
weights.models <- opt$weights_models
CHR.list  <- opt$CHR
out.dir  <- opt$output_dir
output.name <- opt$output_name

CHR.list = as.numeric(CHR.list)
if(is.na(CHR.list)) {
    CHR.list = c(1:22)
}


#weights.models = "deCODE"
#CHR.list = c(22)
#sumstatfile = "Stroke_eur_GBMI_CHR22.sumstats"
#gwassum.dir = "/Users/cwu18/Dropbox/MDACC_research/Undergoing/PWAS-cis/BLISS-software/"
#out.dir = "/Users/cwu18/Dropbox/MDACC_research/Undergoing/PWAS-cis/BLISS-software/results/"

#create out directory
if (!dir.exists(out.dir)) {
    dir.create(out.dir)
    cat("Out directory created!\n")
} else {
    cat("Out directory already exists!\n")
}


if(weights.models=="UKB") {
    # read the weights
    load("models/UKB_weights_info.RData")
    ldref = "1000G/1000G.EUR.usedSNP.QC.CHR"
    weights.dir = "models/UKB/"
} else if (weights.models=="deCODE") {
    load("models/deCODE_weights_info.RData")
    ldref = "1000G/1000G.EUR.usedSNP.QC.CHR"
    weights.dir = "models/deCODE/"
} else if (weights.models == "ARIC") {
    load("models/ARIC_weights_info.RData")
    ldref = "1000G/1000G.EUR.usedSNP.QC.CHR"
    weights.dir = "models/ARIC/"
    
} else if (weights.models == "ARIC_AA") {
    load("models/ARIC_weights_info.RData")
    ldref = "1000G/1000G.AFR.usedSNP.QC.CHR"
    weights.dir = "models/ARIC_AA/"
    
} else if (weights.models == "UKB_AFR_std") {
    load("models/UKB_weights_info.RData")
    ldref = "1000G/1000G.AFR.usedSNP.QC.CHR"
    weights.dir = "models/UKB_AFR_std/"
    
} else if (weights.models == "UKB_AFR_super") {
    load("models/UKB_weights_info.RData")
    ldref = "1000G/1000G.AFR.usedSNP.QC.CHR"
    weights.dir = "models/UKB_AFR_super/"
    
} else if (weights.models == "UKB_ASN_std") {
    load("models/UKB_weights_info.RData")
    ldref = "1000G/1000G.EAS.usedSNP.QC.CHR"
    weights.dir = "models/UKB_ASN_std/"
} else if (weights.models == "UKB_ASN_super") {
    load("models/UKB_weights_info.RData")
    ldref = "1000G/1000G.EAS.usedSNP.QC.CHR"
    weights.dir = "models/UKB_ASN_super/"
} else {
    stop("Currently, we provide the following models: ARIC, ARIC_AA, deCODE, UKB, UKB_AFR_std, UKB_AFR_super, UKB_ASN_std, UKB_ASN_super. Please specify your choice for weights_models by selecting one of these available options." )
}

if (is.na(opt$sumstats) || is.na(opt$sumstats_dir)) {
    stop("Both --gwas and --gwassum_dir need to be specified.")
}

# load GWAS summary data
sumstats.org = fread(paste0(gwassum.dir,"/",sumstatfile)) %>% as.data.frame()

# Initialize the header
header.inner = colnames(sumstats.org)
header.inner <- tolower(header.inner)

# SNP
try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps","rsids")
header.inner[header.inner %in% try.snp] <- "SNP"

# Z-score
try.z <- c("zscore", "z-score", "gc_zscore", "z")
header.inner[header.inner %in% try.z] <- "Z"

# P
try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p")
header.inner[header.inner %in% try.p] <- "P"

# A1
try.a1 <- c("a1", "allele1", "allele_1", "nea", "non_effect_allele")
header.inner[header.inner %in% try.a1] <- "A1"

# A2
try.a2 <- c("a2", "allele2", "allele_2", "effect_allele","ea")
header.inner[header.inner %in% try.a2] <- "A2"

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

# Samplesize
try.samplesize <- c("n", "samplesize", "num_samples", "sample")
header.inner[header.inner %in% try.samplesize] <- "N"

# Update the header
colnames(sumstats.org) <- header.inner

# calculate Z score if does not exist
# Missing z-score?
#calculate.z <- FALSE
if (!("Z" %in% header.inner)) {
    do.z <- warining("No z-score column is found. We highly recommend using APSS.R function to pre-processing the data. We calculate the Z score based on the information we had.\n")
    
    if ("BETA" %in% header.inner & "SE" %in% header.inner) {
        sumstats.org["Z"] <- sumstats.org$BETA / sumstats.org$SE
        calculate.z <- TRUE
    } else if ("ODDS_RATIO" %in% header.inner & "SE" %in% header.inner) {
        sumstats.org["Z"] <- log(sumstats.org$ODDS_RATIO) / sumstats.org$SE
        calculate.z <- TRUE
    } else if ("LOG_ODDS" %in% header.inner & "SE" %in% header.inner) {
        sumstats.org["Z"] <- sumstats.org$LOG_ODDS / sumstats.org$SE
        calculate.z <- TRUE
    } else if ("BETA" %in% header.inner & "P" %in% header.inner) {
        sumstats.org["Z"] <- sign(sumstats.org$BETA) * abs(qnorm(sumstats.org$P / 2))
        calculate.z <- TRUE
    } else if ("ODDS_RATIO" %in% header.inner & "P" %in% header.inner) {
        sumstats.org["Z"] <- sign(log(sumstats.org$ODDS_RATIO)) * abs(qnorm(sumstats.org$P / 2))
        calculate.z <- TRUE
    } else if ("LOG_ODDS" %in% header.inner & "P" %in% header.inner) {
        sumstats.org["Z"] <- sign(sumstats.org$ODDS_RATIO) * abs(qnorm(sumstats.org$P / 2))
        calculate.z <- TRUE
    } else {
        cat("I can't calculate z-score based on the information I have. SAD FACE EMOJI.", sep = "\n")
    }
    
    if (sum(is.na(sumstats.org$Z)) != 0) {
        n.start <- nrow(sumstats.org)
        sumstats.org <- sumstats.org[!is.na(sumstats.org$Z),]
        n.end <- nrow(sumstats.org)
        cat(paste0(n.start - n.end, " rows removed for having invalid z-score!"), sep = "\n")
    }
    cat(separator, sep = "\n")
}




if(sum(grepl("chr",sumstats.org[,"CHROMOSOME"]))>100) {
    
} else {
    sumstats.org[,"CHROMOSOME"] = paste0("chr",sumstats.org[,"CHROMOSOME"])
}

sumstat.cols = colnames(sumstats.org)

# Check if necessary columns are present in the GWAS summary data
# Improved language for error messages and notes

if (!"CHROMOSOME" %in% sumstat.cols) {
    stop("The column 'CHROMOSOME' is missing from the GWAS summary data. Please reprocess the data to include this column.")
}

if (!"SNP" %in% sumstat.cols) {
    stop("The column 'SNP' is missing from the GWAS summary data. Please reprocess the data to include this column.")
}

if (!"A1" %in% sumstat.cols) {
    stop("The column 'A1' is missing from the GWAS summary data. Please reprocess the data to include this column.")
} else {
    cat("Note: In this data, 'A1' represents the non-effect allele.\n")
}

if (!"A2" %in% sumstat.cols) {
    stop("The column 'A2' is missing from the GWAS summary data. Please reprocess the data to include this column.")
} else {
    cat("Note: In this data, 'A2' represents the effect allele.\n")
}

if (!"N" %in% sumstat.cols) {
    cat("'N' is not present in the GWAS summary data.\n")
    
    if (is.na(n.sumstats)) {
        stop("The sample size ('N') is not provided. Please include it in the GWAS summary data or specify it as an argument.")
    }
}

if (!"Z" %in% sumstat.cols) {
    stop("The column 'Z' for Z-scores is missing from the GWAS summary data. Please reprocess the data to include this column.")
}

# obtain all proteins:
proteins = list.files(weights.dir)

outres = as.data.frame(matrix(NA,7000,12))

colnames(outres) = c("chr","p0","p1","gene","R2","Zscore.classic","p.classic",   "beta_BLISS","se_BLISS","p_BLISS","n_used_snp","n_snp")

ref = ldref
outindx = 1
for(chromosome in CHR.list) {
    reference.bim <- paste0(ref, chromosome, ".bim") %>% fread(., data.table = FALSE)
    reference.bed <- paste0(ref, chromosome) %>% BEDMatrix(., simple_names = TRUE)
    
    sumstats <- sumstats.org[sumstats.org$CHROMOSOME == paste0("chr", chromosome), ]
    sumstats = sumstats[rowSums(is.na(sumstats))==0, ]
    
    used.lookup = lookup[lookup[,"chr"] ==   chromosome,]
    used.proteins = proteins[proteins %in% paste0(used.lookup[,"gene"],".RData")]
    
    reference.dup <- reference.bim$V2 %>% duplicated()
    if (reference.dup %>% sum() != 0) {
        reference.bim <- reference.bim[!reference.dup, ]
        reference.bed <- reference.bed[, !reference.dup]
    }
    
    for(j in 1:length(used.proteins)) {
        load(paste0(weights.dir,used.proteins[j]))
        
        ref = ldref
        
        protein = gsub(".RData","",used.proteins[j])
        outres[outindx,1:4] = used.lookup[used.lookup[,"gene"]==protein,c("chr","p.0","p.1","gene")]
        outres[outindx,5] = R2
        outres[outindx,12] = sum(weight!=0)
        
        tryCatch({
            outres[outindx,6:11] = TestAssociation(sumstats, weight, ref, n.sumstats,reference.bim, reference.bed, skip.robust = TRUE)
        }, error=function(e){
            cat("Warning: No overlapping SNPs found in protein:", protein, "\n")
            cat("Error Details:", conditionMessage(e), "\n")
            
        })
        outindx = outindx + 1
        if(outindx %% 20 ==0 ) {
            cat("Finish Indx ",outindx,"\n")
            
        }
    }
    
    cat("Finish CHR",chromosome,"\n")
}

outres = outres[!is.na(outres[,1]),]

write.table(outres,quote = FALSE, row.names = FALSE, file = paste0(out.dir,"/",output.name))
