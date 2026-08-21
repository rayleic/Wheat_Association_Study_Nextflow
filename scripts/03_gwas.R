# ============================================================
# 03_gwas.R
#
# Genome-wide association study of PlantHeight using:
#   - RAW phenotype file
#   - RAW HapMap genotype file
#   - MLM
#   - 5 PCs
#   - kinship calculated internally by GAPIT
#   - SNP MAF >= 0.05
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript 03_gwas.R <phenotype.txt> <genotype.hmp.txt>"
  )
}

phenotype_file <- args[1]
genotype_file  <- args[2]

# ---- GAPIT initialization ----
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


# ---- RAW phenotype ----
#
# GAPIT receives the original phenotype file directly.
# We deliberately do not use the QC/metadata table here.

myY <- read.delim(phenotype_file,
  header = TRUE
)

myY <- myY[, c("Taxa", "PlantHeight")]

# ---- RAW HapMap ----

myG <- read.table(genotype_file,
  header = F
)

# ---- GWAS ----
#
# GAPIT performs the genotype/phenotype matching and SNP filtering.
# Kinship is calculated internally, as in the successful standalone run.

myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 5,
  SNP.MAF = 0.05,
  model = "MLM",
  file.output = TRUE
)

# ---- Compact results for the repository ----

gwas_files <- list.files(
  pattern = "^GAPIT.Association.GWAS_Results.*\\.csv$",
  full.names = TRUE
)

if (length(gwas_files) == 0) {
  stop("GWAS result file was not found.")
}

gwas <- read.csv(
  gwas_files[1],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  gwas,
  "GWAS_results.csv",
  row.names = FALSE
)

p_col <- grep(
  "P.value",
  names(gwas),
  ignore.case = TRUE,
  value = TRUE
)[1]

if (is.na(p_col)) {
  stop("Could not identify GWAS p-value column.")
}

gwas$P_value <- as.numeric(gwas[[p_col]])

top <- gwas[
  order(gwas$P_value),
  ,
  drop = FALSE
]

write.csv(
  head(top, 20),
  "top_20_SNPs.csv",
  row.names = FALSE
)

n_tested <- sum(!is.na(gwas$P_value))
bonferroni_threshold <- 0.05 / n_tested

significant <- gwas[
  !is.na(gwas$P_value) &
    gwas$P_value < bonferroni_threshold,
  ,
  drop = FALSE
]

write.csv(
  significant,
  "significant_SNPs_Bonferroni.csv",
  row.names = FALSE
)

write.table(
  data.frame(
    n_tested = n_tested,
    bonferroni_threshold = bonferroni_threshold,
    n_significant = nrow(significant),
    PCs_used = 5,
    MAF_threshold = 0.05,
    model = "MLM",
    kinship = "GAPIT internal"
  ),
  "GWAS_summary.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n===== GWAS SUMMARY =====\n")
cat("SNPs tested:", n_tested, "\n")
cat("Bonferroni threshold:", bonferroni_threshold, "\n")
cat("Significant SNPs:", nrow(significant), "\n")
cat("PCs used: 5\n")
cat("Model: MLM\n")
cat("Kinship: calculated internally by GAPIT\n")
cat("MAF threshold: 0.05\n")