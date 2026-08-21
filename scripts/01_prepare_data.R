# 01_prepare_data.R
#
# Phenotype quick statistical analysis before using it for the GWAS.
#
# Genotype filtering is deliberately left to GAPIT.
# Country of origin is optional metadata and does not determine whether an accession is retained for GWAS.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {stop("Usage: Rscript 01_prepare_data.R <phenotype> <genotype> <country>")}

phenotype_file <- args[1]
genotype_file  <- args[2]
country_file   <- args[3]

output_metadata <- "gwas_metadata.csv"
output_qc       <- "qc_summary.csv"
output_summary  <- "plant_height_summary.csv"
output_plot     <- "plant_height_distribution.png"

# ------------------------------------------------------------
# Read phenotype and country data
# ------------------------------------------------------------

pheno <- read.delim(
  phenotype_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)

country <- read.csv(
  country_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Read only the HapMap header.
# This is enough to identify the accessions that have genotype data
# without loading the ~35,000-marker matrix at this stage.
geno_header <- read.delim(
  genotype_file,
  header = TRUE,
  nrows = 0,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

geno_taxa <- names(geno_header)[12:ncol(geno_header)]

# ------------------------------------------------------------
# Identify accessions usable for the GWAS
# ------------------------------------------------------------

# An accession is retained when:
#   - it is present in the phenotype data,
#   - it is present in the genotype data,
#   - PlantHeight is not missing.
# Country of origin is not required.

usable <- !is.na(pheno$PlantHeight) &
  pheno$Taxa %in% geno_taxa

gwas_metadata <- pheno[
  usable,
  c("Taxa", "PlantHeight")
]

# Add country where available.
# Missing countries remain NA and do not lead to exclusion.
gwas_metadata <- merge(
  gwas_metadata,
  country[, c("acc_no", "COUNTRY.of.origin")],
  by.x = "Taxa",
  by.y = "acc_no",
  all.x = TRUE,
  sort = FALSE
)

# Preserve the original phenotype order after merging.
gwas_metadata <- gwas_metadata[
  match(
    pheno$Taxa[usable],
    gwas_metadata$Taxa
  ),
]

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

phenotype_taxa <- unique(pheno$Taxa)
genotype_taxa  <- unique(geno_taxa)

common_taxa <- intersect(
  phenotype_taxa,
  genotype_taxa
)

common_with_phenotype <- pheno$Taxa[
  pheno$Taxa %in% genotype_taxa &
    !is.na(pheno$PlantHeight)
]

qc_summary <- data.frame(
  metric = c(
    "Phenotype accessions",
    "Genotype accessions",
    "Accessions with phenotype and genotype",
    "Common accessions with non-missing PlantHeight",
    "Country missing among retained accessions"
  ),
  value = c(
    length(phenotype_taxa),
    length(genotype_taxa),
    length(common_taxa),
    length(unique(common_with_phenotype)),
    sum(is.na(gwas_metadata$COUNTRY.of.origin))
  )
)

# ------------------------------------------------------------
# Plant height summary statistics
# ------------------------------------------------------------

plant_height_summary <- data.frame(
  n = length(gwas_metadata$PlantHeight),
  mean = mean(gwas_metadata$PlantHeight),
  sd = sd(gwas_metadata$PlantHeight),
  min = min(gwas_metadata$PlantHeight),
  median = median(gwas_metadata$PlantHeight),
  max = max(gwas_metadata$PlantHeight)
)

# ------------------------------------------------------------
# Phenotype distribution
# ------------------------------------------------------------

png(
  output_plot,
  width = 1800,
  height = 1200,
  res = 300
)

hist(
  gwas_metadata$PlantHeight,
  breaks = 30,
  main = "Plant height distribution",
  xlab = "Plant height",
  ylab = "Number of accessions"
)

dev.off()

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

write.csv(
  gwas_metadata,
  output_metadata,
  row.names = FALSE
)

write.csv(
  qc_summary,
  output_qc,
  row.names = FALSE
)

write.csv(
  plant_height_summary,
  output_summary,
  row.names = FALSE
)

# Console summary

cat("\n===== DATA PREPARATION =====\n")
cat("Phenotype accessions:", length(phenotype_taxa), "\n")
cat("Genotype accessions:", length(genotype_taxa), "\n")
cat("Common phenotype/genotype accessions:", length(common_taxa), "\n")
cat("Accessions retained for GWAS:", nrow(gwas_metadata), "\n")
cat(
  "Retained accessions without country:",
  sum(is.na(gwas_metadata$COUNTRY.of.origin)),
  "\n"
)
