# 02_population_structure.R
#
# population structure and kinship estimation using GAPIT.
#
# Inputs:
#   1. HapMap genotype file
#   2. GWAS metadata file (taxa + phenotype + country)
#
# Outputs:
#   PCA coordinates
#   PCA variance explained
#   PCA plots by geographic region
#   Kinship matrix
#   Kinship heatmap
#
# GAPIT handles SNP filtering (MAF >= 0.05).

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript 02_population_structure.R <genotype.hmp.txt> <gwas_metadata.csv>"
  )
}

genotype_file <- args[1]
metadata_file <- args[2]

# GAPIT installation
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

library(ggplot2)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# GAPIT only needs the taxa identifiers for this stage.
gwas_y <- metadata[, c("Taxa", "PlantHeight")]

# GAPIT's genotype module will calculate:
#   - PCA
#   - kinship
#   - filtered marker statistics
#
# We initially calculate 10 PCs so that we can inspect the
# population structure before deciding how many PCs to retain
# in the GWAS model.

g <- GAPIT(
  Y = gwas_y,
  G = read.table(genotype_file),
  PCA.total = 10,
  SNP.MAF = 0.05,
  SNP.test = FALSE,
  file.output = TRUE
)

# ---- Locate GAPIT outputs ----

if (!file.exists("GAPIT.Genotype.PCA.csv")) {
  stop("GAPIT PCA output was not generated.")
}

if (!file.exists("GAPIT.Genotype.PCA_eigenvalues.csv")) {
  stop("GAPIT PCA eigenvalue output was not generated.")
}

if (!file.exists("GAPIT.Genotype.Kin_Zhang.csv")) {
  stop("GAPIT kinship output was not generated.")
}

pca <- read.csv(
  "GAPIT.Genotype.PCA.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eigen <- read.csv(
  "GAPIT.Genotype.PCA_eigenvalues.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kin <- as.matrix(
  read.csv(
    "GAPIT.Genotype.Kin_Zhang.csv",
    row.names = 1,
    check.names = FALSE
  )
)

# ---- PCA / country matching ----

metadata$Taxa <- as.character(metadata$Taxa)
pca$taxa <- as.character(pca$taxa)

pca_meta <- merge(
  pca,
  metadata[, c("Taxa", "COUNTRY.of.origin")],
  by.x = "taxa",
  by.y = "Taxa",
  all.x = TRUE
)

# Broad geographic groups are used only to make PCA plots easier to interpret. They are not used as GWAS covariates.

pca_meta$Region <- "Other/Unknown"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("France", "Spain", "Portugal", "Italy", "Greece")
] <- "Western/Southern Europe"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("Poland", "Bulgaria", "Yugoslavia")
] <- "Eastern Europe"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("Iran", "Afghanistan", "India")
] <- "South/Central Asia"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("Turkey", "USSR")
] <- "Eurasia"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("Morocco", "Tunisia", "Canary Islands")
] <- "North Africa"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("China")
] <- "East Asia"

pca_meta$Region[
  pca_meta$COUNTRY.of.origin %in%
    c("Australia")
] <- "Oceania"

# ---- PCA variance explained ----

eig_numeric <- eigen[
  ,
  sapply(eigen, is.numeric),
  drop = FALSE
]

eig_values <- eig_numeric[[1]]

pca_variance <- data.frame(
  PC = seq_along(eig_values),
  Variance_explained = eig_values / sum(eig_values) * 100
)

pca_variance$Cumulative_variance <-
  cumsum(pca_variance$Variance_explained)

# ---- PCA figures ----

make_pca_plot <- function(x, y, filename, title) {
  
  p <- ggplot(
    pca_meta,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      colour = Region
    )
  ) +
    geom_point(size = 2, alpha = 0.75) +
    theme_classic() +
    labs(
      title = title,
      x = x,
      y = y,
      colour = "Geographic region"
    )
  
  ggsave(
    filename,
    p,
    width = 12,
    height = 8,
    dpi = 300
  )
}

make_pca_plot(
  "PC1", "PC2",
  "PCA_PC1_PC2_region.png",
  "Population structure: PC1 vs PC2"
)

make_pca_plot(
  "PC1", "PC3",
  "PCA_PC1_PC3_region.png",
  "Population structure: PC1 vs PC3"
)

make_pca_plot(
  "PC2", "PC3",
  "PCA_PC2_PC3_region.png",
  "Population structure: PC2 vs PC3"
)

# ---- PCA variance figure ----

p_var <- ggplot(
  pca_variance[1:min(15, nrow(pca_variance)), ],
  aes(x = PC, y = Variance_explained)
) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(
    breaks = 1:min(15, nrow(pca_variance))
  ) +
  theme_classic() +
  labs(
    title = "Variance explained by principal components",
    x = "Principal component",
    y = "Variance explained (%)"
  )

ggsave(
  "PCA_variance_explained.png",
  p_var,
  width = 10,
  height = 7,
  dpi = 300
)

# ---- Kinship heatmap ----
#
# This is a visualization of pairwise genetic relatedness as estimated by GAPIT. It will be used in the MLM.

png(
  "kinship_heatmap.png",
  width = 1800,
  height = 1800,
  res = 200
)

image(
  1:nrow(kin),
  1:ncol(kin),
  kin,
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Kinship matrix (Zhang method)"
)

box()

dev.off()

# ---- Save clean outputs ----

write.csv(
  pca_meta,
  "PCA_coordinates_with_region.csv",
  row.names = FALSE
)

write.csv(
  pca_variance,
  "PCA_variance_explained.csv",
  row.names = FALSE
)

write.csv(
  kin,
  "kinship_Zhang.csv"
)

cat("\n===== POPULATION STRUCTURE =====\n")
cat("Genotyped accessions:", nrow(pca), "\n")
cat("Number of PCs calculated:", length(eig_values), "\n")
cat("\nVariance explained:\n")
print(
  round(
    pca_variance[1:min(10, nrow(pca_variance)), ],
    3
  )
)
