 Reproducing GWAS pipeline using Watkin's wheat diversity panel and the plant height trait.

Data used are phenotypic data, HapMap SNP genotype and geographic metadata (country where seed was bought in 1930)

## Workflow

Phenotype + genotype
→ data/QC
→ population structure (PCA + kinship)
→ MLM GWAS
→ association results

## Methods

- Phenotype/genotype matching
- Phenotype quality control
- SNP filtering: MAF >= 0.05
- PCA for population structure
- Zhang kinship matrix
- Mixed Linear Model (MLM) GWAS
- Manhattan and QQ plots

## Software

- R
- GAPIT 3
- Nextflow
- Git
- Slurm-compatible configuration

## Reproducibility

The analysis is orchestrated with Nextflow, using R 4.2 and GAPIT 3.0. 
A further improvement of the reproducibility is expected to be carried out soon using a containerized env. (with Docker).

Local execution:

```bash
nextflow run main.nf -profile local