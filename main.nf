nextflow.enable.dsl=2

process PREPARE_DATA {

    tag "prepare phenotype/genotype data"

    publishDir "${projectDir}/results/phenotype", mode: 'copy'

    input:
    path phenotype
    path genotype
    path country

    output:
    path "gwas_metadata.csv", emit: metadata
    path "qc_summary.csv", emit: qc
    path "plant_height_summary.csv", emit: summary
    path "plant_height_distribution.png", emit: phenotype_plot

    script:
    """
    Rscript ${projectDir}/scripts/01_prepare_data.R \
        ${phenotype} \
        ${genotype} \
        ${country}
    """
}


process POPULATION_STRUCTURE {

    tag "PCA + kinship"

    publishDir "${projectDir}/results/population_structure", mode: 'copy'

    input:
    path genotype
    path metadata

    output:
    path "PCA_coordinates_with_region.csv", emit: pca
    path "PCA_variance_explained.csv", emit: eigen
    path "PCA_PC1_PC2_region.png"
    path "PCA_PC1_PC3_region.png"
    path "PCA_PC2_PC3_region.png"
    path "PCA_variance_explained.png"
    path "kinship_Zhang.csv", emit: kinship
    path "kinship_heatmap.png"

    script:
    """
    Rscript ${projectDir}/scripts/02_population_structure.R \
        ${genotype} \
        ${metadata}
    """
}


process GWAS {

    tag "PlantHeight MLM GWAS"

    publishDir "${projectDir}/results/gwas", mode: 'copy'

    input:
    path genotype
    path phenotype

    output:
    path "GWAS_results.csv"
    path "top_20_SNPs.csv"
    path "significant_SNPs_Bonferroni.csv"
    path "GWAS_summary.txt"
    path "GAPIT.Association.QQ.*.pdf", optional: true
    path "GAPIT.Association.Manhattan.*.pdf", optional: true

    script:
    """
    Rscript ${projectDir}/scripts/03_gwas.R \
        ${phenotype} \
        ${genotype}
    """
}


workflow {

    phenotype = file("${projectDir}/data/phenotype.txt")
    genotype  = file("${projectDir}/data/genotype.hmp.txt")
    country   = file("${projectDir}/data/country.csv")

    prepared = PREPARE_DATA(
        phenotype,
        genotype,
        country
    )

    structure = POPULATION_STRUCTURE(
        genotype,
        prepared.metadata
    )

    GWAS(
        genotype,
        phenotype
    )
}