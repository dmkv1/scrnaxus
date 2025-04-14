workflow QC {
    take:
    counts
    seed

    main:
    DROPLETS_TO_CELLS(counts, seed)

    DOUBLET_DETECTION(DROPLETS_TO_CELLS.out.sce, seed)

    CELL_QC(DOUBLET_DETECTION.out.sce, seed)

    emit:
    sce = CELL_QC.out.sce
}

process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/droplets_to_cells/", mode: 'copy'

    container "quay.io/biocontainers/bioconductor-dropletutils:1.26.0--r44h77050f0_1"

    input:
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path(counts_dir)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_droplet_metrics.json"), emit: metrics
    tuple val(sample_id), path("${sample_id}_droplets.sce")

    script:
    """
    droplets_to_cells.R ${seed} "${sample_id}" "${expected_cells}" "${patient_id}" "${timepoint}" "${compartment}" "${counts_dir}"
    """
}

// skipped for now
process AMBIENT_RNA_REMOVAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/ambient_rna/", mode: 'copy'

    container "quay.io/biocontainers/bioconductor-decontx:1.4.0--r44he5774e6_0"

    input:
    tuple val(sample_id), path(cells_sce), path(droplets_sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_decont.sce"), emit: sce

    script:
    """
    ambient_RNA.R ${seed} "${sample_id}" "${cells_sce}" "${droplets_sce}"
    """
}

process DOUBLET_DETECTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/doublets/", mode: 'copy'

    container "quay.io/biocontainers/bioconductor-scdblfinder:1.20.2--r44hdfd78af_0"

    input:
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_singlets.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_doublet_metrics.json"), emit: metrics

    script:
    """
    doublets.R ${seed} "${sample_id}" "${sce}"
    """
}

process CELL_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/cell_QC/", mode: 'copy'

    // container "docker://satijalab/seurat:5.0.0"

    input:
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_clean.sce"), emit: sce

    script:
    """
    cell_QC.R ${seed} "${sample_id}" "${sce}"
    """
}