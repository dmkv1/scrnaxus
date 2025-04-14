workflow QC {
    take:
    counts
    seed

    main:
    DROPLETS_TO_CELLS(counts, seed)

    combined_sce = DROPLETS_TO_CELLS.out.cells_sce.join(DROPLETS_TO_CELLS.out.droplets_sce)

    AMBIENT_RNA_REMOVAL(combined_sce, seed)
    
    DOUBLET_DETECTION(AMBIENT_RNA_REMOVAL.out.decont_sce, seed)

    emit:
    cells_sce = DROPLETS_TO_CELLS.out.cells_sce
    droplets_sce = DROPLETS_TO_CELLS.out.droplets_sce
}

process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/droplets_to_cells/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path(counts_dir)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: cells_sce
    tuple val(sample_id), path("${sample_id}_droplets.sce"), emit: droplets_sce
    path("FDRplot_${sample_id}.png")

    script:
    """
    droplets_to_cells.R ${seed} "${sample_id}" "${expected_cells}" "${patient_id}" "${timepoint}" "${compartment}" "${counts_dir}"
    """
}

process AMBIENT_RNA_REMOVAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/ambient_rna/${sample_id}/", mode: 'copy'

    container "quay.io/biocontainers/bioconductor-decontx:1.4.0--r44he5774e6_0"

    input:
    tuple val(sample_id), path(cells_sce), path(droplets_sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_decont.sce"), emit: decont_sce

    script:
    """
    ambient_RNA.R ${seed} "${sample_id}" "${cells_sce}" "${droplets_sce}"
    """
}

process DOUBLET_DETECTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/doublets/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), path(decont_sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_decont.sce"), emit: decont_sce

    script:
    """
    doublets.R ${seed} "${sample_id}" "${decont_sce}"
    """
}