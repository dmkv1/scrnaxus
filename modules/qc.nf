workflow QC {
    take:
    counts
    seed

    main:
    DROPLETS_TO_CELLS(counts, seed)

    DOUBLET_DETECTION(
        file("assets/doublets.Rmd"),
        DROPLETS_TO_CELLS.out.sce,
        seed)

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

process DOUBLET_DETECTION {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/doublets/", mode: 'copy'

    // container "quay.io/biocontainers/bioconductor-scdblfinder:1.20.2--r44hdfd78af_0"

    input:
    path('doublets.Rmd')
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_singlets.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_doublet_metrics.json"), emit: metrics
    tuple val(sample_id), path("${sample_id}_doublets.nb.html"), emit: report

    script:
    """
    Rscript -e "rmarkdown::render('doublets.Rmd', \\
                output_file = '${sample_id}_doublets.nb.html', \\
                params = list(
                               sample_id = '${sample_id}', \\
                               path_sce_input = '${sce}', \\
                               path_sce_output = '${sample_id}_singlets.sce', \\
                               seed = ${seed}
                             ))"
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