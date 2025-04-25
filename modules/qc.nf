process DROPLETS_TO_CELLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/droplets_to_cells/", mode: 'copy'

    // container "quay.io/biocontainers/bioconductor-dropletutils:1.26.0--r44h77050f0_1"

    input:
    path('droplets_to_cells.Rmd')
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path(counts_dir)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_droplets_to_cells.nb.html"), emit: report
    tuple val(sample_id), path("${sample_id}_cells.sce"), emit: sce
    tuple val(sample_id), path("${sample_id}_droplets.sce")
    tuple val(sample_id), path("${sample_id}_droplet_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('droplets_to_cells.Rmd',
                output_file = '${sample_id}_droplets_to_cells.nb.html',
                params = list(
                               sample_id = '${sample_id}',
                               expected_cells = '${expected_cells}',
                               patient_id = '${patient_id}',
                               timepoint  = '${timepoint}',
                               compartment = '${compartment}',
                               counts_dir = '${counts_dir}',
                               FDR_thresh = '${params.qc.FDR_thresh}',
                               path_sce_output = '${sample_id}_cells.sce',
                               seed = ${seed}
                             ))"
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
    Rscript -e "rmarkdown::render('doublets.Rmd',
                output_file = '${sample_id}_doublets.nb.html',
                params = list(
                               sample_id = '${sample_id}',
                               path_sce_input = '${sce}',
                               path_sce_output = '${sample_id}_singlets.sce',
                               seed = ${seed}
                             ))"
    """
}

process CELL_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/cell_QC/", mode: 'copy'

    // container "docker://satijalab/seurat:5.0.0"

    input:
    path('cell_qc.Rmd')
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_cell_QC.nb.html"), emit: report
    // tuple val(sample_id), path("${sample_id}_clean.sce"), emit: sce

    script:
    """
    Rscript -e "rmarkdown::render('cell_qc.Rmd',
                output_file = '${sample_id}_cell_QC.nb.html',
                params = list(
                    sample_id = '${sample_id}',
                    path_sce_input = '${sce}',
                    path_sce_output = '${sample_id}_cells.sce',
                    nUMI_thresh = '${params.qc.nUMI_thresh}',
                    nGenes_thresh = '${params.qc.nGenes_thresh}',
                    mitochondrial_thresh = '${params.qc.mitochondrial_thresh}',
                    seed = ${seed}
                    ))"
    """
}