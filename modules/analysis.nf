process SEURAT_CLUSTERING {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/annotations/", mode: 'copy'

    input:
    path('seurat_clustering.Rmd')
    tuple val(sample_id), path(sce)
    val seed

    output:
    tuple val(sample_id), path("${sample_id}_seurat_clustering.nb.html"), emit: report
    // tuple val(sample_id), path("${sample_id}_annotated.sce"), emit: sce
    // tuple val(sample_id), path("${sample_id}_cell_types_metrics.json"), emit: metrics

    script:
    """
    Rscript -e "rmarkdown::render('seurat_clustering.Rmd',
                output_file = '${sample_id}_seurat_clustering.nb.html',
                params = list(
                    sample_id = '${sample_id}',
                    path_sce_input = '${sce}',
                    path_sce_output = '${sample_id}.sce',
                    seed = ${seed}
                    ))"
    """
}