workflow QC {
    take:
    counts_ch

    main:
    CREATE_SCE(counts_ch)

    emit:
    sce_object = CREATE_SCE.out.sce_object
}

process CREATE_SCE {
    tag "${sample_id}"

    container 'https://depot.galaxyproject.org/singularity/bioconductor-dropletutils:1.22.0--r43hf17093f_0'

    input:
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path(counts_dir)

    output:
    tuple val(sample_id), path("${sample_id}_unfiltered.sce"), emit: sce_object

    script:
    """
    droplets_load.R "${sample_id}" "${expected_cells}" "${patient_id}" "${timepoint}" "${compartment}" "${counts_dir}"
    """
}
