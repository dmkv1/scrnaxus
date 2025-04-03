workflow QC {
    take:
    counts_ch

    main:
    QC_PROCESS(counts_ch)

    emit:
    filtered_counts = QC_PROCESS.out.filtered_counts
}

process QC_PROCESS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(expected_cells), path(counts_dir)

    output:
    tuple val(sample_id), path("${sample_id}_filtered"), emit: filtered_counts

    script:
    """
    mkdir -p ${sample_id}_filtered
    cp -r ${counts_dir}/Gene/filtered/* ${sample_id}_filtered/
    echo "QC process completed for ${sample_id}" > ${sample_id}_filtered/qc_log.txt
    """
}
