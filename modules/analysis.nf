// modules/analysis.nf

process ANALYSIS_PROCESS {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(filtered_counts)
    
    output:
    path "${sample_id}_results"
    
    script:
    """
    mkdir -p ${sample_id}_results
    echo "Sample ID: ${sample_id}" > ${sample_id}_results/summary.txt
    echo "Analysis completed successfully" >> ${sample_id}_results/summary.txt
    """
}

workflow ANALYSIS {
    take:
    filtered_counts_ch
    
    main:
    ANALYSIS_PROCESS(filtered_counts_ch)
}