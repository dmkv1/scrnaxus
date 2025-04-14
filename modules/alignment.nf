workflow ALIGNMENT {
    take:
    ch_input // Channel: [ sample_id, expected_cells, fq_r1, fq_r2 ]
    ref_dir // Path: STAR reference directory
    gtf_file // Path: GTF file
    cb_whitelist // Path: Cell barcode whitelist

    main:
    STAR_ALIGNMENT(ch_input, ref_dir, gtf_file, cb_whitelist)

    emit:
    counts = STAR_ALIGNMENT.out.counts
    versions = STAR_ALIGNMENT.out.versions
}
process STAR_ALIGNMENT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    input:
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path(fq_r1), path(fq_r2)
    path ref_dir
    path gtf_file
    path cb_whitelist

    output:
    tuple val(sample_id), val(expected_cells), val(patient_id), val(timepoint), val(compartment), path("${sample_id}_Solo.out"), emit: counts
    path "versions.yml", emit: versions

    script:
    def gtf_filename = gtf_file.getName()
    def unzip_gtf = gtf_filename.endsWith('.gz') ? "gunzip -c ${gtf_file} > uncompressed.gtf && GTF_FILE=uncompressed.gtf" : "GTF_FILE=${gtf_file}"

    """
    # Decompress GTF file if needed
    ${unzip_gtf}
    
    STAR \\
        --genomeLoad ${params.starsolo.genomeLoad} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand ${params.starsolo.readFilesCommand} \\
        --outFilterType ${params.starsolo.outFilterType} \\
        --limitBAMsortRAM ${params.starsolo.limitBAMsortRAM} \\
        --outSAMtype ${params.starsolo.outSAMtype} \\
        --outSAMunmapped ${params.starsolo.outSAMunmapped} \\
        --outSAMattributes ${params.starsolo.outSAMattributes} \\
        --outFilterScoreMin ${params.starsolo.outFilterScoreMin} \\
        --outFilterMultimapNmax ${params.starsolo.outFilterMultimapNmax} \\
        --outSAMmultNmax ${params.starsolo.outSAMmultNmax} \\
        --soloMultiMappers ${params.starsolo.soloMultiMappers} \\
        --outMultimapperOrder ${params.starsolo.outMultimapperOrder} \\
        --clipAdapterType ${params.starsolo.clipAdapterType} \\
        --soloType ${params.starsolo.soloType} \\
        --soloStrand ${params.starsolo.soloStrand} \\
        --soloCBstart ${params.starsolo.soloCBstart} \\
        --soloCBlen ${params.starsolo.soloCBlen} \\
        --soloUMIstart ${params.starsolo.soloUMIstart} \\
        --soloUMIlen ${params.starsolo.soloUMIlen} \\
        --soloBarcodeReadLength ${params.starsolo.soloBarcodeReadLength} \\
        --soloUMIdedup ${params.starsolo.soloUMIdedup} \\
        --soloUMIfiltering ${params.starsolo.soloUMIfiltering} \\
        --soloCBmatchWLtype ${params.starsolo.soloCBmatchWLtype} \\
        --soloCellFilter ${params.starsolo.soloCellFilter} \\
        --soloFeatures ${params.starsolo.soloFeatures} \\
        --genomeDir ${ref_dir} \\
        --sjdbGTFfile \$GTF_FILE \\
        --soloCBwhitelist ${cb_whitelist} \\
        --readFilesIn ${fq_r2} ${fq_r1} \\
        --outFileNamePrefix "${sample_id}_"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
