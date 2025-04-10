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
    publishDir "${params.outdir}/alignment", mode: 'copy'

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
        --genomeLoad NoSharedMemory \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --outFilterType BySJout \\
        --limitBAMsortRAM 128000000000 \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \\
        --outFilterScoreMin 30 \\
        --outFilterMultimapNmax 10 \\
        --outSAMmultNmax 1 \\
        --soloMultiMappers EM \\
        --outMultimapperOrder Random \\
        --clipAdapterType CellRanger4 \\
        --soloType CB_UMI_Simple \\
        --soloStrand Forward \\
        --soloCBstart 1 \\
        --soloCBlen 16 \\
        --soloUMIstart 17 \\
        --soloUMIlen 10 \\
        --soloBarcodeReadLength 0 \\
        --soloUMIdedup 1MM_CR \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloCellFilter None \\
        --soloFeatures Gene GeneFull \\
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
