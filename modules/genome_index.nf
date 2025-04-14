workflow GENOME_INDEX {
    take:
    genome_fasta // Path: Genome FASTA file
    gtf_file // Path: GTF file for annotation
    read_length // Value: Read length for sjdbOverhang calculation

    main:
    STAR_INDEX(genome_fasta, gtf_file, read_length)

    emit:
    index_dir = STAR_INDEX.out.index_dir
    versions = STAR_INDEX.out.versions
}
process STAR_INDEX {
    tag "Genome indexing"
    publishDir "${params.outdir}/genome_index", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.11a--h0033a41_0'

    input:
    path genome_fasta
    path gtf_file
    val read_length

    output:
    path "star_index", emit: index_dir
    path "versions.yml", emit: versions

    script:
    def sjdbOverhang = read_length - 1
    def gtf_filename = gtf_file.getName()
    def unzip_gtf = gtf_filename.endsWith('.gz') ? "gunzip -c ${gtf_file} > uncompressed.gtf && GTF_FILE=uncompressed.gtf" : "GTF_FILE=${gtf_file}"

    """
    mkdir -p star_index
    
    # Decompress GTF file if needed
    ${unzip_gtf}
    
    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile \$GTF_FILE \\
        --sjdbOverhang ${sjdbOverhang}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
