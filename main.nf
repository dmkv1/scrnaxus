#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
// Alignment
include { GENOME_INDEX } from './modules/genome_index'
include { ALIGNMENT } from './modules/alignment'
// Droplet processing
include { DROPLETS_TO_CELLS } from './modules/qc'
include { DOUBLET_DETECTION } from './modules/qc'
include { CELL_QC } from './modules/qc'
// Cell annotations
include { SEURAT_CLUSTERING } from './modules/analysis'

def helpMessage() {
    log.info(
        """
    =========================================
    scRNA-seq Pipeline
    =========================================
    Usage:
    nextflow run main.nf --input samplesheet.csv --genome_fasta /path/to/genome.fa --gtf_file /path/to/gtf --cb_whitelist /path/to/barcodes --run_index

    Mandatory arguments:
      --input                Path to samplesheet CSV file
      --cb_whitelist         Path to cell barcode whitelist
      
    For genome indexing (required if --run_index is specified):
      --genome_fasta         Path to genome FASTA file
      --gtf_file             Path to GTF file for annotation
      --read_length          Read length for calculating sjdbOverhang (default: 150)
      --run_index            Flag to run genome indexing (default: false)
      
    If using pre-built index (when --run_index is not specified):
      --ref_dir              Path to pre-built STAR reference directory
      --gtf_file             Path to GTF file for annotation

    Optional arguments:
      --outdir               Output directory (default: ./results)
      --help                 Show this message and exit
    """.stripIndent()
    )
}

// Define the workflow
workflow {
    seed = Channel.value(params.seed)

    // Show help message if needed
    if (params.help) {
        helpMessage()
        exit(0)
    }

    // Validate input parameters
    if (params.input == null) {
        exit(1, "Input samplesheet not specified!")
    }

    // Ensure genome FASTA is provided when indexing is needed
    if (params.run_index && params.genome_fasta == null) {
        exit(1, "Genome FASTA file not specified for indexing!")
    }

    // Create channel from input samplesheet
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def expected_cells = row.expected_cells.toInteger()
            def patient_id = row.patient
            def timepoint = row.timepoint
            def compartment = row.compartment
            def fq_r1 = file(row.fq_r1)
            def fq_r2 = file(row.fq_r2)

            if (!fq_r1.exists()) {
                exit(1, "ERROR: Read 1 fastq file does not exist: ${fq_r1}")
            }
            if (!fq_r2.exists()) {
                exit(1, "ERROR: Read 2 fastq file does not exist: ${fq_r2}")
            }

            return [sample_id, expected_cells, patient_id, timepoint, compartment, fq_r1, fq_r2]
        }

    // Run genome indexing if requested
    if (params.run_index) {
        // Run the indexing step and use its output for alignment
        GENOME_INDEX(params.genome_fasta, params.gtf_file, params.read_length)

        // Execute alignment module with the newly created index
        ALIGNMENT(ch_input, GENOME_INDEX.out.index_dir, params.gtf_file, params.cb_whitelist)
    }
    else {
        // Use pre-built index
        if (params.ref_dir == null) {
            exit(1, "Pre-built STAR reference directory not specified!")
        }

        // Execute alignment module with pre-built index
        ALIGNMENT(ch_input, Channel.value(file(params.ref_dir)), params.gtf_file, params.cb_whitelist)
    }

    // Execute QC modules
    DROPLETS_TO_CELLS(
        file("templates/droplets_to_cells.Rmd"),
        ALIGNMENT.out.counts,
        seed,
    )

    DOUBLET_DETECTION(
        file("templates/doublets.Rmd"),
        DROPLETS_TO_CELLS.out.sce,
        seed,
    )

    CELL_QC(
        file("templates/cell_qc.Rmd"),
        DOUBLET_DETECTION.out.sce,
        seed
    )

    // Annotate cells
    SEURAT_CLUSTERING(
        file("templates/cell_types.Rmd"),
        CELL_QC.out.sce,
        seed
    )
}
