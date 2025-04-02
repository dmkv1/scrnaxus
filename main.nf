#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { ALIGNMENT } from './modules/alignment'
include { QC }        from './modules/qc'
include { ANALYSIS }  from './modules/analysis'

// Define the workflow
workflow {
    // Show help message if needed
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Validate input parameters
    if (params.input == null) {
        exit 1, "Input samplesheet not specified!"
    }
    
    // Create channel from input samplesheet
    ch_input = Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> 
            def sample_id = row.sample
            def expected_cells = row.expected_cells.toInteger()
            def fq_r1 = file(row.fq_r1)
            def fq_r2 = file(row.fq_r2)
            
            if (!fq_r1.exists()) exit 1, "ERROR: Read 1 fastq file does not exist: ${fq_r1}"
            if (!fq_r2.exists()) exit 1, "ERROR: Read 2 fastq file does not exist: ${fq_r2}"
            
            return [ sample_id, expected_cells, fq_r1, fq_r2 ]
        }
    
    // Execute alignment module
    ALIGNMENT(ch_input, params.ref_dir, params.gtf_file, params.cb_whitelist)
    
    // Execute QC module
    QC(ALIGNMENT.out.counts)
    
    // Execute analysis module
    ANALYSIS(QC.out.filtered_counts)
}

def helpMessage() {
    log.info"""
    =========================================
    scRNA-seq Pipeline
    =========================================
    Usage:
    nextflow run main.nf --input samplesheet.csv --ref_dir /path/to/reference --gtf_file /path/to/gtf --cb_whitelist /path/to/barcodes

    Mandatory arguments:
      --input                Path to samplesheet CSV file
      --ref_dir              Path to STAR reference directory
      --gtf_file             Path to GTF file for annotation
      --cb_whitelist         Path to cell barcode whitelist

    Optional arguments:
      --outdir               Output directory (default: ./results)
      --help                 Show this message and exit
    """.stripIndent()
}