#!/bin/bash
set -e

nextflow run main.nf \
  -resume \
  --input samplesheet.csv \
  --genome_fasta /mnt/data/NGS/refs/10X_sc/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
  --gtf_file /mnt/data/NGS/refs/10X_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
  --read_length 150 \
  --cb_whitelist /mnt/data/NGS/refs/10X_sc/10X_barcodes/3M-february-2018.txt \
  --outdir results \
  --run_index
  