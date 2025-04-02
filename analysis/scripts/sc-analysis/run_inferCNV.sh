#!/bin/bash
set -e

# set current dir to the script's dir
OUTDIR="/media/data/NGS/Projects/Clonal2/scRNAseq/analysis/outputs/inferCNV_trees_kclusters"
mkdir -p $OUTDIR
cd $OUTDIR

SCRIPT="/media/data/NGS/Projects/Clonal2/scRNAseq/analysis/scripts/latest/inferCNV.R"
cp $SCRIPT $OUTDIR

Rscript $SCRIPT $OUTDIR >> $OUTDIR/inferCNV.log 2>&1
