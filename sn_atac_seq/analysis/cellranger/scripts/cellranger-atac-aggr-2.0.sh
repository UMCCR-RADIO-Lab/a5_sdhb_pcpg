#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/pq08
#PBS -l mem=160GB
#PBS -l ncpus=48
#PBS -l wd

CELLR_ATAC=/g/data/pq08/software/cellranger-atac/cellranger-atac-2.0.0/cellranger-atac
REFERENCE=/g/data/pq08/reference/cellranger-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

AGGR_CSV=/g/data/pq08/projects/bowenb/a5_project/snATAC-seq/hg38/aggregate_atac_gadi.csv

echo "$CELLR_ATAC aggr --id=A5_snATAC_AGG_hg38 --reference=$REFERENCE --csv=$AGGR_CSV --nosecondary --localmem=160 --localcores=48"

$CELLR_ATAC aggr --id=A5_snATAC_AGG_hg38 --reference=$REFERENCE --csv=$AGGR_CSV --nosecondary --localmem=160 --localcores=48
