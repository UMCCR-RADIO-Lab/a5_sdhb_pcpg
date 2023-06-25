#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/pq08
#PBS -l mem=128GB
#PBS -l ncpus=32
#PBS -l wd

REFERENCE=/g/data/pq08/reference/cellranger-atac/refdata-cellranger-atac-hg19-1.2.0
CellR_ATAC=/g/data/pq08/software/cellranger-atac/cellranger-atac-1.2.0/cellranger-atac

AGGR_CSV=/g/data/pq08/projects/bowenb/a5_project/snATAC-seq/cellranger_outs/aggregated/aggregate_atac_gadi.csv

# aggregates the samples excluding NPG-103 
$CellR_ATAC aggr --id=A5_snATAC_AGG_02 --reference=$REFERENCE --csv=$AGGR_CSV --nosecondary --localmem=128 --localcores=32


