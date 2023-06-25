#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/pq08
#PBS -l mem=128GB
#PBS -l ncpus=32
#PBS -l wd

PROJECT_HOME=/g/data/pq08/projects/bowenb/a5_project/snATAC-seq
REFERENCE=/g/data/pq08/reference/cellranger-atac/refdata-cellranger-atac-hg19-1.2.0
CELLR_ATAC=/g/data/pq08/software/cellranger-atac/cellranger-atac-1.2.0/cellranger-atac
FASTQS=$PROJECT_HOME/atac_fqs
MEM=128
CORES=32

# samples E156-1,E156_topup=E156-1 E166-1,E166_topup=E166-1 NAM018 NPG-103 S126=E146-1 TMC-41=E197-1 NCCS322T=E200-1  NZ-WOLI4T=E201-1  TEX-558=E188-1  V-PH-58=E1230-1

$CELLR_ATAC count --id=E156-1 --reference=$REFERENCE --fastqs=$FASTQS/E156-1,$FASTQS/E156_topup --localmem=$MEM --localcores=$CORES --sample=E156-1,E156_topup
