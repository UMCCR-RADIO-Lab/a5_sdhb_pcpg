#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/pq08
#PBS -l mem=160GB
#PBS -l ncpus=48
#PBS -l wd

CELLR_ATAC=/g/data/pq08/software/cellranger-atac/cellranger-atac-2.0.0/cellranger-atac
REFERENCE=/g/data/pq08/reference/cellranger-atac/refdata-cellranger-atac-hg19-2.0.0/ensembl/hg19
FASTQS=/g/data/pq08/projects/bowenb/a5_project/snATAC-seq/atac_fqs
MEM=160
CORES=48

# samples: E156-1,E156_topup=E156-1 E166-1,E166_topup=E166-1 NAM018 NPG-103 S126=E146-1 TMC-41=E197-1 NCCS322T=E200-1 NZ-WOLI4T=E201-1 TEX-558=E188-1 V-PH-58=E123-1 E140-1 E171-1 E225-1

# --id is the name to call the output file 
# --fastqs is the filepath of the folder(s) with the fastqs inside (comma-separated list used to specify multiple folders
# --sample is the name(s) prepended to the fastqs, comma-sep list is ok

$CELLR_ATAC count --id=E225-1 --reference=$REFERENCE --fastqs=$FASTQS/E225-1 --localmem=$MEM --localcores=$CORES --sample=E225
