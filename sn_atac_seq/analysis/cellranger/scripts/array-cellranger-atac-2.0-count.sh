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
FQS_DIR=/g/data/pq08/projects/bowenb/a5_project/snATAC-seq/atac_fqs
MEM=160
CORES=48

# NOTES: the -J is used to specify the number of jobs in the array, each job will have a different PBS_ARRAYID variable 
# the -r y specifies that the script is rerunnable 

# Get parameters from input.txt file using $PBS_ARRAY_INDEX as the line number 
# params=`sed -n "${PBS_ARRAY_INDEX} p" input.txt`

# store the line in an array, then pull out each parameter 
# paramsArray=($params)

# ID=${paramsArray[0]}
# FASTQS=${FQS_DIR}/${paramsArray[1]}
# SAMPLE=${paramsArray[2]}

# echo "$CELLR_ATAC count --id=$ID --reference=$REFERENCE --fastqs=$FQS_dir/E156-1,$FQS_dir/E156_topup  --localmem=$MEM --localcores=$CORES --sample=E156-1,E156"echo "$CELLR_ATAC count --id=$ID --reference=$REFERENCE --fastqs=$FQS_DIR/$FASTQS --localmem=$MEM --localcores=$CORES --sample=$SAMPLE" > ${ID}_${PBS_JOBID}.log

$CELLR_ATAC count --id=E166-1 --reference=$REFERENCE --fastqs=$FQS_DIR/E166-1,$FQS_DIR/E166_topup  --localmem=$MEM --localcores=$CORES --sample=E166-1,E166
