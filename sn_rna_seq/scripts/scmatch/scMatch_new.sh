#!/bin/bash
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --nodes=1
#SBATCH -J scmatch
#SBATCH -p vccc
#SBATCH --mem=400G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=120:00:00
#SBATCH --mail-user=bowenb@student.unimelb.edu.au
#SBATCH --mail-type=ALL

. /data/gpfs/projects/punim0648/Projects/Blake/conda/miniconda3/etc/profile.d/conda.sh
conda activate scMatch

RAW_COUNTS_MAT=/data/gpfs/projects/punim0648/Projects/Blake/scATAC_R_Project/Data/scMatch/input/rna_all_rawcounts_mat.csv
SCMATCH=/data/gpfs/projects/punim0648/Projects/Blake/scMatch_files/scMatch

python $SCMATCH/scMatch.py --refDS $SCMATCH/refDB/FANTOM5 --dFormat csv --testDS $RAW_COUNTS_MAT --coreNum 16
