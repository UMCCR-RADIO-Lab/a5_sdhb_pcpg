#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -J sinto
#SBATCH -p snowy
#SBATCH --mem=200G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=12:00:00
#SBATCH --mail-user=bowenb@student.unimelb.edu.au
#SBATCH --mail-type=ALL

module purge

. /data/gpfs/projects/punim0648/Projects/Blake/conda/miniconda3/etc/profile.d/conda.sh
conda activate sintoenv

PROJECT_HOME=/data/gpfs/projects/punim0648/Projects/Blake/A5_R_project
BAMS=$PROJECT_HOME/results/snRNA-seq/neo_only_bams/input_bams
CELLS_FILES=$PROJECT_HOME/results/snRNA-seq/neo_only_bams/cells_tsvs
CPUS=12


for i in $(echo $BAMS/*.bam)
do
	BAM=$i	
	CELLS_TSV=$(basename $i _possorted_genome_bam.bam)_neo_only_cells.tsv
	CELLS=$CELLS_FILES/$CELLS_TSV
	
	echo "sinto filterbarcodes -b $BAM -c $CELLS -p $CPUS"
	sinto filterbarcodes -b $BAM -c $CELLS -p $CPUS
done

