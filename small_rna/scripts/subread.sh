#!/bin/bash
#SBATCH --ntasks-per-node 6
#SBATCH --cpus-per-task 5
#SBATCH --nodes=1
#SBATCH -J subread
#SBATCH -p vccc
#SBATCH --mem=450G
#SBATCH -o %J.err
#SBATCH -e %J.out
#SBATCH --time=120:00:00
#SBATCH --mail-user=andrew.pattison@unimelb.edu.au
#SBATCH --mail-type=ALL

# Run subread on A5 small RNA-Seq

. /data/cephfs/punim0010/extras/Pattison/miniconda2/etc/profile.d/conda.sh
conda activate subread

run_subread(){

fastq=$1
sample=$(basename $fastq | cut -d "_" -f 1)
index=reference/Homo_sapiens.GRCh38.dna.primary_assembly 

# Trim adapters
cutadapt -a TGGAATTCTCGGGTGCCAAGG --minimum-length 23 -o trimmed_fastqs/$sample.trimmed.fastq $fastq

# Trim off first and last randomised 4 bases
cutadapt -u 4 -u -4 -o trimmed_fastqs/$sample.trimmed_2.fastq trimmed_fastqs/$sample.trimmed.fastq

# Run subread with 5 threads
subread-align -T 5 -i $index -t 0 -r trimmed_fastqs/$sample.trimmed_2.fastq -o bams/$sample.subread_results.bam

# Samtools sort and index results
samtools sort -@ 5 -o sorted_bams/$sample.subread_results.sorted.bam bams/$sample.subread_results.bam 

samtools index sorted_bams/$sample.subread_results.sorted.bam

}

export -f run_subread

mkdir -p trimmed_fastqs
mkdir -p bams
mkdir -p sorted_bams
mkdir -p counts
inputs=$(ls Merged_fastqs/*)

parallel --jobs 6 run_subread ::: $inputs

to_count=$(ls bams/*.bam)

#featureCounts -s 1 -T 5 -a reference/Homo_sapiens.GRCh38.97.chr.gtf.gz -o counts/A5_subread_counts.tsv $to_count

# Count with a mirbase ref 
featureCounts -F GFF -t miRNA -g Name -s 1 -T 5 -a reference/Mirbase_hsa_hg38.gff3 -o counts/A5_subread_counts_miRbase_ref.tsv $to_count


