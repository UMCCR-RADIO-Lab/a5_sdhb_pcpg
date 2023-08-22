#!/usr/bin/bash

fa_ref="/g/data/pq08/local/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
base_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples"
qsub_dir="${base_dir}/qsub"
positions_dir="${base_dir}/positions"
pileup_dir="${base_dir}/pileups"

patients=(E122 E128 E132 E143 E146 E158 E159 E166 E167 E169 E225 E229)
for p in "${patients[@]}";
do 
    
    qsub_file="${qsub_dir}/${p}.qsub"

    echo "#!/bin/bash
    #PBS -N mpileup_${p}
    #PBS -l ncpus=1
    #PBS -l walltime=24:00:00
    #PBS -l mem=12GB
    #PBS -l storage=scratch/pq08+gdata/pq08
    #PBS -l jobfs=5GB
    #PBS -q normal
    #PBS -m a
    #PBS -M aidan.flynn@unimelb.edu.au

    module load parallel samtools/1.12
    bams=\$(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/${p}/final/${p}-*/ -iname '*-ready.*am' | sort | xargs)
    chr_pos_file=${positions_dir}/${p}.txt


    samtools mpileup -BQ0 -d10000000 --positions <(cut -f1,2 \${chr_pos_file}) --fasta-ref ${fa_ref} \${bams} > ${pileup_dir}/${p}.pileup
    " > "${qsub_file}"

done
