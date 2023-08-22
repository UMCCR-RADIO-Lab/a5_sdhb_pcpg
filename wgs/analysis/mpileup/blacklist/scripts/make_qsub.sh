#!/usr/bin/bash

fa_ref="/g/data/pq08/local/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa"
base_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup"
qsub_dir=${base_dir}/qsub
positions_dir=${base_dir}/positions
pileup_dir=${base_dir}/pileups

threads=7

chroms=($(seq 1 1 22) X Y)
for c in "${chroms[@]}";
do 
    chr="chr${c}"
    qsub_file="${qsub_dir}/${chr}.qsub"

    echo "#!/bin/bash
    #PBS -N mpileup_${chr}
    #PBS -l ncpus=${threads}
    #PBS -l walltime=24:00:00
    #PBS -l mem=32GB
    #PBS -l storage=scratch/pq08+gdata/pq08
    #PBS -l jobfs=5GB
    #PBS -q normal
    #PBS -m a
    #PBS -M aidan.flynn@unimelb.edu.au

    module load parallel samtools/1.12

    bams=\$(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E*/final/E*-B01/ -iname '*-B01-ready.bam' | xargs)
    mapfile chr_pos_files < <(find ${positions_dir} -iname \"${chr}_*\")


    parallel -j ${threads} \"samtools mpileup -BQ0 -d10000000 -r ${chr} --positions <(cut -f1,2 {}) --fasta-ref ${fa_ref} \${bams} > ${pileup_dir}/{/.}.pileup\" ::: \"\${chr_pos_files[@]}\"
    " > "${qsub_file}"

done
