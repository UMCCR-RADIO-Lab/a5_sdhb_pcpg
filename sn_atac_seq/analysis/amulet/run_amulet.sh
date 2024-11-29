#!/bin/bash
source "/g/data/pq08/people/flynna/software/miniconda3/etc/profile.d/conda.sh"
conda activate /g/data/pq08/projects/ppgl/a5/software/amulet

module load parallel
module load java/jdk-17.0.2 

samples=(E123-1 E156-1 E166-1 E188-1 E197-1 E200-1 E201-1 NAM018)

run_amulet () {
    base_dir="/g/data/pq08/projects/ppgl"
    counts_dir="${base_dir}/a5/sn_atac_seq/analysis/cellranger/counts_hg38"
    script_dir="${CONDA_PREFIX}/share/amulet"

    autosomes_list="${script_dir}/human_autosomes.txt"
    repeat_regions="${script_dir}/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed"

    sample=$1

    bam=${counts_dir}/${sample}/outs/possorted_bam.bam
    sc_csv=${counts_dir}/${sample}/outs/singlecell.csv
    out_dir="${base_dir}/a5/sn_atac_seq/analysis/amulet/${sample}"
    mkdir -p "${out_dir}" 

    log_stderr="${out_dir}/amulet.err"
    log_stdout="${out_dir}/amulet.out"

    amulet --forcesorted "${bam}" "${sc_csv}" "${autosomes_list}" "${repeat_regions}" "${out_dir}"  "${script_dir}" > "${log_stdout}" 2> "${log_stderr}" 
}

export -f run_amulet

parallel -j 8 run_amulet  ::: "${samples[@]}"