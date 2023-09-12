#!/bin/bash
#PBS -N A5_gistic
#PBS -l ncpus=3
#PBS -l walltime=3:00:00
#PBS -l mem=24gb
#PBS -l storage=scratch/pq08+gdata/pq08 
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

threads=3

## run GISTIC analysis on GADI (modified from Andrew's Spartan scripts)

conda activate /g/data/pq08/projects/ppgl/a5/software/gistic2
module load parallel

ref_gene_file="$CONDA_PREFIX/share/gistic2*/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

base_dir="/g/data/pq08/projects/ppgl/a5"
gistic_dir="${base_dir}/wgs/analysis/gistic"
input_dir="${gistic_dir}/input"
output_dir="${gistic_dir}/output"


## input file definitions
declare -A seg_files
seg_files['all']="${input_dir}/a5_purple_segs_gistic_format.seg"
seg_files['non_metastatic_atpgl']="${input_dir}/a5_purple_segs_gistic_format_non_metastatic_atpgl.seg"
seg_files['metastatic_atpgl']="${input_dir}/a5_purple_segs_gistic_format_metastatic_atpgl.seg"
seg_files['tert']="${input_dir}/a5_purple_segs_gistic_format_tert_atpgl.seg"
seg_files['atrx']="${input_dir}/a5_purple_segs_gistic_format_atrx_atpgl.seg"
seg_files['hn']="${input_dir}/a5_purple_segs_gistic_format_all_hnpgl.seg"

#runner function for parallel
run_gistic() { #outdir, seg_group, seg_file
    output_dir=$1 
    seg_group=$2
    seg_file=$3
    ref_gene_file=$4
  
    echo "Processing ${seg_group}: ${seg_file} ..."
    echo "--- creating output directory ---"
    mkdir -p "${output_dir}/${seg_group}"
    
    echo "--- running GISTIC ---"
    gistic2 -b "${output_dir}/${seg_group}" \
    -seg "${seg_file}" \
    -refgene "${ref_gene_file}" \
    -genegistic 1 \
    -smallmem 0 \
    -broad 1 \
    -brlen 0.5 \
    -conf 0.90 \
    -armpeel 1 \
    -savegene 1 \
    -gcm extreme \
    -rx 0 > "${output_dir}/${seg_group}/${seg_group}_gistic.log"
    
}

export -f run_gistic

parallel --link -j $threads "run_gistic $output_dir {1} {2} $ref_gene_file" ::: "${!seg_files[@]}" ::: "${seg_files[@]}"