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

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/A5/software/gistic2
module load parallel

ref_gene_file="$CONDA_PREFIX/share/gistic2*/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

base_dir="/g/data/pq08/projects/A5"
output_dir="${base_dir}/WGS/analysis/gistic/output"
input_dir="${base_dir}/WGS/analysis/gistic/input"

## input file definitions
declare -A seg_files
seg_files['all_samples']="${input_dir}/A5_purple_segs_gistic_format.seg"
seg_files['non_malignant_primary_nohn']="${input_dir}/A5_purple_segs_gistic_format_non_malignant_primary_nohn.seg"
seg_files['malignant_all']="${input_dir}/A5_purple_segs_gistic_format_malignant_all_nohn.seg"
seg_files['malignant_tert_atrx']="${input_dir}/A5_purple_segs_gistic_format_malignant_tert_atrx_nohn.seg"

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