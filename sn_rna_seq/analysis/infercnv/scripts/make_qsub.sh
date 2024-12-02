
#!/bin/bash

qsub_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/analysis/infercnv/qsub"

threads=10
patients=(E018 E019-1 E123-1 E140-1 E143-1 E146-1 E156-1 E166-1 E171-1 E197-1 E225-1 P018-PGL1 P018-PGL3)

 
for patient_id in ${patients[@]}
do

qsub_file="${qsub_dir}/${patient_id}.qsub"

echo "
#PBS -N ${patient_id}
#PBS -l ncpus=${threads}
#PBS -l walltime=24:00:00
#PBS -l mem=256GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=50GB
#PBS -q hugemem
#PBS -m a
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/people/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/ppgl/a5/software/infer_cnv_conda

Rscript /g/data/pq08/projects/ppgl/a5/sn_rna_seq/analysis/infercnv/scripts/run_infercnv.r ${patient_id} ${threads}" >  "${qsub_file}"
done