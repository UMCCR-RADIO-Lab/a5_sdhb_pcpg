#!/bin/bash
#PBS -N pairtree
#PBS -l ncpus=12
#PBS -l walltime=12:00:00
#PBS -l mem=90GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=50GB
#PBS -q normalbw
#PBS -m a
#PBS -M aidan.flynn@unimelb.edu.au

pairtree_dir="/g/data/pq08/projects/ppgl/a5/software/pairtree"
pairtree_utils_dir="${pairtree_dir}/util"
pairtree_bin_dir="${pairtree_dir}/bin"

source "/g/data/pq08/people/flynna/software/miniconda3/etc/profile.d/conda.sh"
conda activate ${pairtree_dir}/conda


input_dir_base="/g/data/pq08/projects/ppgl/a5/wgs/analysis/pairtree/input"

patients=(E122 E128 E132 E136 E143 E146 E158 E159 E167 E169 E225 E229)
patients=(E225)
for patient_id in ${patients[@]}
do
python3 ${pairtree_utils_dir}/fix_bad_var_read_prob.py \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json" \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json" \
--action add_to_garbage \
--verbose &
done
wait

# threads=14
# for patient_id in ${patients[@]}
# do
# echo "
# #PBS -N ${patient_id}
# #PBS -l ncpus=${threads}
# #PBS -l walltime=12:00:00
# #PBS -l mem=510GB
# #PBS -l storage=scratch/pq08+gdata/pq08
# #PBS -l jobfs=50GB
# #PBS -q hugemem
# #PBS -m a
# #PBS -M aidan.flynn@unimelb.edu.au

# source /g/data/pq08/people/flynna/software/miniconda3/etc/profile.d/conda.sh
# conda activate /g/data/pq08/projects/ppgl/a5/software/pairtree/conda

# python3 ${pairtree_bin_dir}/removegarbage \
# --parallel ${threads} \
# --verbose \
# ${input_dir_base}/${patient_id}/${patient_id}.ssm \
# ${input_dir_base}/${patient_id}/${patient_id}.params.json \
# ${input_dir_base}/${patient_id}/${patient_id}.params.json" | qsub

# done


for patient_id in ${patients[@]}
do
python3 ${pairtree_bin_dir}/clustervars \
--parallel 12 \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json"
done


for patient_id in ${patients[@]}
do
python3 ${pairtree_bin_dir}/pairtree \
--parallel 12 \
--seed=5555 \
--params "${input_dir_base}/${patient_id}/${patient_id}.params.json" \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.npz"
done

for patient_id in ${patients[@]}
do
python3 ${pairtree_bin_dir}/summposterior \
--runid "${patient_id}" \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json" \
"${input_dir_base}/${patient_id}/${patient_id}.npz" \
"${input_dir_base}/${patient_id}/${patient_id}.summposterior.html" &
done
wait

for patient_id in ${patients[@]}
do
python3 ${pairtree_bin_dir}/plottree \
--runid ${patient_id} \
"${input_dir_base}/${patient_id}/${patient_id}.ssm" \
"${input_dir_base}/${patient_id}/${patient_id}.params.json" \
"${input_dir_base}/${patient_id}/${patient_id}.npz" \
"${input_dir_base}/${patient_id}/${patient_id}.plottree.html" &
done
wait