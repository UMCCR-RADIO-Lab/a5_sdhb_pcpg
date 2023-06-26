#!/bin/bash


base_dir="/g/data/pq08/projects/ppgl/a5/wgs"
qsub_dir="${base_dir}/analysis/bcbio/qsub"
fastq_dir="${base_dir}/fastq"
config_manifest_fullpath=$1

run_name=$(sed -E 's/(_premerge)?.csv//' <<< "$(basename ${config_manifest_fullpath})")

echo "#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -l software=bcbio
#PBS -l wd
#PBS -l storage=gdata/pq08

export PATH=/g/data/pq08/local/bin:\$PATH

base_dir=${base_dir}
"  > "${qsub_dir}/run_bcbio_${run_name}.qsub"

config_manifest_fullpath_final=${config_manifest_fullpath}
if [[ "${config_manifest_fullpath}" =~ .*premerge.csv ]];
then
echo "cd ${fastq_dir}" >> "${qsub_dir}/run_bcbio_${run_name}.qsub"
echo "bcbio_prepare_samples.py --out ${fastq_dir}/merged --csv ${config_manifest_fullpath}" >> "${qsub_dir}/run_bcbio_${run_name}.qsub"
config_manifest_fullpath_merged=${config_manifest_fullpath//premerge.csv/premerge-merged.csv}
config_manifest_fullpath_final=${config_manifest_fullpath_merged/_premerge-merged/}
echo "mv ${config_manifest_fullpath_merged} ${config_manifest_fullpath_final}"  >> "${qsub_dir}/run_bcbio_${run_name}.qsub"
fastq_dir="${fastq_dir}/merged"
fi

echo "
cd \"${base_dir}/analysis/bcbio\"

# Sets the system params for Gadi. Edit to include your raw samples
gadi_config_yaml=\"${base_dir}/analysis/bcbio/config/bcbio_system_normalgadi.yaml\"
# The tools you want bcbio to run. Standards are usually set for this
run_yaml=\"${base_dir}/analysis/bcbio/config/std_workflow_cancer_hg38_ensemble.yaml\"
#List of samples for BCBio to see
samples=\$(find ${fastq_dir} -iname '*.fastq.gz' | grep -f <(cut -f1 -d',' ${config_manifest_fullpath_final} | tail -n+2))

# Tell bcbio how to set up for your run. The run folder will be named after your samples csv
bcbio_nextgen.py -w template --only-metadata \${run_yaml} ${config_manifest_fullpath_final} \${samples}

# cd into the work directory that bcbio just made
cd ${base_dir}/analysis/bcbio/${run_name}/work

# Run BCBio from the work directory. BCBio will have made test_run.yaml in the config dir based
# on the name of the samples csv. First give the global config (for gadi) then give the run config
bcbio_nextgen.py \${gadi_config_yaml}  ${base_dir}/analysis/bcbio/${run_name}/config/${run_name}.yaml -n 48 -q normal -s pbspro -t ipython -r 'walltime=24:00:00;noselect;jobfs=100GB;storage=gdata/pq08' -r conmem=24 --retries 1 --timeout 900" >> "${qsub_dir}/run_bcbio_${run_name}.qsub"

