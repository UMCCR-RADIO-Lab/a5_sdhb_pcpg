#!/bin/bash
basedir=/g/data/pq08/projects/ppgl/a5/wgs/analysis

patient_ids=()
tumour_ids=()

qsub_dir="${basedir}/linx/qsub"

while read -r bamfileprefix; 
do
patient_id=$(sed -E "s/(E...)-T0(.)/\1/" <<< "${bamfileprefix}");
tumour_id=$(sed -E "s/(E...)-T0(.)/\2/" <<< "${bamfileprefix}");
patient_ids+=("${patient_id}");
tumour_ids+=("${tumour_id}");
done < <(find ${basedir}/bcbio/E*/final -iname '*-T0*-ready.bam' -exec basename {} \; | sed "s/-ready.bam//")

for i in $(seq 0 1 $((${#patient_ids[@]}-1 )));
do

patient_id="${patient_ids[i]}"
tumour_id="${tumour_ids[i]}"
threads=5

echo "
#!/bin/bash
#PBS -N ${patient_id}-${tumour_id}_linx
#PBS -l ncpus=${threads}
#PBS -l walltime=00:10:00
#PBS -l mem=32GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=5GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/ppgl/a5/software/circos_conda


hartwig_ref_dir=/g/data/pq08/reference/hartwig/hmf5/GRCh38

purple_out_dir=${basedir}/purple/${patient_id}-T0${tumour_id}

linx_out_dir=${basedir}/linx/${patient_id}-T0${tumour_id}
mkdir -p \${linx_out_dir}

java -Xms8G -Xmx31G -jar /g/data/pq08/software/linx/1.16/linx.jar \
-sample \"${patient_id}-T0${tumour_id}\" \
-sv_vcf \${purple_out_dir}/${patient_id}-T0${tumour_id}.purple.sv.vcf.gz \
-purple_dir \${purple_out_dir} \
-ref_genome_version HG38 \
-output_dir \${linx_out_dir} \
-fragile_site_file \${hartwig_ref_dir}/knowledgebases/fragile_sites_hmf.csv \
-line_element_file \${hartwig_ref_dir}/knowledgebases/line_elements.csv \
-replication_origins_file \${hartwig_ref_dir}/knowledgebases/heli_rep_origins.bed \
-viral_hosts_file \${hartwig_ref_dir}/knowledgebases/viral_host_ref.csv \
-gene_transcripts_dir \${hartwig_ref_dir}/../../ensembl/ensembl_data_cache_hg38 \
-check_fusions \
-known_fusion_file \${hartwig_ref_dir}/knowledgebases/known_fusion_data.csv \
-check_drivers \
-driver_gene_panel \${hartwig_ref_dir}/knowledgebases/DriverGenePanel.tsv \
-chaining_sv_limit 0 \
-write_vis_data 2> \${linx_out_dir}/${patient_id}-T0${tumour_id}.linx.error.log 1> \${linx_out_dir}/${patient_id}-T0${tumour_id}.linx.output.log

mkdir \"\${linx_out_dir}/plot\"
mkdir \"\${linx_out_dir}/data\"

java -cp /g/data/pq08/software/linx/1.16/linx.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
-sample ${patient_id}-T0${tumour_id} \
-gene_transcripts_dir \${hartwig_ref_dir}/../../ensembl/ensembl_data_cache_hg38 \
-plot_out \"\${linx_out_dir}/plot\" \
-data_out \"\${linx_out_dir}/data\" \
-vis_file_dir \"\${linx_out_dir}\" \
-circos /g/data/pq08/projects/ppgl/a5/software/circos_conda/bin/circos \
-threads ${threads}  2> \${linx_out_dir}/${patient_id}-T0${tumour_id}.linxvis.error.log 1> \${linx_out_dir}/${patient_id}-T0${tumour_id}.linxvis.output.log
" > "${qsub_dir}/run_linx_${patient_id}-T0${tumour_id}.qsub"
done
