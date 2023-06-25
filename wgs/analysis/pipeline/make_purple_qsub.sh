#!/bin/bash
basedir=/g/data/pq08/projects/A5/WGS/analysis

patient_ids=()
tumour_ids=()

qsub_dir=${basedir}/purple/qsub

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
#PBS -N ${patient_id}-${tumour_id}_purple
#PBS -l ncpus=${threads}
#PBS -l walltime=00:30:00
#PBS -l mem=32GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=5GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/A5/software/circos_conda


hartwig_ref_dir=/g/data/pq08/reference/hartwig/hmf5/GRCh38

purple_out_dir=${basedir}/purple/${patient_id}-T0${tumour_id}
mkdir -p \$purple_out_dir

cobalt_result_dir=${basedir}/bcbio/${patient_id}/umccrised/work/${patient_id}-${tumour_id}__${patient_id}-T0${tumour_id}/purple/cobalt/
amber_result_dir=${basedir}/bcbio/${patient_id}/umccrised/work/${patient_id}-${tumour_id}__${patient_id}-T0${tumour_id}/purple/amber/

mkdir -p \$purple_out_dir/amber
mkdir -p \$purple_out_dir/cobalt

##Symlink UMCCRise output to remove \"${patient_id}-${tumour_id}__\" prefix which messes up purple's file dectection
while read -r f; 
do
linkname=\$(basename \${f} | sed \"s/${patient_id}-${tumour_id}__//\")
ln -s \"\${f}\" \"\${purple_out_dir}/cobalt/\${linkname}\"
done < <(find \"\${cobalt_result_dir}\")

while read -r f; 
do
linkname=\$(basename \${f} | sed \"s/${patient_id}-${tumour_id}__//\")
ln -s \"\${f}\" \"\${purple_out_dir}/amber/\${linkname}\"
done < <(find \"\${amber_result_dir}\")

strelka_vcf=${basedir}/bcbio/${patient_id}/final/${patient_id}-T0${tumour_id}/${patient_id}-T0${tumour_id}-strelka2.vcf.gz
strelka_vcf_ad=\"\${purple_out_dir}/${patient_id}-T0${tumour_id}-strelka2_ad.vcf.gz\"

java -Xmx4G -cp /g/data/pq08/software/purple/3.1/purple.jar com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth \
-in \${strelka_vcf} \
-out \${strelka_vcf_ad} 2> \"\${purple_out_dir}/${patient_id}-T0${tumour_id}.purple.error.log\" 1> \"\${purple_out_dir}/${patient_id}-T0${tumour_id}.purple.output.log\"

gripss_vcf=\"/g/data/pq08/projects/A5/WGS/analysis/gridss/${patient_id}-T0${tumour_id}/${patient_id}-T0${tumour_id}.gripss.somatic.vcf.gz\"
gripss_filtered_vcf=\"/g/data/pq08/projects/A5/WGS/analysis/gridss/${patient_id}-T0${tumour_id}/${patient_id}-T0${tumour_id}.gripss.somatic.filtered.vcf.gz\"

java -Xms8G -Xmx31G  -jar /g/data/pq08/software/purple/3.1/purple.jar \
-reference ${patient_id}-B01 \
-tumor ${patient_id}-T0${tumour_id} \
-output_dir \"\${purple_out_dir}\" \
-amber \"\${purple_out_dir}/amber\" \
-cobalt \"\${purple_out_dir}/cobalt\" \
-somatic_vcf \"\${strelka_vcf_ad}\" \
-structural_vcf \"\${gripss_filtered_vcf}\" \
-sv_recovery_vcf \"\${gripss_vcf}\" \
-circos /g/data/pq08/projects/A5/software/circos_conda/bin/circos \
-ref_genome /g/data/pq08/projects/flynna/acquired_resist/reference/hg38.fa \
-somatic_hotspots \"\${hartwig_ref_dir}/knowledgebases/KnownHotspots.vcf.gz\" \
-gc_profile \"\${hartwig_ref_dir}/gc/GC_profile.1000bp.cnp\" \
-driver_gene_panel \"\${hartwig_ref_dir}/knowledgebases/DriverGenePanel.tsv\" \
-driver_catalog \
-threads ${threads} 2> \"\${purple_out_dir}/${patient_id}-T0${tumour_id}.purple.error.log\" 1> \"\${purple_out_dir}/${patient_id}-T0${tumour_id}.purple.output.log\"
" > "${qsub_dir}/run_purple_${patient_id}-T0${tumour_id}.qsub"

done