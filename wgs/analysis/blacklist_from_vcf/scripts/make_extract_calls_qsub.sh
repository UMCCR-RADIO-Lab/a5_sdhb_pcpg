#!/bin/bash

position_file_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/input/positions"
mapfile position_files < <(find "${position_file_dir}" -iname '*.txt' | sort)

qsub_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/qsub"

outdir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/output/variant_observations"

for pf in "${position_files[@]}"
do

pf=$(sed -E 's/[[:space:]]+$//' <<< "${pf}")

pf_base=$(basename "${pf}")
qsub_file="${qsub_dir}"/"${pf_base/.txt/.qsub}"
out_file="${outdir}"/"${pf_base/.txt/.tsv}"

echo "
#PBS -N varcall_${pf_base}
#PBS -l ncpus=1
#PBS -l walltime=2:00:00
#PBS -l mem=8GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=50GB
#PBS -q normalbw
#PBS -m a
#PBS -M aidan.flynn@unimelb.edu.au

module load bcftools/1.12

for caller in mutect2 strelka2 vardict
do
    echo \"Processing \${caller}...\"

    mapfile vcfs < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/*/final -iname \"*-\${caller}.vcf.gz\")
    for vcf in \"\${vcfs[@]}\"
    do
        vcf=\$(sed -E 's/[[:space:]]+$//' <<< \${vcf})
        vcf_basename=\$(basename \"\${vcf}\")
        vcf_sample=\$(sed \"s/-\${caller}.vcf.gz//\" <<< \"\${vcf_basename}\")

        echo \"Sub-processing \${vcf_sample}...\"
        
        if [ \"\${caller}\" = \"strelka2\" ]; then
        bcftools view --regions-file \"${pf}\" \"\${vcf}\" | \
        bcftools query -s \"\$vcf_sample\" --format \"%CHROM\t%POS\t%ALT\t%FILTER\t[%AU:%TU:%CU:%GU:%TIR:%TAR\t%DP]\n\" | \
        awk -v caller=\${caller} -v vcf_sample=\"\${vcf_sample}\" 'BEGIN{FS=\"\t\"; OFS=\"\t\"}{print(\$1, \$2, \$3, \$4, \$5, \$6, caller, vcf_sample)}' >> ${out_file}
        else
        bcftools view --regions-file \"${pf}\" \"\${vcf}\" | \
        bcftools query -s \"\$vcf_sample\" --format \"%CHROM\t%POS\t%ALT\t%FILTER\t[%AD\t%DP]\n\" | \
        awk -v caller=\"\${caller}\" -v vcf_sample=\"\${vcf_sample}\" 'BEGIN{FS=\"\t\"; OFS=\"\t\"}{print(\$1, \$2, \$3, \$4, \$5, \$6, caller, vcf_sample)}' >> ${out_file}
        fi
    
    done
done
    " > "${qsub_file}"
done

