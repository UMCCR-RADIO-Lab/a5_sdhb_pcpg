#!/bin/bash
PatientIDs=()
TumourIDs=()

qsub_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/gridss/qsub"

while read -r bamfileprefix; 
do
PatientID=$(sed -E "s/(E...)-(T..)/\1/" <<< "${bamfileprefix}");
TumourID=$(sed -E "s/(E...)-(T..)/\2/" <<< "${bamfileprefix}");
PatientIDs+=("${PatientID}");
TumourIDs+=("${TumourID}");
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E*/final -iname '*-T0*-ready.bam' | xargs -I {} basename {} | sed "s/-ready.bam//")

for i in $(seq 0 1 $((${#PatientIDs[@]}-1 )));
do

PatientID="${PatientIDs[i]}"
TumourID="${TumourIDs[i]}"
NormalID="B01"
threads=24

qsub_file="${qsub_dir}/${PatientID}-${TumourID}_gridss.qsub"

echo "
#PBS -N ${PatientID}-${TumourID}_gridss
#PBS -l ncpus=${threads}
#PBS -l walltime=24:00:00
#PBS -l mem=64GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=200GB
#PBS -l wd
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au


source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/ppgl/a5/software/gridss_conda

gridss_refdir=/g/data/pq08/reference/hartwig/hmf5/GRCh38
ref_genome_fa=/g/data/pq08/local/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa

NormalBam=/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/${PatientID}/final/${PatientID}-${NormalID}/${PatientID}-${NormalID}-ready.bam
TumourBam=/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/${PatientID}/final/${PatientID}-${TumourID}/${PatientID}-${TumourID}-ready.bam

outDir=/g/data/pq08/projects/ppgl/a5/wgs/analysis/gridss/${PatientID}-${TumourID}
outVCF=\${outDir}/${PatientID}-${TumourID}.gridss.driver.vcf.gz
outVCFrepeatmasker=\${outDir}/${PatientID}-${TumourID}.gridss.repeatmasker.vcf.gz
outVCFunfiltered=\${outDir}/${PatientID}-${TumourID}.gridss.unfiltered.vcf.gz

outBAM=\${outDir}/${PatientID}-${TumourID}.assembly.bam
mkdir -p \${outDir}
gridss \
-t ${threads} \
-o \${outVCF} \
-a  \${outBAM} \
-w \${outDir} \
-r \${ref_genome_fa} \
-j \${CONDA_PREFIX}/share/gridss-2.12.2-0/gridss.jar \
-b \${gridss_refdir}/gridss/ENCFF001TDO.bed \
-c \${gridss_refdir}/gridss/gridss.properties \
\${NormalBam} \${TumourBam}

tabix \${outVCF} -p vcf

#There is an issue with repeat masker where it sometimes uses ~/ as the working directory (possibly depending on '-l wd' PBS flag). The home directory has a 10GB limit, 
# so if too many samples run this stage at one time then all will fail due to the home directory filling up. This is 
# unlikely to happen due to the naturally staggered finishing times of GRIDSS, but should be kept in mind if rerunning from this point. 

gridss_annotate_vcf_repeatmasker \
--output \${outVCFrepeatmasker} \
--jar \${CONDA_PREFIX}/share/gridss-2.12.2-0/gridss.jar \
-w /g/data/pq08/projects/ppgl/a5/wgs/analysis/gridss/${PatientID}-${TumourID} \
--rm \${CONDA_PREFIX}/bin/RepeatMasker \${outVCF}

java -Xmx8G -Dsamjdk.create_index=true -Dsamjdk.use_async_io_read_samtools=true \
-Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=true \
-Dsamjdk.buffer_size=4194304 \
-cp \${CONDA_PREFIX}/share/gridss-2.12.2-0/gridss.jar \
gridss.AnnotateInsertedSequence \
REFERENCE_SEQUENCE=\"\${gridss_refdir}/refgenomes/human_virus/human_virus.fa\" \
INPUT=\${outVCFrepeatmasker} \
OUTPUT=\${outVCFunfiltered} \
ALIGNMENT=APPEND \
WORKER_THREADS=8

in_vcf=\${outVCFunfiltered}
out_vcf=\"\${outDir}/${PatientID}-${TumourID}.gripss.somatic.vcf.gz\"
gripss -Xms24G -Xmx64G \
-tumor ${PatientID}-${TumourID} \
-reference ${PatientID}-${NormalID} \
-ref_genome \${ref_genome_fa} \
-breakend_pon \${gridss_refdir}/gridss_pon/gridss_pon_single_breakend.bed \
-breakpoint_pon \${gridss_refdir}/gridss_pon/gridss_pon_breakpoint.bedpe \
-breakpoint_hotspot \${gridss_refdir}/knowledgebases/KnownFusionPairs.bedpe \
-input_vcf \${in_vcf} \
-output_vcf \${out_vcf}

in_vcf=\"\${out_vcf}\"
out_vcf=\"\${outDir}/${PatientID}-${TumourID}.gripss.somatic.filtered.vcf.gz\"
java -Xms24G -Xmx64G -cp \$CONDA_PREFIX/share/hmftools-gripss-1.11-0/gripss.jar com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
-input_vcf \"\${in_vcf}\" \
-output_vcf \"\${out_vcf}\" 

if [ -e \"\${out_vcf}\" ]
then
rm \${outVCFrepeatmasker}
fi

" > "${qsub_file}"

done
