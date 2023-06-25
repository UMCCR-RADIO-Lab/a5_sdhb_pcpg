#!/bin/bash
basedir="/g/data/pq08/projects/A5/wts"
fastq_dir="$basedir/data/fastq"
gencode_gtf="/g/data/pq08/reference/GRCh38/gencode/gencode.v36.primary_assembly.annotation.gtf"
gencode_fa="/g/data/pq08/reference/GRCh38/gencode/GRCh38_full_analysis_set_plus_decoy_hla.fa"
threads=8

readarray -t Samples < <(tail -n+2 "$1")

for Sample in "${Samples[@]}";
do

IFS=$'\t' read -r A5_ID A5_Patient_ID prepkit SampleID LibraryID DescriptiveName FASTQ_Filename_R1 FASTQ_Filename_R2 <<< "${Sample}"

FASTQ_path_R1=$(sed "s:^:${fastq_dir}/:" <<< ${FASTQ_Filename_R1}) #append path to first fastq file
FASTQ_path_R1=$(sed "s:,:,${fastq_dir}/:g" <<< ${FASTQ_path_R1}) #append path to topup fastq files if any

FASTQ_path_R2=$(sed "s:^:${fastq_dir}/:" <<< ${FASTQ_Filename_R2}) #append path to first fastq file
FASTQ_path_R2=$(sed "s:,:,${fastq_dir}/:g" <<< ${FASTQ_path_R2}) #append path to topup fastq files if any

star_out_dir_fp="${basedir}/analysis/star/${prepkit}/${A5_ID}/firstpass"
mkdir -p "${star_out_dir_fp}"

star_fp_qsubfile="$basedir/scripts/qsub/star/${A5_ID}_${prepkit}_star_firstpass.qsub"
mkdir -p "$(dirname ${star_fp_qsubfile})"

#######################
### STAR FIRST PASS ###
#######################

echo "#PBS -N ${LibraryID}_${A5_ID}_star_firstpass
#PBS -l ncpus=${threads}
#PBS -l walltime=03:00:00
#PBS -l mem=48gb
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=150GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

module load star/2.7.9a

echo \"Starting STAR first pass...\"

#Parameters taken from https://arriba.readthedocs.io/en/latest/workflow/#demo-script and https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

STAR \
--genomeDir /g/data/pq08/reference/GRCh38/star/2.7.9 \
--readFilesCommand zcat \
--runThreadN ${threads} \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--alignSplicedMateMapLminOverLmate 0.5 \
--chimJunctionOverhangMin 10 \
--chimMainSegmentMultNmax 1 \
--chimMultimapNmax 50 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentMin 10 \
--chimSegmentReadGapMax 3 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 0 \
--limitSjdbInsertNsj 1200000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 50 \
--outFilterMultimapScoreRange 1 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMismatchNoverLmax 0.1 \
--outSAMtype None \
--outSAMmode None \
--peOverlapNbasesMin 10 \
--sjdbOverhang 150 \
--sjdbScore 2 \
--outFileNamePrefix ${star_out_dir_fp}/${A5_ID}_ \
--outTmpDir \${PBS_JOBFS}/star_tmp_${LibraryID}_firstpass \
--readFilesIn ${FASTQ_path_R1} ${FASTQ_path_R2}" > "${star_fp_qsubfile}"


# STAR \
# --runMode genomeGenerate \
# --genomeDir ${basedir}/reference/star_second_pass_genome \
# --genomeFastaFiles /g/data/pq08/reference/GRCh38/gencode/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# --sjdbOverhang 150 \
# --runThreadN ${threads} \
# --sjdbFileChrStartEnd ${basedir}/reference/star_second_pass_genome/A5_pooled_SJ.out.tab

########################
### STAR SECOND PASS ###
########################

star_out_dir_sp="${basedir}/analysis/star/${prepkit}/${A5_ID}"
star_sp_qsubfile="$basedir/scripts/qsub/star/${A5_ID}_${prepkit}_star_secondpass.qsub"
mkdir -p "$(dirname ${star_sp_qsubfile})"

echo "
#PBS -N ${LibraryID}_${A5_ID}_star_secondpass
#PBS -l ncpus=${threads}
#PBS -l walltime=3:00:00
#PBS -l mem=48gb
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=150GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

module load star/2.7.9a
echo \"Starting STAR second pass...\"

#Parameters taken from https://arriba.readthedocs.io/en/latest/workflow/#demo-script and https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

STAR \
--genomeDir ${basedir}/reference/star_second_pass_genome_${prepkit} \
--readFilesCommand zcat \
--runThreadN ${threads} \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--alignSoftClipAtReferenceEnds Yes \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--alignSplicedMateMapLminOverLmate 0.5 \
--chimJunctionOverhangMin 10 \
--chimMainSegmentMultNmax 1 \
--chimMultimapNmax 50 \
--chimOutType WithinBAM HardClip \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentMin 10 \
--chimSegmentReadGapMax 3 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 0 \
--limitSjdbInsertNsj 1200000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 50 \
--outFilterMultimapScoreRange 1 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMismatchNoverLmax 0.1 \
--outSAMattributes NH HI NM MD AS XS \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--peOverlapNbasesMin 10 \
--sjdbOverhang 150 \
--sjdbScore 2 \
--outFilterIntronMotifs None \
--outSAMtype BAM SortedByCoordinate \
--outSAMheaderHD @HD VN:1.4 \
--outSAMattrRGline ID:${A5_Patient_ID}_${LibraryID} SM:${A5_ID} \
--outFileNamePrefix ${star_out_dir_sp}/${A5_ID}_ \
--outTmpDir \${PBS_JOBFS}/star_tmp_${LibraryID}_secondpass \
--readFilesIn ${FASTQ_path_R1} ${FASTQ_path_R2}


echo \"Starting BAM indexing ...\"

module load samtools
samtools index ${star_out_dir_sp}/${A5_ID}_Aligned.sortedByCoord.out.bam" > "${star_sp_qsubfile}"

#############
### HTSEQ ###
#############

htseq_qsubfile="$basedir/scripts/qsub/htseq/${A5_ID}_${prepkit}_htseq.qsub"
mkdir -p "$(dirname ${htseq_qsubfile})"

htseqoutdir="${basedir}/analysis/htseq/${prepkit}"

mkdir -p "${htseqoutdir}"
mkdir -p "${htseqoutdir}/gene"
mkdir -p "${htseqoutdir}/exon"

echo "#PBS -N ${LibraryID}_${A5_ID}_htseq
#PBS -l ncpus=2
#PBS -l walltime=10:00:00
#PBS -l mem=24gb
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=10GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh

conda activate /g/data/pq08/software/htseq_conda

echo \"Starting Gene counts ...\"

printf \"gene_id\\tgene_symbol\\tcount\n\" > ${htseqoutdir}/gene/${A5_ID}.gene.counts 

htseq-count \
-f bam \
-r pos \
-s no \
-a 10 \
-t exon \
-i gene_id \
-m intersection-nonempty \
--additional-attr=gene_name \
${star_out_dir_sp}/${A5_ID}_Aligned.sortedByCoord.out.bam \
${gencode_gtf} >> ${htseqoutdir}/gene/${A5_ID}.gene.counts &

echo \"Starting Exon counts ...\"

printf \"exon_id\\tgene_symbol\\tgene_id\\ttranscript_id\\texon_number\\tcount\n\" > ${htseqoutdir}/exon/${A5_ID}.exon.counts

htseq-count \
-f bam \
-r pos \
-s no \
-a 10 \
--additional-attr=gene_id \
--additional-attr=transcript_id \
--additional-attr=gene_name \
--additional-attr=exon_number \
-t exon \
-i exon_id \
-m intersection-nonempty \
${star_out_dir_sp}/${A5_ID}_Aligned.sortedByCoord.out.bam \
${gencode_gtf} >> ${htseqoutdir}/exon/${A5_ID}.exon.counts &

wait" > "${htseq_qsubfile}"

##############
### ARRIBA ###
##############

arriba_qsubfile="$basedir/scripts/qsub/arriba/${A5_ID}_${prepkit}_arriba.qsub"
mkdir -p "$(dirname ${arriba_qsubfile})"

arribaoutdir="$basedir/analysis/arriba/${prepkit}/${A5_ID}"
mkdir -p "${arribaoutdir}"

echo "#PBS -N ${LibraryID}_${A5_ID}_arriba
#PBS -l ncpus=1
#PBS -l walltime=12:00:00
#PBS -l mem=24gb
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=150GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh

conda activate /g/data/pq08/software/arriba_conda

arriba \
-x ${star_out_dir_sp}/${A5_ID}_Aligned.sortedByCoord.out.bam \
-o ${arribaoutdir}/${A5_ID}_fusions.tsv \
-O ${arribaoutdir}/${A5_ID}_fusions.discarded.tsv \
-a ${gencode_fa} \
-g ${gencode_gtf} \
-b \$CONDA_PREFIX/var/lib/arriba/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \
-k \$CONDA_PREFIX/var/lib/arriba/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz \
-t \$CONDA_PREFIX/var/lib/arriba/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz \
-p \$CONDA_PREFIX/var/lib/arriba/protein_domains_hg38_GRCh38_v2.1.0.gff3" > /dev/null #"$arriba_qsubfile"

done

