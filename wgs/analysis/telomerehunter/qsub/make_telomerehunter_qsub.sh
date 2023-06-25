for s in E019-T01 E120-T01 E121-T01 E122-T02 E122-T01 E123-T01 E124-T01 E125-T01 E126-T01 E127-T01 E128-T01 E128-T02 E129-T01 E130-T01 E131-T01 E132-T02 E132-T01 E133-T01 E134-T01 E135-T01 E136-T01 E138-T01 E140-T01 E141-T01 E142-T01 E143-T01 E143-T02 E143-T03 E144-T01 E145-T01 E146-T01 E146-T02 E147-T01 E148-T01 E149-T01 E150-T01 E151-T01 E152-T01 E153-T01 E154-T01 E155-T01 E156-T01 E157-T01 E158-T02 E158-T01 E159-T03 E159-T04 E159-T01 E159-T02 E160-T01 E161-T01 E162-T01 E163-T01 E164-T01 E165-T01 E166-T02 E166-T01 E167-T02 E167-T01 E168-T01 E169-T02 E169-T01 E170-T01 E171-T01 E179-T01 E180-T01 E182-T01 E183-T01 E184-T01 E185-T01 E186-T01 E188-T01 E189-T01 E190-T01 E192-T01 E193-T01 E194-T01 E195-T01 E196-T01 E197-T01 E198-T01 E199-T01 E200-T01 E201-T01 E202-T01 E203-T01 E204-T01 E223-T01 E224-T01 E225-T02 E225-T01 E226-T01 E227-T01 E228-T01 E229-T01 E229-T02 E230-T01 E231-T01 E233-T01;
do

basedir="/g/data/pq08/projects/A5/WGS/analysis"
patient=$(sed 's/-T..//' <<< $s)
tumour=$(sed 's/E...-//' <<< $s)
tumour_bam="${basedir}/bcbio/${patient}/final/${patient}-${tumour}/${patient}-${tumour}-ready.bam"
control_bam="${basedir}/bcbio/${patient}/final/${patient}-B01/${patient}-B01-ready.bam"
outdir="${basedir}/telomerehunter/"
cytobands="/g/data/pq08/projects/A5/software/telomerehunter_conda/reference/GRCh38.cytoBand_primary.txt"
echo "
#PBS -N ${patient}-${tumour}_telomerehunter
#PBS -l ncpus=2
#PBS -l walltime=8:00:00
#PBS -l mem=32GB
#PBS -l storage=scratch/pq08+gdata/pq08
#PBS -l jobfs=10GB
#PBS -q normal
#PBS -m ae
#PBS -M aidan.flynn@unimelb.edu.au

source /g/data/pq08/projects/flynna/software/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/pq08/projects/A5/software/telomerehunter_conda/

telomerehunter -ibt ${tumour_bam} -ibc ${control_bam} -o "${outdir}" \
-p ${patient}-${tumour} -b $cytobands -pl \
-pff png -p1 -p2 -p3 -p4 -p5 -p6 -p7 -prc
" > ${patient}-${tumour}_telomerehunter.qsub
done
