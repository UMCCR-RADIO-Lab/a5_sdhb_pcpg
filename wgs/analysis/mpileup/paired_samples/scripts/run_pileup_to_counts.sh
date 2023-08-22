#pileup_to_counts script requires python>3.8
module load python3/3.11.0 parallel

base_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup"
pileup2counts_script="${base_dir}/blacklist/scripts/pileup_to_counts.py"

patients=$(find "${base_dir}/paired_samples/positions" -iname 'E*.txt' | sort | xargs -I {} basename {} | sed -E 's/.txt//' | xargs)

for patient in ${patients[@]};
do

samples=$(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/${patient}/final/ -iname '*-ready.*am' | sort | xargs -I {} basename {} | sed -E 's/-ready.(cr|b)am//' | xargs)
chr_pos_file="${base_dir}/paired_samples/positions/${patient}.txt"
pileup_file="${base_dir}/paired_samples/pileups/${patient}.pileup"
count_file="${base_dir}/paired_samples/count_summaries/${patient}.txt"

echo "Processing ${patient}..."
python3 ${pileup2counts_script} --pileup ${pileup_file} --outputValue AltRefTotal --variantCalls <(cut ${chr_pos_file} -f1,2,4,5) --pileupSamples ${samples} > ${count_file} 

done
