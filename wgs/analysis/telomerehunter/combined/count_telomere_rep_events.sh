echo "A5_ID-n_telomere_insert_events-total_svs" | tr "-" $'\t' > /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter/combined/telomere_rep_inserts_count.tsv
while read -r vcf; 
do
vcfbase=$(basename $vcf)
sample=${vcfbase/.purple.sv.vcf.gz/}
echo $sample$'\t'$(zgrep -P 'INSRMRT=.(TTAGGG|TCAGGG|TGAGGG|TTGGGG|TTCGGG|TTTGGG|ATAGGG|CATGGG|CTAGGG|GTAGGG|TAAGGG).n' $vcf | wc -l)$'\t'$(zgrep -v '^#' $vcf | wc -l);
done < <(find /g/data/pq08/projects/A5/WGS/analysis/purple/ -iname '*.purple.sv.vcf.gz') >> /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter/combined/telomere_rep_inserts_count.tsv
