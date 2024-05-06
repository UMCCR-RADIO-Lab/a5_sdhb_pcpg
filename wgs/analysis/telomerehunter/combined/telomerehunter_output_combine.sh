basedir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter"

mkdir -p $basedir/combined

head -n1 $(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_summary.tsv' | head -n1) > $basedir/combined/telomerehunter_a5_summary.tsv
while read -r f; 
do 
tail -n+2 $f >> $basedir/combined/telomerehunter_a5_summary.tsv; 
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_summary.tsv')

head -n1 $(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_singletons.tsv' | head -n1) > $basedir/combined/telomerehunter_a5_singletons.tsv
while read -r f; 
do 
tail -n+2 $f >> $basedir/combined/telomerehunter_a5_singletons.tsv; 
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_singletons.tsv')

head -n1 $(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_TVR_top_contexts.tsv' | head -n1) > $basedir/combined/telomerehunter_a5_TVR_top_contexts.tsv
while read -r f; 
do 
tail -n+2 $f >> $basedir/combined/telomerehunter_a5_TVR_top_contexts.tsv; 
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_TVR_top_contexts.tsv')

head -n1 $(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_normalized_TVR_counts.tsv' | head -n1) | sed "s/^/PID\t/" > $basedir/combined/telomerehunter_a5_normalized_TVR_counts.tsv
while read -r f; 
do 
filename=$(basename $f); 
samplename=${filename/%_normalized_TVR_counts.tsv/}; 
tail -n+2 $f | sed "s/^/$samplename\t/" >> $basedir/combined/telomerehunter_a5_normalized_TVR_counts.tsv; 
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter -iname '*T0?_normalized_TVR_counts.tsv')
