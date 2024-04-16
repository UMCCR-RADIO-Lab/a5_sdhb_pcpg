module load bcftools 
while read -r vcf; do 

if [ ! -e ${vcf}.old ]; then
mv ${vcf} ${vcf}.old
fi

zcat ${vcf}.old | sed 's/##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Mean Mapping Quality">/##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">/' | bcftools view -O z -o ${vcf}
bcftools index -t ${vcf}
done < <(find /g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E*/final/2021* -iregex '.+/E...-.-ensemble-annotated.vcf.gz')
