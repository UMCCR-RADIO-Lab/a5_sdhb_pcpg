  #!/usr/bin/bash


base_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/analysis/cellranger_hg38/"
qsub_dir=${base_dir}/qsub
out_dir=${base_dir}/intronic_counted

while read -r sample;
do 
    qsub_file="${qsub_dir}/${sample}_cranger_count_introns.qsub"

    echo "#!/bin/bash
    #PBS -N cellranger_${sample}
    #PBS -l ncpus=16
    #PBS -l walltime=12:00:00
    #PBS -l mem=128GB
    #PBS -l storage=scratch/pq08+gdata/pq08
    #PBS -l jobfs=50GB
    #PBS -q normal
    #PBS -m a
    #PBS -M aidan.flynn@unimelb.edu.au

    export PATH=\"${PATH}:/g/data/pq08/projects/ppgl/a5/software/cellranger/cellranger-7.2.0/\"
    
	cd ${out_dir}
    
	 " > "${qsub_file}"
	 
	
	if [ "${sample}" = "E140-1" ] || [ "${sample}" = "E143-1" ]  || [ "${sample}" = "E171-1" ] || [ "${sample}" = "E225-1" ];
	then	
	echo "cellranger count --id=${sample} --r1-length 26 --fastqs=/g/data/pq08/projects/ppgl/a5/sn_rna_seq/raw_data/fastq/${sample} --include-introns true --sample=${sample} --transcriptome=/g/data/pq08/reference/cellranger/refdata-gex-GRCh38-2020-A/" >> "${qsub_file}"
	else
	echo "cellranger count --id=${sample} --fastqs=/g/data/pq08/projects/ppgl/a5/sn_rna_seq/raw_data/fastq/${sample} --include-introns true --sample=${sample} --transcriptome=/g/data/pq08/reference/cellranger/refdata-gex-GRCh38-2020-A/" >> "${qsub_file}";
	fi;
	
 done < <(ls /g/data/pq08/projects/ppgl/a5/sn_rna_seq/raw_data/fastq)
