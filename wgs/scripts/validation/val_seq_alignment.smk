# Ensure conda is sourced in bashrc if submitting with "--profile pbs"
# Load appropriate config file with --configfile

import pandas as pd


configfile: "config.yml"

sample_info = pd.read_table(config["samples_info_file"])
sample_info['key']=sample_info[['a5_id', 'gene']].agg('_'.join, axis=1)
sample_info=sample_info.set_index("key", drop=False)

rule all:
    input:
        "output/bam/a5_variant_validation.bam"
        

def make_fq_name(key, sample_info_pd, read_number):
    pd_prefix=sample_info_pd.loc[key,"prefix"]    
    fq="../fastqs/{prefix}_R{read}_001.fastq.gz".format(prefix=pd_prefix, read=read_number)
    return(fq)

rule :
    input:
        fq1=lambda wildcards: make_fq_name(key=wildcards.key, sample_info_pd=sample_info, read_number="1"),
        fq2=lambda wildcards: make_fq_name(key=wildcards.key, sample_info_pd=sample_info, read_number="2")
    output:
        bam="output/bam_intermediate/{key}.bam"
    params:
        dna_id=lambda wildcards: sample_info.loc[wildcards.key,"dna_id"],
        a5_id=lambda wildcards: sample_info.loc[wildcards.key,"a5_id"],
        gene=lambda wildcards: sample_info.loc[wildcards.key,"gene"],
        prefix=lambda wildcards: sample_info.loc[wildcards.key,"prefix"],
        ref_genome=config["ref_genome_fa"]
    resources:
        walltime="18000",
        mem="90GB",
        storage="scratch/pq08+gdata/pq08",
        jobfs="10GB"
    conda:
        "A5_variant_validation_conda_env.yaml"
    threads: 8        
    log:
        stdout="logs/{key}_alignment.stdout",
        stderr="logs/{key}_alignment.stderr"
    shell:
        "bwa mem -R \"@RG\\tID:{params.prefix}\\tLB:{params.dna_id}\\tPL:ILLUMINA\\tPU:NA\\tSM:{params.a5_id}\" "
        "-t {threads} "
        "{params.ref_genome} "
        "{input.fq1} {input.fq2} | " 
        "sambamba view -f bam -S -l0 /dev/stdin | "
        "sambamba sort -o {output.bam} /dev/stdin"
        

rule bam_merge:
    input:
       expand("output/bam_intermediate/{key}.bam", key=sample_info["key"])
    output:
        bam="output/bam/a5_variant_validation.bam"
    resources:
        walltime="86400",
        mem="32GB",
        storage="scratch/pq08+gdata/pq08",
        jobfs="10GB"
    conda:
        "A5_variant_validation_conda_env.yaml"
    log:
        stdout="logs/final_merge.stdout",
        stderr="logs/final_merge.stderr"
    threads: 24
    shell:
        "sambamba merge -t {threads} -l 9 {output.bam} {input}"


# rule xx:
#     input:
        
#     output:
#         temp()
#     params:
#         x=config["y"]
#     resources:
#         walltime="86400",
#         mem="32GB",
#         storage="scratch/pq08+gdata/pq08",
#         jobfs="10GB"
#     conda:
#         "A5_variant_validation_conda_env.yaml"
#     log:
#         stdout="logs/{prefix}_{rule}.stdout",
#         stderr="logs/{prefix}_{rule}.stderr"
#     shell:
#         ""