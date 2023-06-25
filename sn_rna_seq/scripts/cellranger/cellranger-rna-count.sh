#!/bin/bash
#PBS -P pq08
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/pq08
#PBS -l mem=160GB
#PBS -l ncpus=48
#PBS -l wd

CELLR=/g/data/pq08/software/cellranger/6.1.1/cellranger-6.1.1/cellranger
REFERENCE=/g/data/pq08/reference/cellranger/refdata-gex-GRCh38-2020-A
FASTQS=/g/data/pq08/projects/bowenb/a5_project/snRNA-seq/rna_fqs
MEM=160
CORES=48

# samples: E140-1=LPRJ210351, E143-1=LPRJ210349, E156-1=PRJ190411_E156-1, E166-1=PRJ190412_E166-1, E171-1=LPRJ210350, E225-1=LPRJ210348
# NAM021=PRJ180574_NAM021, NAM025=PRJ180573_NAM025, NPG103=PRJ190607_NPG103, E146-1=S126=PRJ190608_S126, E197-1=TMC41=PRJ190609_TMC41
# E123-1=VPH58=PRJ190610_VPH58, 
# E019=VPH23T=PRJ180542_VPH23T, P018-PGL1=PGL1=10x_nc_5prime_Lbry_PGL1, P018-PGL3=PGL3=PRJ180544_PGL3

$CELLR count --id=P018-PGL3 \
             --transcriptome=$REFERENCE \
             --fastqs=$FASTQS/PGL3 \
             --sample=PRJ180544_PGL3 \
             --localcores=$CORES \
             --localmem=$MEM \
	         --r1-length=26 \
             --r2-length=91

# Notes: when processing the 3' snRNA-seq chemistry samples (this was E140-1, E143-1, E171-1, E225-1) set --r1-length=26 to trim the read length to 26bp
# without this, got ERROR: We detected a mixture of different R1 lengths ([26-151])
# The libraries sequenced at UMCCR  were sequenced with 150bp paired-end sequencing, therefore I trim r2-length to 91bp as recommended by 10x genomics at https://kb.10xgenomics.com/hc/en-us/articles/360016221712-Can-I-sequence-longer-than-91-98bp-on-Read-2-for-Single-Cell-3-and-5-Gene-Expression-libraries-
# EDIT: doing all samples (5 and 3') with  --r1-length=26 --r2-length=91
