## Purpose
This repository contains the code for all analysis presented in the publication **Multi-omic analysis of SDHB-deficient pheochromocytomas and paragangliomas identifies metastasis and treatment-related molecular profiles**.

## Structure
The main folders are named for a data-type including: 
- Whole Genome Sequencing (wgs) 
- Whole Transcriptome Sequencing (wts)
- Illumina EPIC methylation profiling (methylation)
- Small-RNA sequencing (small_rna)
- Single-nuclei ATAC sequencing (sn_atac_seq)
- Single-nuclei RNA sequencing (sn_rna_seq)
- Sample annotation (sample_annotation)
- Tertiary analysis involving integration of multiple data types (tertiary)

## Sub-Structure
Within each data-type folder there can be the following subfolders:
- analysis: Contains code for primary analysis of data, such as sequence alignment or feature counting
- scripts/
  - data_loaders: Dataloader scripts encapsulate data loading and preprocessing steps, such as reading in files, normalisation, and filtering
  - data_mergers: Data merger scripts combine multiple datasets of the same type from different sources, e.g publicly available datasets 
  - figures: these scripts contain the code to reproduce the figures as presented in the manuscript
 
#Additional folders
- external_dataloaders: Contains data loaders to load and preprocess publicly available datasets
- utility_scripts: Helper scripts used by one or more other scripts 
