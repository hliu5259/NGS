# Seurat Pipeline for analysis GE matrix
Updated Nov.22 2020

## Seurat1-1.R
### Command-line:
	Rscript Seurat1-1.R -s <sample_list>

### Description

This script is to generate the Seurat original datset to visualization the feature distribution for downstream analysis. 


### Input 
A list containing the sample_id 
A directory containing the gene expression files (one per sample) 

### sample input matrix file name
	<sample_id>_wide_counts.tsv 

### Required Argument
	-s Sample list contains sample_id
	

## Seurat1-2.R
### Command-line:
	Rscript Seurat1-2.R -s <sample_feauture_list>

### Description

This script is to filter the Seurat original datset to do cluster.


### Input 
A list containing the sample_id, feature_min, feature_max, percent of Mtio
A directory containing the original Seurat files (_beforefilter.rds)

### sample input matrix file name
	<sample_id>_wide_counts.tsv 

### Required Argument
	-s Sample_feature_list contains sample_id, feature_min, feature_max, percent of Mito
