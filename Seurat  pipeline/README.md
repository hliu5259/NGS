# Seurat Pipeline for analysis GE matrix
Updated Nov.22 2020

## [Seurat1-1.R](https://github.com/hliu5259/singlecell/blob/master/Seurat%20%20pipeline/seurat1-1.R)
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
	<sample_id>_beforefilter.rds

### Required Argument
	-s Sample_feature_list contains sample_id, feature_min, feature_max, percent of Mito
	
## Seurat2_batch_cc.R
### Command-line:
	Rscript Seurat2_batch_cc.R -s <sample_list> -p <pattern for VAF matrix>

### Description

This script is to intergate the Seurat filtered datset to reduce batch effect, then reduce cycle effect for each sample and do cluster annotation.


### Input 
A list containing the sample_id, feature_min, feature_max, percent of Mtio as ordered
A directory containing the filtered Seurat files (_seurat_clusterd_singleR.rds). 

### sample input matrix file name
	<sample_id>_seurat_clustered_singleR.rds

### Required Argument
	-s Sample_list contains sample_id
	-p pattern for VAF matrix without sample_id (optional)

