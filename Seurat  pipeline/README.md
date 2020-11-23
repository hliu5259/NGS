#Seurat Pipeline for analysis GE matrix
Updated Nov.22 2020


### Command-line:
	Rscript seurat1-1.R -s <sample_list>

## Description

This script is to generate the Seurat original datset to visualization the feature distribution for downstream analysis. 

## Required Argument
	-s Sample list contains sample_id

## Input 
A list containing the sample_id 
A directory containing the gene expression files (one per sample) 

### sample input matrix file name
	<sample_id>_wide_counts.tsv 

	

