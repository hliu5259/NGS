# Seurat_batch_cc.R
# LAST UPDATED ON DEC.22 2020
# interget samples together to reduce batch effect and cycle cell effect; split the GE, PC, VAF matirx by cell type and filtered

# set for mgpc environment
options(bitmapType="cairo")

# load package
print('loading required package data.table, tidyverse, Seurat, SingleR, factoextra')
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
suppressMessages(library(SingleR, quietly = TRUE))
suppressMessages(library(factoextra, quietly = TRUE))


# manages command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %%2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>% mutate_all(as.character)
    
  # sample list
  sample <<- arg_df$value[arg_df$flag == "-s"]
  sample_list <<- fread(sample_id, header = F)
  sample_id <<- sample_list$V1
  # check the argument
  if (nrow(sample_list) ==0) stop('sample id supplied incorrectly')
  # check the import dataset
  for (i in 1:nrow(sample_list)){
    if (!file.exists(paste0(sample_id[i], '_seurat_clustered_singleR.rds'))
            print(paste0(sample_id[i], '_seurat_clustered_singleR.rds NOT exist'))
    else print(paste0('gonna process ', sample_id[i]))
    }
  # vaf matrix pattern
  vaf_pattern <<- arg_df$value[arg_df$flag == "-p"]
  vaf_matrix_load <- paste0(sample_list$V1, vaf_pattern)
  if ((arg_df$value[arg_df$flag == "-p"]) !=0) {print ('with provided VAF matirx')}
  else print('without provide VAF matrix')
  
  for (i in 1:length(vaf_matrix_load)){
    if (!file.exists(vaf_matrix_load[i])) print(paste0(vaf_matrix_load[i], ' NOT exist'))
    else print(paste0('gonna process ', vaf_matrix_load[i]))
  }
}

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# loading _seurat_clustered_singleR.rds file
print('import the radom _seurat_clustered_singleR.rds file to reduce batch effect')

sample_list_name <- paste0(sample_list$V1, '_seurat_clustered_singleR.rds')
list_data <- list()
for (i in (1:length(sample_list_name))){
  list_data[[i]] <- readRDS(sample_list_name[[i]])
  print(paste0('loading ...\n', sample_list_name[[i]]))
}

# set up intergate feature across all samples as 2000
features <- SelectIntegrationFeatures(object.list = list_data, nfeatures = 2000)
options(future.globals.maxSize = 9768*1024^2)
list_data <- PrepSCTIntegration(object.list = list_data, anchor.features = features,
                           verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list_data, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
sample_comb <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

# analysis the intergate dataset
DefaultAssay(sample_comb) <- "integrated"
sample_comb <- RunPCA(sample_comb, verbose = FALSE)
sample_comb <- RunUMAP(sample_comb, verbose = FALSE, dim = 1:30)
sample_comb<- FindNeighbors(sample_comb, dims = 1:30) #,  k.param = 10
sample_comb<- FindClusters(sample_comb, resolution = 0.2)

# annotation cell type by using SingleR
rna_re <- BlueprintEncodeData()
cluster <- sample_comb@active.ident
b <- GetAssayData(sample_comb)
result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
sample_comb[["SingleR.cluster.labels"]] <-
  result_cluster$labels[match(sample_comb[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]

p = DimPlot(sample_comb, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
png('sample_comb_singleR.png',width = 450, height = 400)
print(p)
dev.off()

p = DimPlot(sample_comb, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
png('sample_comb_cluster_seurat_umap.png',width = 450, height = 400)
print(p)
dev.off()

# plot heatmap
png('sample_comb_filtered_heatmap.png'), width = 450, height = 300)
print(plotScoreHeatmap(result_cluster))
dev.off()
# save combine sample
saveRDS(sample_comb, 'sample_comb_afterPCA.rds')

# split the intergated dataset
print('split the intergated samples')
list_new <- SplitObject(sample_comb, split.by = "orig.ident")
print(paste0('split lenght ',length(list_new)))
# loading cell cycle marker gene
if (file.exists('stage_s.txt')){
    print('loading S stage marker genes')
    s.genes <- fread('stage_s.txt')
    s.genes <- s.genes$Gene
} else { print('No S stage marker genes attached')}

if (file.exists('stage_G2M.txt')){
    print('loading G2M stage marker genes')
    g2m.genes <- fread('stage_G2M.txt')
    g2m.genes <- g2m.genes$Gene
} else {print('No G2M stage marker genes attached')}

# function
lentop20 <- function(x){
  for (i in (2:(length(x)-1))) {
    len <- length(which(x[i] < x[-1]))/(length(x)-2)
    return(len)
  }}

lenbotton20 <- function(x){
  for (i in (2:(length(x)-1))) {
    len <- length(which(x[i] > x[-1]))/(length(x)-2)
    return(len)
  }}

# database choose for SingleR
print('loading dataset for SingleR')
rna_re <- BlueprintEncodeData()

# reduce batch and cc for all three samples
print('processing three samples')
for (i in (1:length(list_new))){
  sample_choose <- unique(list_new[[i]]@meta.data$orig.ident)
  # write outputs
  if (!dir.exists(sample_choose)) {
    print('Creating output directory...\n')
    dir.create(sample_choose, showWarnings = FALSE)
  # set up working direction
  setwd(sample_choose)
  print(paste0('plotting batch effect for ', sample_choose))
  sample_1 <- list_new[[i]]
  DefaultAssay(sample_1) <- "integrated"
  
  b <- GetAssayData(sample_1)
  cluster <- sample_1@active.ident
  result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
  sample_1[["SingleR.cluster.labels"]] <-
    result_cluster$labels[match(sample_1[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]
  
  p = DimPlot(sample_1, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE)
  png(paste0(sample_choose,'_afterbatch_singleR.png'),width = 450, height = 400)
  print(p)
  dev.off()
  
  # plot heatmap
  png(paste0(sample_choose,"_afterbatch_heatmap.png"), width = 450, height = 300)
  print(plotScoreHeatmap(result_cluster))
  dev.off()
  
  print(paste0('start processing cell cycle effect for ', sample_choose))
  # cell cycle calculate
  
  sample_1 <- CellCycleScoring(sample_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  sample_1 <- RunPCA(sample_1)
  
  print('plotting cell cycle effect')
  p= ElbowPlot(object = sample_1)
  png(paste0(sample_choose,'_pca.png'),width = 450, height = 400)
  print(p)
  dev.off()
  
  sample_1 <- RunUMAP(sample_1, dims = 1:20)
  
  p = DimPlot(sample_1, group.by = "Phase", reduction = "umap",label = TRUE)
  png(paste0(sample_choose, '_afterbatchwithcc.png'),width = 450, height = 400)
  print(p)
  dev.off()
  
  sample_1 <- ScaleData(sample_1, vars.to.regress = c("S.Score", "G2M.Score"),
                        use.umi = FALSE, do.scale = FALSE, do.center =  FALSE )
  sample_1 <- RunPCA(sample_1)
  sample_1 <- RunUMAP(sample_1,dims = 1:20)
  
  p = DimPlot(sample_1, group.by = "Phase", reduction = "umap",label = TRUE)
  png(paste0(sample_choose,'_afterbatchcc.png'),width = 450, height = 400)
  print(p)
  dev.off()
  
  print('cell type annotation after batch and cc effect')
  sample_1<- FindNeighbors(sample_1, dims = 1:30)#, k.param = 10)
  sample_1<- FindClusters(sample_1, resolution = 0.2)
  cluster <- sample_1@active.ident
  #annotation cell type
  b <- GetAssayData(sample_1)
  
  result_cluster <- SingleR(test = b, ref = rna_re, labels = rna_re$label.fine, method="cluster", clusters = cluster)
  sample_1[["SingleR.cluster.labels"]] <-
    result_cluster$labels[match(sample_1[[]]["seurat_clusters"]$seurat_clusters, rownames(result_cluster))]
  
  p = DimPlot(sample_1, group.by =  "SingleR.cluster.labels", reduction = "umap", label = TRUE) + labs(title = sample_1)
  png(paste0(sample_choose,'_afterbatchcc_singleR.png'),width = 450, height = 400)
  print(p)
  dev.off()
  
  png(paste0(sample_choose,"_afterbatchcc_heatmap.png"), width = 450, height = 300)
  print(plotScoreHeatmap(result_cluster))
  dev.off()
  
  ide <- data.frame(sample_1@active.ident)
  row.names(ide) <- sapply(row.names(ide), function(x) strsplit(x, '_')[[1]][1])
  fwrite(ide, paste0(sample_choose,'_cluster_Seurat.txt'),sep = "\t",quote = F, row.names = T, col.names = T)
  saveRDS(sample_1, file= paste0(sample_choose,'_done.rds'))
  
  cell_list <- sample_1[["SingleR.cluster.labels"]]
  row.names(cell_list) <- sapply(row.names(cell_list), function(x) strsplit(x, '_')[[1]][1])
  cell_list$SingleR.cluster.labels <- gsub(' ', '', cell_list$SingleR.cluster.labels)
  cell_list$SingleR.cluster.labels <- gsub('\\+', '_', cell_list$SingleR.cluster.labels)
  fwrite(cell_list, paste0(sample_choose, '_cell_list.txt'),sep = "\t",quote = F, row.names = T, col.names = T)
  
  GE_matrix <- data.frame(sample_1@assays$integrated@data)
  colnames(GE_matrix) <- sapply(colnames(GE_matrix), function(x) strsplit(x, '_')[[1]][1])
  GE_matrix_output <- GE_matrix %>% rownames_to_column()
  names(GE_matrix_output)[1] <- 'gene_id'
  fwrite(GE_matrix_output, paste0(sample_choose,'_int_GE_matrix.txt'),sep = "\t",quote = F, row.names = F, col.names = T)
  
  print('split GEbycelltype')
  # write output for GE by cell type
  if (!dir.exists(paste0(sample_choose,'_cell_type_GE'))) {
    print('Creating output directory...\n')
    dir.create(paste0(sample_choose,'_cell_type_GE'), showWarnings = FALSE)

  cell_list <- cell_list %>% rownames_to_column() 
  # load cell type list
  cell_type <- as.list(unique(cell_list$SingleR.cluster.labels))
  
  for (m in (1:length(cell_type))){
    cell_bycluster <- cell_list[cell_list$SingleR.cluster.labels==cell_type[[m]],]
    cell_matrix <- GE_matrix[,names(GE_matrix) %in% cell_bycluster$rowname]
    #colnames(cell_matrix) <- sapply(colnames(cell_matrix), function(x) strsplit(x, '_')[[1]][1])
    cell_matrix <- cell_matrix %>% rownames_to_column()
    names(cell_matrix)[1] <- 'gene_id'
    fwrite(cell_matrix, paste0(sample_choose,'_cell_type_GE/', sample_choose,'_', cell_type[m],'_GE_matrix_celltype.txt'), quote = F, row.names = F, sep = '\t')
    
    print(paste0('filter GE matrix for ', sample_choose))
    ##filter the 80% of gene across cells in the top20 percent of the range
    cell_matrix$top20 <- apply(cell_matrix[,-1], 1, function(x) max(x) - 0.2*(max(x) - min(x)))
    cell_matrix$percent20 <- apply(cell_matrix[,-1], 1, function(x) lentop20(x))
    
    cell_matrix <- cell_matrix %>% filter(percent20 < 0.8)
    cell_matrix <- cell_matrix[ ,!(names(cell_matrix) %in% c('top20','percent20'))]
    
    ##filter the 80% of gene across cells in the botton20 percent of the range
    cell_matrix$botton20 <- apply(cell_matrix[,-1], 1, function(x) min(x) + 0.2*(max(x) - min(x)))
    cell_matrix$percent20 <- apply(cell_matrix[,-1], 1, function(x) lenbotton20(x))
    
    cell_matrix <- cell_matrix %>% filter(percent20 < 0.8)
    cell_matrix <- cell_matrix[ ,!(names(cell_matrix) %in% c('botton20','percent20'))]
    fwrite(cell_matrix, paste0(sample_choose, '_cell_type_GE/',sample_choose,'_', cell_type[m],'_GE_matrix_filtered20.txt'), row.names = F, quote = F, sep = '\t')
    
    print(paste0('top 15pc from GE Matrixs by cell type for ', sample_choose))
    # write output for PC
    if (!dir.exists(paste0(sample_choose,'_PC'))) {
      print('Creating output directory...\n')
      dir.create(paste0(sample_choose,'_PC'), showWarnings = FALSE)
    # build up PC
    pca <- prcomp(na.omit(cell_matrix[,-1]))
    pc <- pca$rotation
    pc <- t(pc) %>% as.data.frame()
    pc <- pc %>% rownames_to_column()
    names(pc)[1] <- 'id'
    fwrite(pc[1:15,], paste0(sample_choose,'_PC/', sample_choose, '_', cell_type[m],'_top15_pcs.txt'), row.names = F, sep = '\t')
    #pc <- pca$rotation
    p=fviz_eig(pca, addlabels = TRUE) + labs(title=cell_type[m])
    pdf(paste0(sample_choose,'_PC/', sample_choose, '_' , cell_type[i], '_top15pcs.pdf'))
    print(p)
    dev.off()
    print(paste0('GE processing done for ', sample_choose))
  }
  
  ## split VAF matrix
  vaf_matrix_load <- paste0('../', sample_choose, vaf_pattern)
  if (file.exists(vaf_matrix_load)) {
      print(paste0('loading VAF matrix for ', sample_choose))
      vaf_matrix <- data.frame(fread(vaf_matrix_load), row.names = 1)
      # build up VAF_loc
      print(paste0('build up VAF loc matrix for ', sample_choose))
      df_loc <- str_split_fixed(row.names(vaf_matrix), ':', 2) %>% as.data.frame()
      tmp <- str_split_fixed(df_loc$V2, '_',2) %>% as.data.frame()
      df_loc <- cbind(row.names(vaf_matrix), df_loc$V1)
      df_loc <- cbind(df_loc, tmp$V1) %>% as.data.frame()
      df_loc$V2 <- paste0('chr', df_loc$V2)
      names(df_loc) <- c('SNV','CHROM','POS')
      fwrite(df_loc, paste0(sample_choose, "_VAF_loc.txt"), quote = F, row.names = F, sep = '\t', na ='NA')
      for (m in (1:length(cell_type))){
        print(paste0('split VAF matrix by cell type for ', sample_choose))
        # write output for VAF by cell type
        if (!dir.exists(paste0(sample_choose,'_cell_type_VAF'))) {
          print('Creating output directory...\n')
          dir.create(paste0(sample_choose,'_cell_type_VAF'), showWarnings = FALSE)

        # split vaf matrix by cell type
        cell_bycluster <- cell_list[cell_list$SingleR.cluster.labels==cell_type[[m]],]
        vaf_cell_matrix <- vaf_matrix[,names(vaf_matrix) %in% cell_bycluster$rowname]
        vaf_cell_matrix <- vaf_cell_matrix %>% rownames_to_column()
        names(vaf_cell_matrix)[1] <- 'SNV'
        fwrite(vaf_cell_matrix, paste0(sample_choose, '_cell_type_VAF/', sample_choose,'_', cell_type[m],'_VAF_matrix_barcode_celltype.txt'), quote = F, row.names = F, sep = '\t', na ='NA')
        # filter VAF matrix
        vaf_cell_matrix <- vaf_cell_matrix %>% column_to_rownames('SNV')
        # add columns with
        ##1.1) number of nonNA cells per SNV
        vaf_cell_matrix$nonNA <- rowSums(!is.na(vaf_cell_matrix))
        
        ##1.2) mean VAF
        vaf_cell_matrix$mean <- rowMeans(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA"))], na.rm = T)
        
        ##1.3) median VAF
        vaf_cell_matrix$median <- rowMedians(as.matrix(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean"))]), na.rm = T)
        
        ##1.4) percentage of nonNA cells w VAF<0.25
        vaf_cell_matrix$nonNA25 <- rowSums(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median"))]<0.25, na.rm = T)/nrow(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median"))])
        
        ##1.5) percentage of nonNA cells w VAF>0.75
        vaf_cell_matrix$nonNA75 <- rowSums(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median","nonNA25"))]>0.75, na.rm = T)/nrow(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median","nonNA25"))])
        
        ##1.6) percentage of nonNA cells w VAF between 0.4 and 0.6
        vaf_cell_matrix$nonNA46 <- rowSums(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median","nonNA","nonNA75"))] >0.4 &
                                        vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median","nonNA","nonNA75"))] <0.6, na.rm = T)/nrow(vaf_cell_matrix[,-which(colnames(vaf_cell_matrix) %in% c("nonNA","mean","median","nonNA","nonNA75"))])
        
        
        #2) Filter the matrix:
        ##2.1) remove SNVs with mean < 0.25 or mean >0.75
        vaf_cell_matrix <- vaf_cell_matrix %>% filter(mean >0.25 & mean <0.75)
        
        ##2.2) remove SNVs with median < 0.25 or median >0.75
        vaf_cell_matrix <- vaf_cell_matrix %>% filter(median >0.25 & median <0.75)
        
        ##2.2) remove SNVs with >75% nonNA cells w VAF<0.25 OR VAF>0.75 OR VAF between 0.4-0.6
        vaf_cell_matrix <- vaf_cell_matrix %>% filter(nonNA25 < 0.75 & nonNA75 <0.75 & nonNA46 <0.75)
        
        ##filter out SNVs with less than 20 cells
        vaf_cell_matrix <- vaf_cell_matrix %>% filter(nonNA > 20)
        vaf_cell_matrix <- vaf_cell_matrix[, !(names(vaf_cell_matrix) %in% c('nonNA','mean','median','nonNA25','nonNA75', 'nonNA46'))]
        vaf_cell_matrix <- vaf_cell_matrix %>% rownames_to_column()
        names(vaf_cell_matrix)[1] <- 'SNV'
        fwrite(vaf_cell_matrix, paste0(sample_choose,'_cell_type_VAF/', sample_choose, '_', cell_type[m], '_VAF_matrix_barcode_celltype_filtered.txt'), quote = F, sep = '\t', na = 'NA', row.names = F)
      }
    }
    else {
      print(paste0('without processing VAF matrix for ', sample_choose))
    }
    
  dir_work <- getwd()
  setwd(gsub(sample_choose,'',dir_work))
  
}
