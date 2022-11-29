# merFISH pipeline
library(Seurat)
library(ggplot2)
library(data.table)
library(SeuratData)
library(dplyr)
library(DoubletFinder)

merFISH_pipeline = function(dir, proj){
  
  # generate expr matrix
  expr1 = fread(paste0(dir, "/countsPerCellCorrectedIn.csv"),  header = F, data.table = FALSE)
  expr2 = fread(paste0(dir, "/countsPerCellExactIn.csv"),  header = F, data.table = FALSE)
  expr = expr1 + expr2

  coln = fread(paste0(dir, "/featureNames.csv"), header = F)
  rown = fread(paste0(dir, "/geneNames.csv"), header = F)
  
  colnames(expr) = coln$V1
  rownames(expr) = rown$V1
  
  meta_data = read.csv(paste0(dir, "/feature_metadata.csv"), header = T, row.names = 1)
  
  # convert to seurat
  merFISH <- CreateSeuratObject(counts = expr, meta.data = meta_data[,1:6], project = proj, assay = "RNA")
  
  # remove Blank "Genes"
  merFISH <- subset(merFISH, features = rown$V1[1:416])
  
  # filter cells by cell volume
  merFISH <- subset(merFISH, subset = abs_volume >= 100 & abs_volume <= 4000)
  
  # remove cells with lower/higher RNA counts
  count_quantile = quantile(merFISH$nCount_RNA, probs= 0.98)
  merFISH <- subset(merFISH, subset = nCount_RNA >= 10 & nFeature_RNA >= 10 & nCount_RNA <= count_quantile)
  
  # add spatial infor
  coord.df = merFISH@meta.data[,c(7,8)]
  merFISH@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
  
  # normalization
  # First normalize cell size
  norm_data = t(apply(merFISH@assays$RNA@counts, 1, function(x) x / merFISH@meta.data$abs_volume)) * mean(merFISH@meta.data$abs_volume) 
  # Normalize the total RNA counts
  norm_data = t(apply(norm_data, 1, function(x) x / colSums(norm_data) * 500))
  
  merFISH@assays$RNA@data = as(log(norm_data + 1), "dgCMatrix")
  
  # Identify highly variable genes
  merFISH <- FindVariableFeatures(merFISH, selection.method = "vst")
  
  all.genes <- rownames(merFISH)
  merFISH <- ScaleData(merFISH, features = all.genes)
  
  # PCA
  merFISH <- RunPCA(merFISH)
  
  # DoubletFinder
  #merFISH <- doubletFinder_v3(merFISH, pN = 0.25, pK = 0.09, nExp = ncol(merFISH)*0.12, PCs = 1:30, sct = T)
  #merFISH@meta.data[10] <- NULL
  #colnames(merFISH@meta.data)[10] <- "DoubletFinder"
  
  # Scrublet
  scrublet = reticulate::import("scrublet")
  
  scr <- scrublet$Scrublet(counts_matrix = t(merFISH@assays$RNA@counts),
                           expected_doublet_rate = 0.12)
  scr_out = scr$scrub_doublets()
  
  merFISH$doulet_scores = scr_out[[1]]
  merFISH$predicted_doulets = scr_out[[2]]
  
  return(merFISH)
}


