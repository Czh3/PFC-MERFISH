# integrate scRNA and merFISH data
# 3/28/2022
library(pals)

# load data
load("../scRNA/all_cell_inDD_sobj.RData")
DimPlot(all_cell_inDD_sobj, reduction = "tsne", group.by= "CellType",
        label = F) + NoAxes() +
scale_color_manual(values =  as.vector(alphabet()))


merFISH_integrated = readRDS(file = "./RDS/merFISH_integrated_L1.seurat.rds")

pdf("../figure_new/Correspondence_scRNA_merFISH_CellType.pdf", width = 3.5, height = 3)
DimPlot(merFISH_integrated, reduction="umap", #label = F,label.size = 6, pt.size = 0.3,
        repel = T,group.by = "L1_cluster", raster=FALSE) +
  scale_color_manual(values =  c(as.vector(alphabet()), as.vector(alphabet2()))) & NoAxes()
dev.off()

#
genes_select = intersect(rownames(merFISH_integrated), rownames(all_cell_inDD_sobj))

merFISH_integrated_avg = AverageExpression(merFISH_integrated, slot = "scale.data", features = genes_select)
Idents(all_cell_inDD_sobj) = "CellType"
scRNA_avg = AverageExpression(all_cell_inDD_sobj, group.by = "CellType", slot = "scale.data", features = genes_select)

MERFISH_scRNA_cor = cor(merFISH_integrated_avg$RNA, scRNA_avg$RNA)
pheatmap::pheatmap(MERFISH_scRNA_cor[,-6], cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c( "gray90", "white", "red2"))(50), border_color = NA)



#
set.seed(3)
merFISH_integrated_L3_sample = merFISH_integrated_L3[,sample(colnames(merFISH_integrated_L3), ncol(all_cell_inDD_sobj))]

merFISH_integrated_L3_sample = subset(merFISH_integrated_L3, orig.ident == "Pain8")

#
genes_select = intersect(rownames(merFISH_integrated_L3_sample), rownames(all_cell_inDD_sobj))

integrated = merge(all_cell_inDD_sobj[genes_select, ], merFISH_integrated_L3_sample[genes_select,])
integrated@meta.data$orig.ident[is.na(integrated@meta.data$orig.ident)] <- "scRNA"
integrated@meta.data$tech = ifelse(is.na(integrated@meta.data$abs_volume), "scRNA", "MERFISH")

integrated = FindVariableFeatures(integrated) 
integrated = ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated = RunPCA(integrated, npcs = 50, verbose = FALSE)

DimPlot(object = integrated, reduction = "pca", pt.size = .1, group.by = "tech")


# harmony
require(harmony)
integrated <- RunHarmony(
  object = integrated,
  group.by.vars = 'tech',
  theta = 3, lambda = 0.05,
  plot_convergence = TRUE
)

# re-compute the UMAP 
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')


p1 = DimPlot(subset(integrated, subset = tech == "scRNA"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "CellType", label = F) &
  NoAxes()
p2 = DimPlot(subset(integrated, subset = tech == "MERFISH"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "CellType", label = F) &
  NoAxes()

p3 = DimPlot(integrated, reduction = "umap", group.by = "tech",
             repel = TRUE, shuffle = TRUE) &
  scale_color_manual(values =  c(ggsci::pal_d3()(10))) & NoAxes()

p3|p1|p2
ggsave("../figure_new/Coembed_scRNA_merFISH_CellType.png", width = 15, height = 5, dpi = 600)



integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)



# Correspondence between merFISH and scRNA
scRNA2merFISH_cluster = function(integrated_obj, cluster){

  integrated_scRNA = subset(integrated_obj, subset = tech == "scRNA")
  integrated_merFISH = subset(integrated_obj, subset = tech != "scRNA")
  
  map_scRNA_merFISH = function(obj_a, obj_b) {
    scRNA2merFISH = data.frame(matrix(nrow = 0, ncol = 3))
    for(i in unique(obj_a[[cluster]][,1])){
      scRNA_cell = colnames(obj_a)[obj_a[[cluster]] == i]
      scRNA_cell_Neighbor = TopNeighbors(integrated_obj@neighbors$RNA.nn, scRNA_cell, n=30)
      
      merFISH_cell = colnames(obj_b)
      scRNA_cell_Neighbor = scRNA_cell_Neighbor[scRNA_cell_Neighbor %in% merFISH_cell]
      
      merFISH_cell_cluster = obj_b[[cluster]][scRNA_cell_Neighbor, ]
      merFISH_cell_cluster =  as.data.frame(table(merFISH_cell_cluster))
      merFISH_cell_cluster$Freq = merFISH_cell_cluster$Freq / sum(merFISH_cell_cluster$Freq)
      merFISH_cell_cluster$scRNA_cluster = i
      scRNA2merFISH = rbind(scRNA2merFISH, merFISH_cell_cluster)
    }
    scRNA2merFISH = reshape2::dcast(scRNA2merFISH, merFISH_cell_cluster~scRNA_cluster, value.var = "Freq")
    rownames(scRNA2merFISH) = scRNA2merFISH$merFISH_cell_cluster
    scRNA2merFISH = scRNA2merFISH[, -1]
    scRNA2merFISH[is.na(scRNA2merFISH)] <- 0
    scRNA2merFISH
  }
  scRNA2merFISH = map_scRNA_merFISH(integrated_scRNA, integrated_merFISH)
  merFISH2scRNA = map_scRNA_merFISH(integrated_merFISH, integrated_scRNA)
  Correspondence_scRNA_merFISH = 0.5 * (merFISH2scRNA + t(scRNA2merFISH)[rownames(merFISH2scRNA), colnames(merFISH2scRNA)])
  #Correspondence_scRNA_merFISH = sqrt(merFISH2scRNA * t(scRNA2merFISH)[rownames(merFISH2scRNA), colnames(merFISH2scRNA)])
  Correspondence_scRNA_merFISH

} 

Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "CellType")

pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)

pdf("../figure_new/Correspondence_scRNA_merFISH_CellType.pdf", width = 3.5, height = 3)
pheatmap::pheatmap(Correspondence_scRNA_merFISH[c(1:5,8,7),-8], cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()

saveRDS(integrated, file = "./RDS/merFISH_scRNA_integrated.seurat.rds")
integrated = readRDS( file = "./RDS/merFISH_scRNA_integrated.seurat.rds")



# Excitatory
set.seed(3)
merFISH_integrated_L3_sample = merFISH_integrated_L3[,sample(colnames(merFISH_integrated_L3), ncol(all_cell_inDD_sobj))]

#merFISH_integrated_L3_sample = subset(merFISH_integrated_L3, orig.ident == "Pain8")

genes_select = intersect(rownames(merFISH_integrated_L3), rownames(all_cell_inDD_sobj))
integrated = merge(all_cell_inDD_sobj[genes_select, ], merFISH_integrated_L3_sample[genes_select,])

integrated_Exc = subset(integrated, cells = colnames(integrated)[integrated$CellType %in% c("Ext", "Excitatory")])
integrated_Exc = integrated_Exc[,!grepl("Inh",integrated_Exc$L1_clusters, fixed=TRUE)]

integrated_Exc = integrated_Exc[, colSums(integrated_Exc@assays$RNA@counts) > 50]
integrated_Exc = NormalizeData(integrated_Exc)
integrated_Exc = FindVariableFeatures(integrated_Exc) 
integrated_Exc$abs_volume[is.na(integrated_Exc$abs_volume)] = mean(integrated_Exc$abs_volume, na.rm=T)
integrated_Exc$orig.ident[is.na(integrated_Exc$orig.ident)] = "scRNA"
integrated_Exc = ScaleData(integrated_Exc, vars.to.regress = c("orig.ident", "nCount_RNA"), verbose = FALSE)
integrated_Exc = RunPCA(integrated_Exc, npcs = 50,  verbose = FALSE)

integrated_Exc$tech = ifelse(is.na(integrated_Exc$doulet_scores), "scRNA", "MERFISH")
DimPlot(object = integrated_Exc, reduction = "pca", pt.size = .1, group.by = "tech") 

# harmony
require(harmony)
integrated_Exc <- RunHarmony(
  object = integrated_Exc,
  group.by.vars = c('tech'),
  #theta = 3, lambda = 0.05,
  plot_convergence = TRUE,
) # regress batch effect & tech

DimPlot(object = integrated_Exc, reduction = "harmony", pt.size = .1, group.by = "tech") 

integrated_Exc <- FindNeighbors(integrated_Exc, reduction = "harmony", dims = 1:30, 
                                k.param = 30, return.neighbor = T)

integrated_Exc <- RunUMAP(integrated_Exc, reduction = 'harmony', dims = 1:30)

p1 = DimPlot(subset(integrated_Exc, subset = tech == "scRNA"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "L1_clusters", label = F) &
    NoAxes()
p2 = DimPlot(subset(integrated_Exc, subset = tech == "MERFISH"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "subcluster", label = F) &
   NoAxes()

p3 = DimPlot(integrated_Exc, reduction = "umap", group.by = "tech",
        repel = TRUE, shuffle = TRUE) &
  scale_color_manual(values =  c(ggsci::pal_d3()(10))) & NoAxes()

p3|p1|p2
ggsave("../figure_new/Coembed_scRNA_merFISH_Ect.png", width = 15, height = 5, dpi = 600)



integrated_Exc$cluster = na.omit(c(integrated_Exc$L1_clusters, integrated_Exc$subcluster))

Correspondence_scRNA_merFISH_subCluster = scRNA2merFISH_cluster(integrated_Exc, "cluster")

Correspondence_scRNA_merFISH_subCluster[Correspondence_scRNA_merFISH_subCluster > 0.8] <- 0.8

pheatmap::pheatmap(Correspondence_scRNA_merFISH_subCluster, cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c( "gray97", "red", "darkred"))(30), border_color = NA)

# re-order
pdf("../figure_new/Correspondence_scRNA_merFISH_Ect.pdf", width = 4, height = 3)
pheatmap::pheatmap(Correspondence_scRNA_merFISH_subCluster[paste0("Exc_", c(4,3,2,1,6,10,7,9,8,11,12,13,5)),
                                                           c(1,3,2,4:18)],
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("white", "gray90", "red", "red4"))(30), border_color = NA)
dev.off()

saveRDS(integrated_Exc, file = "./RDS/merFISH_scRNA_integrated_Exc_0719.seurat.rds")
integrated_Exc = readRDS(file = "./RDS/merFISH_scRNA_integrated_Exc.seurat.rds")


FeaturePlot(integrated_Exc, features = "nCount_RNA", reduction="umap") 


# Inhibitory
integrated_Inh = subset(integrated, cells = colnames(integrated)[integrated$CellType %in% c("Inh", "Inhibitory")])
integrated_Inh = integrated_Inh[,!grepl("Exc",integrated_Inh$L1_clusters, fixed=TRUE)]


integrated_Inh = NormalizeData(integrated_Inh)
integrated_Inh = FindVariableFeatures(integrated_Inh) 
integrated_Inh$orig.ident[is.na(integrated_Inh$orig.ident)] = "scRNA"
integrated_Inh = ScaleData(integrated_Inh, vars.to.regress = c("orig.ident", "nCount_RNA"), verbose = FALSE)
integrated_Inh = RunPCA(integrated_Inh, npcs = 50, features = genes_select, verbose = FALSE)

integrated_Inh$tech = ifelse(is.na(integrated_Inh$doulet_scores), "scRNA", "MERFISH")
DimPlot(object = integrated_Inh, reduction = "pca", pt.size = .1, group.by = "L1_clusters", split.by="tech") 


# harmony
require(harmony)
integrated_Inh <- RunHarmony(
  object = integrated_Inh,
  group.by.vars = c('tech'),
  #theta = 3, lambda = 0.05
  plot_convergence = TRUE
  
) # regress batch effect & tech

DimPlot(object = integrated_Inh, reduction = "harmony", pt.size = .1, group.by = "tech") 
ElbowPlot(integrated_Inh, reduction = "harmony", ndims=50)

integrated_Inh <- FindNeighbors(integrated_Inh, reduction = "harmony", dims = 1:30,
                                k.param = 30, return.neighbor = T)
integrated_Inh <- RunUMAP(integrated_Inh, reduction = 'harmony', dims = 1:20)
p1 = DimPlot(subset(integrated_Inh, subset = tech == "scRNA"), reduction = "umap",
             shuffle = TRUE, group.by = "L1_clusters", label = F) &
  NoAxes()
p2 = DimPlot(subset(integrated_Inh, subset = tech == "MERFISH"), reduction = "umap",
             shuffle = TRUE, group.by = "subcluster", label = F) &
  NoAxes()
p1|p2

p3 = DimPlot(integrated_Inh, reduction = "umap", group.by = "tech", pt.size = 0.5,  shuffle = TRUE) &
  scale_color_manual(values =  c(ggsci::pal_d3()(10))) & NoAxes()

p3|p1|p2
ggsave("../figure_new/Coembed_scRNA_merFISH_Inh.png", width = 15, height = 5, dpi = 600)


integrated_Inh$cluster = na.omit(c(integrated_Inh$L1_clusters, integrated_Inh$subcluster))
Correspondence_scRNA_merFISH_subCluster = scRNA2merFISH_cluster(integrated_Inh, "cluster")

Correspondence_scRNA_merFISH_subCluster[Correspondence_scRNA_merFISH_subCluster > 0.8] <- 0.8
pheatmap::pheatmap(Correspondence_scRNA_merFISH_subCluster,
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)

pdf("../figure_new//Correspondence_scRNA_merFISH_Inh.pdf", width = 4, height = 3)
pheatmap::pheatmap(Correspondence_scRNA_merFISH_subCluster[paste0("Inhib_", c(7,5,6,10,11,12,3,4,9,2,1,8)),
                                                           c(1,2,5:7,9,3,8,10,11,17:19, 13,15,16,12,14)],
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()

saveRDS(integrated_Inh, file = "./RDS/merFISH_scRNA_integrated_Inh.seurat.rds")
integrated_Inh = readRDS(file = "./RDS/merFISH_scRNA_integrated_Inh.seurat.rds")




FeaturePlot(integrated_Exc, features = c("Tshz2","Pou3f1", "Syt6"), 
            slot = "scale.data",
            reduction="umap", split.by = "tech", min.cutoff = 1)  & 
  NoAxes()  &
  scale_color_gradientn(colours = c("grey95", "pink", "red", "darkred")) 


FeaturePlot(integrated_Inh, features = c("Cnr1","Pvalb", "Satb1"), 
            slot = "scale.data", pt.size = 0.5,
            reduction="umap", split.by = "tech", min.cutoff = 0.5)  & 
  NoAxes()  &
  scale_color_gradientn(colours = c("grey95", "pink", "red", "darkred")) 











#####
# expr cor
load("../scRNA/all_cell_inDD_sobj.RData")
genes_select = intersect(rownames(merFISH_integrated_L3), rownames(all_cell_inDD_sobj))

#scRNA = all_cell_inDD_sobj[genes_select, colnames(all_cell_inDD_sobj)[grepl("Exc",all_cell_inDD_sobj$L1_clusters)]]
#merFISH = merFISH_integrated_L2[genes_select, colnames(merFISH_integrated_L2)[merFISH_integrated_L2$CellType == "Ext"]]
scRNA = all_cell_inDD_sobj[genes_select, ]
merFISH = merFISH_integrated_L3[genes_select, sample(colnames(merFISH_integrated_L3), 30000)]

scRNA = subset(scRNA, CellType == "Excitatory")
merFISH = subset(merFISH, CellType == "Ext")

scRNA = NormalizeData(scRNA)
VariableFeatures(scRNA) = genes_select
scRNA = ScaleData(scRNA)

merFISH = NormalizeData(merFISH)
merFISH = ScaleData(merFISH, vars.to.regress = c("orig.ident", "abs_volume"))

integrated_Exc.scRNA = AverageExpression(scRNA, group.by = "L1_clusters", slot = "data")$RNA
integrated_Exc.merFISH = AverageExpression(merFISH, group.by = "L3_cluster", slot = "data")$RNA
integrated_Exc.scRNA = t(scale(t(log(integrated_Exc.scRNA))))
integrated_Exc.merFISH = t(scale(t(log(integrated_Exc.merFISH))))

integrated_Exc.scRNA = na.omit(integrated_Exc.scRNA)

Ext_cor = cor(integrated_Exc.scRNA, integrated_Exc.merFISH[rownames(integrated_Exc.scRNA), ])
Ext_cor[Ext_cor < 0] <- 0

pheatmap::pheatmap(Ext_cor, cluster_rows = T, cluster_cols = T, 
                   color = colorRampPalette(c("white", "gray97", "red", "darkred"))(30), border_color = NA)

pheatmap::pheatmap(Ext_cor[c(8,7,1,10,12,11,9,2,4,5,6,3), c(1, 10, 9, 11:16, 3:8,2)], cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("white", "gray97", "red", "darkred"))(30), border_color = NA)


FeaturePlot(all_cell_inDD_sobj, features = c("Ddit4l", "Thsd7a"))



