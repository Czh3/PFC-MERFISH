# compare cell type with allen data


# load Allen Seurat data: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq
load("../scRNA/Seurat.ss.rda")

ss.seurat$subclass_label
table(ss.seurat$region_label)


genes = intersect(rownames(merFISH_integrated), rownames(ss.seurat))

merFISH_integrated_L2_1 = merFISH_integrated
merFISH_integrated_L2_1 = NormalizeData(merFISH_integrated_L2_1[genes,])
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "L3_cluster", slot = "data")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 


# PL-ILA ACA MOp ALM
ss.seurat1 = subset(ss.seurat, subset = region_label %in% c("PL-ILA", "ACA", "ALM", "MOp"))
ss.seurat1 = ss.seurat1[genes, ]

subclass_label1 = table(ss.seurat1$subclass_label)
ss.seurat1 = subset(ss.seurat1, subclass_label %in% names(subclass_label1[subclass_label1 >= 50]))
ss.seurat1 = NormalizeData(ss.seurat1)
Allen_avg_expr = AverageExpression(ss.seurat1, group.by = c("subclass_label"), slot = "data")
Allen_avg_expr = Allen_avg_expr$RNA
Allen_avg_expr = as.data.frame(Allen_avg_expr) 
Allen_avg_expr = Allen_avg_expr[rowMeans(Allen_avg_expr) > 1, ]

genes = intersect(rownames(merFISH_integrated), rownames(Allen_avg_expr))

merFISH_avg_expr.scale = t(scale(t(merFISH_avg_expr[genes,])))
Allen_avg_expr.scale = t(scale(t(Allen_avg_expr[genes, ])))

Expr_cor = cor(merFISH_avg_expr.scale, Allen_avg_expr.scale[rownames(merFISH_avg_expr.scale), ])

Expr_cor = Expr_cor[, c(3,5,6,10, 7,8,9, 11,12,15,16,17,19,1,2,13,14,20)]
pdf("../figure_new/compareALLEN.cluster.pdf", 8, 8 )
pheatmap::pheatmap(Expr_cor, cluster_rows = F, cluster_cols = F, 
                   color =  colorRampPalette(c("white", "white", "white", "pink", "red3"))(30),
                   #color = colorRampPalette(c("lightblue",  "gray97", "pink", "red", "red4"))(30),
                   border_color = NA)
dev.off()



#### subcluster 
merFISH_integrated_L2_1$subcluster = factor(merFISH_integrated_L2_1$subcluster, levels = levels(Idents(merFISH_integrated_L2_1)))
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "subcluster", slot = "data")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 

ss.seurat1 = subset(ss.seurat, subset = region_label %in% c("PL-ILA", "ACA", "ALM", "MOp"))
ss.seurat1 = ss.seurat1[genes, ]

subclass_label1 = table(ss.seurat1$subclass_label)
ss.seurat1 = subset(ss.seurat1, subclass_label %in% names(subclass_label1[subclass_label1 >= 50]))
ss.seurat1 = NormalizeData(ss.seurat1)
Allen_avg_expr = AverageExpression(ss.seurat1, group.by = c( "subclass_label"), slot = "data")
Allen_avg_expr = Allen_avg_expr$RNA
Allen_avg_expr = as.data.frame(Allen_avg_expr) 
Allen_avg_expr = Allen_avg_expr[rowMeans(Allen_avg_expr) > 1, ]

genes = intersect(rownames(merFISH_integrated), rownames(Allen_avg_expr))

merFISH_avg_expr.scale = t(scale(t(merFISH_avg_expr[genes,])))
Allen_avg_expr.scale = t(scale(t(Allen_avg_expr[genes, ])))

Expr_cor = cor(merFISH_avg_expr.scale, Allen_avg_expr.scale[rownames(merFISH_avg_expr.scale), ])

Expr_cor = Expr_cor[c(1:6, 9:11, 17:18, 7:8, 12:16, 19:52), c(3,5,6,10, 7,8,9, 11,12,15,16,17,19,1,2,13,14,20)]
pdf("../figure_new/compareALLEN.subcluster.pdf", 6, 12 )
pheatmap::pheatmap(Expr_cor, cluster_rows = F, cluster_cols = F, 
                   color =  colorRampPalette(c("white", "white", "white", "pink", "red3"))(30),
                   border_color = NA)
dev.off()

#
merFISH_integrated_L2_1$subcluster = factor(merFISH_integrated_L2_1$subcluster, levels = levels(Idents(merFISH_integrated_L2_1)))
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "subcluster", slot = "data")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 

ss.seurat1 = subset(ss.seurat, subset = region_label %in% c("PL-ILA", "ACA", "ALM", "MOp"))
ss.seurat1 = ss.seurat1[genes, ]

subclass_label1 = table(ss.seurat1$cluster_label)
ss.seurat1 = subset(ss.seurat1, cluster_label %in% names(subclass_label1[subclass_label1 >= 30]))
ss.seurat1 = NormalizeData(ss.seurat1)
Allen_avg_expr = AverageExpression(ss.seurat1, group.by = c( "cluster_label"), slot = "data")
Allen_avg_expr = Allen_avg_expr$RNA
Allen_avg_expr = as.data.frame(Allen_avg_expr) 
Allen_avg_expr = Allen_avg_expr[rowMeans(Allen_avg_expr) > 1, ]

genes = intersect(rownames(merFISH_integrated), rownames(Allen_avg_expr))

merFISH_avg_expr.scale = t(scale(t(merFISH_avg_expr[genes,])))
Allen_avg_expr.scale = t(scale(t(Allen_avg_expr[genes, ])))

Expr_cor = cor(merFISH_avg_expr.scale, Allen_avg_expr.scale[rownames(merFISH_avg_expr.scale), ])

Expr_cor = Expr_cor[c(1:6, 9:11, 17:18, 7:8, 12:16, 19:52), c(18:23, 25:76, 2, 11, 15:17, 24, 3:10, 12, 77, 78, 80, 105:122, 88:104, 82, 83, 87,81,79,85:86)]
Expr_cor[Expr_cor>0.8] <- 0.8
pdf("../figure_new/compareALLEN.all.subcluster.pdf", 16, 12 )
pheatmap::pheatmap(Expr_cor, cluster_rows = F, cluster_cols = F,
                   color =  colorRampPalette(c("white", "white", "white", "white", "white", "white", "pink", "red3"))(50),
                   border_color = NA)
dev.off()


### co-embed
ss.seurat1 = subset(ss.seurat, subset = region_label %in% c("PL-ILA", "ACA", "ALM", "MOp"))

genes_select = intersect(rownames(merFISH_integrated_L3), rownames(ss.seurat))

ss.seurat1 = ss.seurat1[genes_select, ]

subclass_label1 = table(ss.seurat1$subclass_label)
ss.seurat1 = subset(ss.seurat1, subclass_label %in% names(subclass_label1[subclass_label1 >= 50]))
ss.seurat1 = NormalizeData(ss.seurat1)

integrated = merge(ss.seurat1, merFISH_integrated_L3)
integrated@meta.data$orig.ident[is.na(integrated@meta.data$orig.ident)] <- "scRNA"
integrated@meta.data$tech = ifelse(is.na(integrated@meta.data$abs_volume), "scRNA", "MERFISH")

integrated = FindVariableFeatures(integrated) 
integrated = ScaleData(integrated,vars.to.regress = c("nFeature_RNA"), verbose = FALSE)
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
integrated <- RunUMAP(integrated, dims = 1:30, reduction = 'harmony')


p1 = DimPlot(subset(integrated, subset = tech == "scRNA"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "subclass_label", label = F) &
  NoAxes()
p2 = DimPlot(subset(integrated, subset = tech == "MERFISH"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "L3_cluster", label = F) &
  NoAxes()

p3 = DimPlot(integrated, reduction = "umap", group.by = "tech",
             repel = TRUE, shuffle = TRUE) &
  scale_color_manual(values =  c(ggsci::pal_d3()(10))) & NoAxes()

p3|p1|p2
ggsave("../figure_new/Allen_Coembed_scRNA_merFISH_CellType.png", width = 15, height = 5, dpi = 600)



###
merFISH_avg_expr1 = merFISH_avg_expr
colnames(merFISH_avg_expr1) = paste0("PFC_", colnames(merFISH_avg_expr1) )

merFISH_avg_expr1.scale = t(scale(t(merFISH_avg_expr1[genes,])))
Allen_avg_expr.scale = t(scale(t(Allen_avg_expr[genes, ])))

merge_expr = cbind(merFISH_avg_expr1.scale, Allen_avg_expr.scale)

merge_expr.cor = cor(merge_expr)
pheatmap::pheatmap(merge_expr.cor, cluster_rows = T, cluster_cols = T, 
                   color = colorRampPalette(c("lightblue",  "gray97", "yellow", "red", "red4"))(30), border_color = NA)


  
pca = prcomp(t((merge_expr)), scale = T, center = T)
pca = as.data.frame(pca$x)

pca$region = stringr::str_replace(rownames(pca), "_.*", "")
pca$cellType = rownames(pca)

ggplot(pca, aes(PC1, PC2, color = region), label = cellType) +
  geom_point( size=3)+
  theme_cowplot() 


dist_matrix <- as.dist(1-merge_expr)
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
plot(dendrogram)

library(ggdendro)
dendrogram_data <- dendro_data(dendrogram)
dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data


#install.packages("dendextend")
library(dendextend)

L2_color1 = brewer.set1(5)
names(L2_color1) = unique((pca$region))

L2_color1 = L2_color1[ stringr::str_replace(dendrogram_data$labels$label, "_.*", "")]
dend110 <- dendrogram %>% color_labels(col = L2_color1) %>% color_branches(col = L2_color1)

par(mar = c(10, 4, 4, 4))
pdf("../figure_new/cellType_hclust_compareAllen.pdf", 20, 5)
dend110 %>%  plot(main = "Original tree") 
dev.off()



# annotation
row_ann = data.frame(region = stringr::str_replace(colnames(merge_expr.cor), "_.*", ""))
rownames(row_ann) = colnames(merge_expr.cor)

ss.seurat.class_lab = unique(ss.seurat@meta.data[,c("class_label", "region_label", "subclass_label")])
rownames(ss.seurat.class_lab) = paste0(ss.seurat.class_lab$region_label, "_", ss.seurat.class_lab$subclass_label)

row_ann[rownames(ss.seurat.class_lab), "class"] <- ss.seurat.class_lab$class_label 
row_ann[row_ann$region == "PFC", "class"] <- c(rep("Glutamatergic", 18), rep("GABAergic", 19), rep("Non-Neuronal", 15))

color1 = brewer.set1(5)
names(color1) = unique(row_ann$region)

color2 = pals::glasbey()[1:3]
names(color2) = unique(row_ann$class)

anno_colors <- list(region = color1,
                    class = color2
                    )


pdf("../figure_new/Compare_Allen_cellType.heatmap.pdf", height = 20, width = 25)
pheatmap::pheatmap(merge_expr.cor, cluster_rows = T, cluster_cols = T, annotation = row_ann, annotation_colors = anno_colors,
                   clustering_method = "ward.D",
                   color = colorRampPalette(c("lightblue", "white", "gray97", "red", "darkred"))(30), border_color = NA)
dev.off()




# only ext
ext = rownames(row_ann[row_ann$class == "Glutamatergic", ])

merge_expr = cbind(merFISH_avg_expr1.scale, Allen_avg_expr.scale)
merge_expr.cor = cor(merge_expr)

merge_expr.cor.ext = merge_expr.cor[ext, ext]

pdf("../figure_new//Compare_Allen_cellType_onlyExt.heatmap.pdf", height = 10, width = 12)
pheatmap::pheatmap(merge_expr.cor.ext, cluster_rows = T, cluster_cols = T, annotation = row_ann[,1,drop=F], annotation_colors = anno_colors,
                   clustering_method = "ward.D",
                   color = pals::gnuplot(30), border_color = NA)
dev.off()



