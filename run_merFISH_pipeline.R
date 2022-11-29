# run merFISH basic analysis

setwd("/nfs4/chaozhang/proj/Neuron/Aritra/merFISH/script")

# init
source("merFISH_pipeline.R")
source("functions.R")

update_geom_defaults("point", list(shape = 16))


### run pipeline per sample
# 
merFISH_mouse0.1 = merFISH_pipeline("../data/20210911_PFCL1", "mouse0.1")
merFISH_mouse0.2 = merFISH_pipeline("../data/20211004_PFCL2", "mouse0.2")

#
merFISH_mouse1.1 = merFISH_pipeline("../data/20211208_PFCL5", "mouse1.1")
merFISH_mouse1.2 = merFISH_pipeline("../data/20211208_PFCL6", "mouse1.2")

#
merFISH_mouse2.1 = merFISH_pipeline("../data/20211024_PFCL1_Results_Pain-d1", "mouse2.1")
merFISH_mouse2.2 = merFISH_pipeline("../data/20211101_PFCL1_Results_Pain-d2", "mouse2.2")

#
merFISH_mouse3.1 = merFISH_pipeline("../data/20220103_PFCL7", "mouse3.1")
merFISH_mouse3.2 = merFISH_pipeline("../data/20220106_PFCL1_Results", "mouse3.2")
merFISH_mouse3.3 = merFISH_pipeline("../data/20220107_PFCL1_Results", "mouse3.3")
merFISH_mouse3.4 = merFISH_pipeline("../data/20220128_PFCL1_Results", "mouse3.4")

#
merFISH_pain1 = merFISH_pipeline("../data/20220218_PFCL1_Results_Pain", "Pain1")
merFISH_pain2 = merFISH_pipeline("../data/20220210_PFCL1_Results_Pain", "Pain2")

#
merFISH_pain3 = merFISH_pipeline("../data/20220227_PFCL1_Results_Cingulate and NAc", "Pain3")

#
merFISH_pain4 = merFISH_pipeline("../data/20220305_RESULTS_Pain and Control_Cingulate and NAc", "Pain4")

merFISH_pain5 = merFISH_pipeline("../data/20220416_PFCL1_Results", "Pain5")

merFISH_pain6 = merFISH_pipeline("../data/20220429_PFCL1_Results", "Pain6")

merFISH_pain7 = merFISH_pipeline("../data/20220430_PFCL1_Results", "Pain7")

merFISH_pain8 = merFISH_pipeline("../data/20220513_PFCL1_Results", "Pain8")

merFISH_pain9 = merFISH_pipeline("../data/20220514_PFCL1_Results", "Pain9")


### integrated
merFISH_integrated = merge(merFISH_mouse0.1, c(merFISH_mouse0.2,
                                               merFISH_mouse1.1, merFISH_mouse1.2,
                                               merFISH_mouse2.1, merFISH_mouse2.2,
                                               merFISH_mouse3.1, merFISH_mouse3.2, merFISH_mouse3.3, merFISH_mouse3.4,
                                               merFISH_pain1, merFISH_pain2, merFISH_pain3, merFISH_pain4, merFISH_pain5,
                                               merFISH_pain6, merFISH_pain7, merFISH_pain8, merFISH_pain9))

saveRDS(merFISH_integrated, file = "./RDS/merFISH_integrated.seurat0627.rds")
rm(merFISH_mouse0.2,
     merFISH_mouse1.1, merFISH_mouse1.2,
     merFISH_mouse2.1, merFISH_mouse2.2,
     merFISH_mouse3.1, merFISH_mouse3.2, merFISH_mouse3.3, merFISH_mouse3.4,
     merFISH_pain1, merFISH_pain2, merFISH_pain3, merFISH_pain4, merFISH_pain5,
     merFISH_pain6, merFISH_pain7, merFISH_pain8, merFISH_pain9)
gc()

merFISH_integrated = readRDS(file = "./RDS/merFISH_integrated.seurat0627.rds")

# remove low quality 2 samples Pain6 Pain7 
merFISH_integrated = subset(merFISH_integrated, subset = orig.ident %in% c("Pain6", "Pain7"), invert = T)

# remove cells from striatum region
cell_remove = gate_points(merFISH_integrated, image = "image.12")
cell_remove1 = gate_points(merFISH_integrated, image = "image.13")

merFISH_integrated = merFISH_integrated[, !colnames(merFISH_integrated) %in% c(cell_remove, cell_remove1)] 

SpatialFeaturePlot(merFISH_integrated, features = "nCount_RNA", images = "image.12")

# scale data
all.genes = rownames(merFISH_integrated)
merFISH_integrated = FindVariableFeatures(merFISH_integrated, selection.method = "vst")
merFISH_integrated = ScaleData(merFISH_integrated, vars.to.regress = c("abs_volume"))
merFISH_integrated = RunPCA(merFISH_integrated, npcs = 50, verbose = FALSE)
ElbowPlot(merFISH_integrated, ndims=50)

# determine how many PCs
determinePCA(merFISH_integrated) # 39

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated <- RunHarmony(
        object = merFISH_integrated,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        dims.use = 1:39,
        plot_convergence = TRUE
)

library(pals)
DimPlot(object = merFISH_integrated, reduction = "pca", 
        pt.size = .1,
        raster=FALSE, shuffle = T,
        group.by = "orig.ident") &
        scale_color_manual(values =  as.vector(alphabet()))


DimPlot(object = merFISH_integrated, reduction = "harmony", 
        pt.size = .1,
        raster=FALSE, shuffle = T,
        group.by = "orig.ident") &
        scale_color_manual(values =  as.vector(alphabet()))


merFISH_integrated <- RunUMAP(merFISH_integrated, dims = 1:39, reduction = 'harmony')

DimPlot(object = merFISH_integrated,
        pt.size = .01,
        raster=FALSE,shuffle = T,
        split.by = "predicted_doulets") &
        scale_color_manual(values =  as.vector(alphabet())) & NoAxes()

FeaturePlot(object = merFISH_integrated,
            features = c("nCount_RNA", "abs_volume"),
            pt.size = .01,
            raster=FALSE) &
        NoAxes()  &
        scale_color_gradientn(colours = rainbow(10)[3:10]) 
        scale_color_gradientn(colours = c( "grey95", "blue"))

FeaturePlot(merFISH_integrated, features = "centroid_3", raster=FALSE) &
        scale_color_gradientn(colours = rainbow(10)) & NoAxes() 

abs_volume_outliar = quantile(merFISH_integrated$abs_volume, probs = c(0.98))
FeaturePlot(subset(merFISH_integrated, subset = abs_volume > abs_volume_outliar[[1]]),
            raster=FALSE, features = "abs_volume") &
        scale_color_gradientn(colours = rainbow(10)) & NoAxes() 

VlnPlot(object = merFISH_integrated, features = c('nCount_RNA', "nFeature_RNA"), pt.size = -1)

merFISH_integrated = subset(merFISH_integrated, subset = abs_volume < abs_volume_outliar[[1]] & predicted_doulets == FALSE) # 490733


#
# clustering
merFISH_integrated <- FindNeighbors(merFISH_integrated, reduction = "harmony", dims = 1:39, k.param = 10)

merFISH_integrated <- FindClusters(merFISH_integrated, 
                                   resolution = 0.3, 
                                   method = "igraph",
                                   algorithm = 4,
                                   group.singletons = FALSE,
                                   verbose = FALSE)
 
saveRDS(merFISH_integrated, file = "./RDS/merFISH_integrated.seurat.rds")
merFISH_integrated = readRDS(file = "./RDS/merFISH_integrated.seurat.rds")


# remove cell cluster that contains less than 10 cells
cluster_select = table(merFISH_integrated$seurat_clusters)
cluster_select = cluster_select[cluster_select > 10]


DimPlot(object = subset(merFISH_integrated, subset = seurat_clusters %in% names(cluster_select)),
        pt.size = .01,
        label = T) &
        scale_color_manual(values =  as.vector(polychrome())) & NoAxes()


VlnPlot(subset(merFISH_integrated, subset = seurat_clusters %in% names(cluster_select)), 
        features = c("doulet_scores", "nCount_RNA", "abs_volume"), pt.size = -1)

FeaturePlot(merFISH_integrated, features = c("centroid_3"), 
            raster=FALSE, pt.size = 0.1) &
        NoAxes() &
        scale_color_gradientn(colours = rainbow(10)) 

FeaturePlot(merFISH_integrated, features = c("Drd1", "Drd2"), 
            pt.size = 0.1) &
        NoAxes()  &
        scale_color_gradientn(colours = c( "grey95", "red", "red4"))


SpatialPlot(subset(merFISH_integrated, idents = c(14)), 
            pt.size.factor=1,crop = T,
            images = c("image.4", "image.5", "image.7"),
            stroke=0)& theme(aspect.ratio = 1)


# remove cluster 4, which are MSN cells
cluster_select = cluster_select[names(cluster_select) != "13"]
# remove cluster 14, no specific marker, maybe incomplete cells
cluster_select = cluster_select[names(cluster_select) != "14"]
# remove singleton
cluster_select = cluster_select[names(cluster_select) != "singleton"]

table(merFISH_integrated$orig.ident, merFISH_integrated$seurat_clusters)[, names(cluster_select)]

merFISH_integrated = subset(merFISH_integrated, idents = names(cluster_select))


DimPlot(object = merFISH_integrated,
        pt.size = .01,
        label = T) &
        scale_color_manual(values =  as.vector(polychrome())) & NoAxes()


FeaturePlot(merFISH_integrated, features = c("Gjb6", "Cldn5", "C1qb", "Mog", "Pdgfra", "Serpinf1", "Gad1", "Grin2a", "Tshz2"), 
            pt.size = 0.1) &
        NoAxes()  &
        scale_color_gradientn(colours = c( "grey95", "red", "red4"))


L1_cluster <- c("Excitatory","Inhibitory", "OPC", "Excitatory", "Excitatory", "Excitatory", "Endo", "Astrocyte", "Excitatory", "Inhibitory", "Oligo",
                "Microglia")

names(L1_cluster) <- levels(merFISH_integrated)
merFISH_integrated <- RenameIdents(merFISH_integrated, L1_cluster)
merFISH_integrated@meta.data$L1_cluster = Idents(merFISH_integrated)

DimPlot(object = merFISH_integrated,
        pt.size = .01,
        raster=FALSE,
        label = T) &
        scale_color_manual(values =  as.vector(alphabet2())) & NoAxes()

# rerun umap
#merFISH_integrated = ScaleData(merFISH_integrated, vars.to.regress = "abs_volume") 
#merFISH_integrated = RunPCA(merFISH_integrated, npcs = 50, verbose = FALSE)

#merFISH_integrated <- RunHarmony(
#        object = merFISH_integrated,
#        group.by.vars = 'orig.ident',
#        assay.use = "RNA",
#        dims.use = 1:39,
#        plot_convergence = TRUE
#)


merFISH_integrated <- RunUMAP(merFISH_integrated, dims = 1:39, n.neighbors = 10L,
                              reduction = 'harmony')


Idents(merFISH_integrated) <- factor(x = Idents(merFISH_integrated), levels = sort(levels(merFISH_integrated)))
merFISH_integrated@meta.data$L1_cluster = Idents(merFISH_integrated)

L1_color = c("#AA0DFE","#993F00",
             pals::kelly()[c(5:8, 10:14)],
             "maroon2", "yellow3")

DimPlot(object = merFISH_integrated,
        pt.size = .01,
        label.size = 5,
        raster=FALSE,
        label = T) &
        scale_color_manual(values = L1_color) & NoAxes()

FeaturePlot(merFISH_integrated,
            features = c("nCount_RNA", "abs_volume", "centroid_3", "doulet_scores"), 
            pt.size = 0.1) &
        NoAxes()  &
        scale_color_gradientn(colours = rainbow(10)[4:10])

SpatialPlot(merFISH_integrated, 
            pt.size.factor=1.5,
            crop = T,
            images = "image.3",
            stroke=0) &  theme(aspect.ratio = 2) &
        scale_fill_manual(values = L1_color) & NoAxes() &
        theme(legend.position = "None")

DimPlot(object = merFISH_integrated, group.by = "orig.ident",
        pt.size = .01,
        label.size = 5,
        raster=FALSE,
        shuffle=T,
        label = F) &
        scale_color_manual(values = as.vector(pals::alphabet())) & NoAxes()

FeaturePlot(merFISH_integrated, features = c("nCount_RNA", "nFeature_RNA"),
            pt.size = .01,
            raster=FALSE) & NoAxes() &
        scale_color_gradientn(colours = c("grey95", "blue"))

saveRDS(merFISH_integrated, file = "./RDS/merFISH_integrated_L1.seurat.rds")
merFISH_integrated = readRDS("./RDS/merFISH_integrated_L1.seurat.rds")

#
merFISH_integrated.markers = FindAllMarkers(merFISH_integrated,
                                                logfc.threshold = 1,
                                                min.pct = 0.5,
                                                only.pos = T)

merFISH_integrated.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated, features = c("Gjb6", "Acsbg1", "Cldn5", "Flt1",
                                                  "Grin2a", "Grin2b", "Gad1", "Gad2" ,
                                                  "Ctss", "C1qb", "Mog", "Ermn", "Pdgfra", "Cacng4"),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = pals::coolwarm(10))




# all cell
plot_cell_location(merFISH_integrated, "L3_cluster", color = L3_color , "../data/20211208_PFCL6", 1:500)

# ext cell
plot_cell_location(merFISH_integrated,
                   "L1_cluster", color = c(rep("gray90", 2), L1_color[3:7], rep("gray90", 10)),
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:500)
# inh cell
plot_cell_location(merFISH_integrated,
                   "L1_cluster", color = c(rep("gray90", 7), L1_color[8:9], rep("gray90", 10)),
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:500)
# non neuron
plot_cell_location(merFISH_integrated,
                   "L1_cluster", color = c(L1_color[1:2], rep("gray90", 7), L1_color[10:13]),
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:500)


SpatialDimPlot(merFISH_integrated,
               cells.highlight = CellsByIdentities(object = merFISH_integrated), 
               images = "image.3",
               stroke=0,
               pt.size.factor = 2,
               ncol = 7,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) & theme(aspect.ratio = 1.5) & coord_flip() & scale_y_reverse()



######### sub types
# ext
merFISH_integrated_ext = subset(merFISH_integrated, idents = "Excitatory") # 256190 samples 

nCount_RNA_quantile = quantile(merFISH_integrated_ext$nCount_RNA, probs = c(0.02, 0.98))
merFISH_integrated_ext = subset(merFISH_integrated_ext, 
                                subset = nCount_RNA >= 50 & nCount_RNA < nCount_RNA_quantile[2]) # 232157


merFISH_integrated_ext = ScaleData(merFISH_integrated_ext, vars.to.regress = c("abs_volume")) 
merFISH_integrated_ext = RunPCA(merFISH_integrated_ext, npcs = 50, verbose = FALSE)
ElbowPlot(merFISH_integrated_ext, ndims=50)
determinePCA(merFISH_integrated_ext) # 34

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_ext <- RunHarmony(
        object = merFISH_integrated_ext,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        dims.use = 1:34,
        theta = 3,
        lambda = 0.05,
        plot_convergence = TRUE
)

merFISH_integrated_ext <- RunUMAP(merFISH_integrated_ext, dims = 1:34,
                              reduction = 'harmony')

DimPlot(merFISH_integrated_ext, group.by = "orig.ident")


FeaturePlot(merFISH_integrated_ext, features = "nCount_RNA",
            pt.size = .01,
            raster=FALSE) & NoAxes() &
        scale_color_gradientn(colours = c("blue", "grey95", "red"))

FeaturePlot(merFISH_integrated_ext, features = "abs_volume") &
        scale_color_gradientn(colours = rainbow(10)) & NoAxes() 


#
KNN_K = determineK(merFISH_integrated_ext, 2, 4, reduction="harmony", dims=1:34)
ggsave("../figure_new//merFISH_integrated_ext_determineK.pdf", width = 4, height = 4)

# clust
merFISH_integrated_ext <- FindNeighbors(merFISH_integrated_ext, k.param = 50,
                                        reduction = "harmony", dims = 1:34)

merFISH_integrated_ext <- FindClusters(merFISH_integrated_ext, 
                                       resolution = 2,
                                       method = "igraph",
                                       algorithm = 4,
                                       group.singletons = F,
                                       verbose = FALSE)


# remove cell cluster that contains less than 10 cells
cluster_select = table(merFISH_integrated_ext$RNA_snn_res.2)
cluster_select = cluster_select[cluster_select > 10]
# remove singleton
cluster_select = cluster_select[names(cluster_select) != "singleton"]

merFISH_integrated_ext = subset(merFISH_integrated_ext, subset = RNA_snn_res.2 %in% names(cluster_select)) #226310

DimPlot(merFISH_integrated_ext, label = T, raster=FALSE)

# remove doublets, check the makers
merFISH_integrated_ext.markers = FindAllMarkers(merFISH_integrated_ext,
                                                logfc.threshold = 1, 
                                                min.pct = 0.2,
                                                only.pos = T)
merFISH_integrated_ext.markers_top <- 
        merFISH_integrated_ext.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

FeaturePlot(merFISH_integrated_ext, features = "Mog") 
FeaturePlot(merFISH_integrated_ext, features = c("Gad1","Gad2"))
FeaturePlot(merFISH_integrated_ext, features = c("C1qb")) 
FeaturePlot(merFISH_integrated_ext, features = c("Cldn5"))
FeaturePlot(merFISH_integrated, features = c("Moxd1", "Kcnab1", "Clic5")) 

FeaturePlot(merFISH_integrated_ext, features = c("Lhx6", "Ndst4")) &
        scale_color_gradientn(colours = rainbow(10))

SpatialPlot(subset(merFISH_integrated_ext, idents = c(5,6,17)), 
            images = "image.4")

#
FindMarkers(merFISH_integrated_ext,
            ident.1 = 19,
               only.pos = T)

table(merFISH_integrated_ext$orig.ident, merFISH_integrated_ext$seurat_clusters)[, names(cluster_select)]


merFISH_integrated_ext = subset(merFISH_integrated_ext, idents = c(19:24,26:33), invert = T) # 193740 samples

KNN_K = determineK(merFISH_integrated_ext, 2, 4, reduction="harmony", dims=1:34)
ggsave("../figure_new//merFISH_integrated_ext_determineK.pdf", width = 4, height = 4)

saveRDS(merFISH_integrated_ext, "./RDS/merFISH_integrated_ext_0708.rds")


# re-clusting
#merFISH_integrated_ext = ScaleData(merFISH_integrated_ext, vars.to.regress = c("abs_volume")) 
#merFISH_integrated_ext = RunPCA(merFISH_integrated_ext, npcs = 30, verbose = FALSE)
#merFISH_integrated_ext <- RunHarmony(
#        object = merFISH_integrated_ext,
#        group.by.vars = 'orig.ident',
#        assay.use = "RNA",
#        plot_convergence = TRUE
#)
merFISH_integrated_ext <- RunUMAP(merFISH_integrated_ext, dims = 1:34, n.neighbors = 50, #min.dist = 0.1,
                                  reduction = 'harmony')


Idents(merFISH_integrated_ext) = "RNA_snn_res.2"            

DimPlot(merFISH_integrated_ext, label = T, raster=FALSE) 
        scale_color_manual(values = as.vector(pals::alphabet(18))) & NoAxes()


#source("stable_cluster.R")
#stable_c = stable_cluster(merFISH_integrated_ext, 2, "RNA_snn_res.2", 1, reduction="harmony", dims = 1:30)
#stable_c = one_bootstrap(merFISH_integrated_ext, 2, "RNA_snn_res.2", 1,reduction="harmony", dims=1:30)

#merFISH_integrated_ext = subset(merFISH_integrated_ext, idents = na.omit(stable_c[stable_c$percent >= 0.5, 1]))  # 136331


# re-label
Ext_subcluster <- c(1:length(levels(merFISH_integrated_ext)))
names(Ext_subcluster) <- levels(merFISH_integrated_ext)
merFISH_integrated_ext <- RenameIdents(merFISH_integrated_ext, Ext_subcluster)
merFISH_integrated_ext@meta.data$Ext_subcluster = Idents(merFISH_integrated_ext)


DimPlot(merFISH_integrated_ext,
        group.by = "Ext_subcluster",
        pt.size = .01,
        label = T,
        label.size = 5,
        raster=FALSE) &
        scale_color_manual(values = as.vector(pals::polychrome())[-c(1:2)]) & NoAxes()


saveRDS(merFISH_integrated_ext, "./RDS/merFISH_integrated_ext_0711.rds")
merFISH_integrated_ext = readRDS( "./RDS/merFISH_integrated_ext_0711.rds")

SpatialDimPlot(merFISH_integrated_ext,
               images = "image.3",
               crop = T,
               stroke=0.1,
               pt.size.factor = 1.6) &
        theme(aspect.ratio = 2) &
        scale_fill_manual(values = as.vector(pals::polychrome())[-c(1:2)]) & NoAxes()



SpatialDimPlot(merFISH_integrated_ext,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_ext), 
               images = "image.3",
               crop = T,
               stroke=0,
               pt.size.factor = 2,
               ncol = 5,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 2)


merFISH_integrated$Ext_subcluster = as.character(merFISH_integrated_ext$Ext_subcluster[colnames(merFISH_integrated)])
merFISH_integrated$Ext_subcluster[is.na(merFISH_integrated$Ext_subcluster)] <- "Other"
merFISH_integrated$Ext_subcluster = factor(merFISH_integrated$Ext_subcluster, levels = c(levels(merFISH_integrated_ext$Ext_subcluster), "Other"))
plot_cell_location(merFISH_integrated, "Ext_subcluster", color = c(as.vector(pals::polychrome())[3:21], "gray90") , "../data/20211208_PFCL6/feature_metadata.csv", 1:500)


merFISH_integrated_ext.markers = FindAllMarkers(merFISH_integrated_ext,
                                                logfc.threshold = 0.5,
                                                min.pct = 0.3,
                                                only.pos = T)
write.csv(merFISH_integrated_ext.markers, "./DEGs/merFISH_integrated_ext.markers0711.csv", quote = F)
merFISH_integrated_ext.markers$pct_diff = merFISH_integrated_ext.markers$pct.1 - merFISH_integrated_ext.markers$pct.2

merFISH_integrated_ext.markers_top <- 
        merFISH_integrated_ext.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_ext, features = unique(merFISH_integrated_ext.markers_top$gene)) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c(pals::coolwarm(10)[3:10]))

FeaturePlot(merFISH_integrated_ext, features = c("Otof", "Cux2", "Rorb", "Hrh3", "Nnat", "Fezf2", "Tshz2", "Syt6"),
            min.cutoff = 0.5, ncol = 4, pt.size = 1.5) &
        scale_color_gradientn(colours = c( "grey95", "pink", "red", "red4")) & NoAxes()



##
# inh
merFISH_integrated_inh = subset(merFISH_integrated, idents = levels(merFISH_integrated)[4]) # 49857 samples 

merFISH_integrated_inh = ScaleData(merFISH_integrated_inh, vars.to.regress = c("abs_volume")) # regress out abs_volume
merFISH_integrated_inh = RunPCA(merFISH_integrated_inh, npcs = 50, verbose = FALSE)
ElbowPlot(merFISH_integrated_inh, ndims=50)
determinePCA(merFISH_integrated_inh) # 31

DimPlot(merFISH_integrated_inh)

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_inh <- RunHarmony(
        object = merFISH_integrated_inh,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        dims.use = 1:31,
        theta = 3,
        lambda = 0.05,
        plot_convergence = TRUE
)
merFISH_integrated_inh <- RunUMAP(merFISH_integrated_inh, dims = 1:31, 
                                  reduction = 'harmony')

#
KNN_K = determineK(merFISH_integrated_inh, 2, 4, reduction="harmony", dims=1:31)
ggsave("../figure_new//merFISH_integrated_inh_determineK.pdf", width = 4, height = 4)


# clust
merFISH_integrated_inh <- FindNeighbors(merFISH_integrated_inh, k.param = 20,
                                        reduction = "harmony", dims = 1:31)


merFISH_integrated_inh <- FindClusters(merFISH_integrated_inh, 
                                       resolution = 2,
                                       method = "igraph",
                                       algorithm = 4,
                                       group.singletons = F,
                                       verbose = FALSE)

DimPlot(merFISH_integrated_inh, group.by = "RNA_snn_res.2")

# remove cell cluster that contains less than 10 cells
cluster_select = table(merFISH_integrated_inh$RNA_snn_res.2)
cluster_select = cluster_select[cluster_select > 10]
# remove singleton
cluster_select = cluster_select[names(cluster_select) != "singleton"]

merFISH_integrated_inh = subset(merFISH_integrated_inh, subset = RNA_snn_res.2 %in% names(cluster_select)) #49847

         
DimPlot(merFISH_integrated_inh, label = T)

#source("stable_cluster.R")
#stable_c = stable_cluster(merFISH_integrated_inh, 2, "RNA_snn_res.2", 1, reduction="harmony", dims = 1:25)
#stable_c = one_bootstrap(merFISH_integrated_inh, 2, "RNA_snn_res.2", 1,reduction="harmony", dims=1:29)
#merFISH_integrated_inh = subset(merFISH_integrated_inh, idents = stable_c[stable_c$percent >= 0.4, 1] ) # 41732 samples 


# remove doublets, check the makers
merFISH_integrated_inh.markers = FindAllMarkers(merFISH_integrated_inh,
                                                logfc.threshold = 0.5, 
                                                min.pct = 0.2,
                                                only.pos = T)
merFISH_integrated_inh.markers_top <- 
        merFISH_integrated_inh.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

FindMarkers(merFISH_integrated_inh,
               ident.1 = 8,
               logfc.threshold = 0.3, 
               min.pct = 0.1,
               only.pos = T)

FeaturePlot(merFISH_integrated_inh, features = c("Mog", "Gjb6")) 
FeaturePlot(merFISH_integrated_inh, features = c("Gad1", "Gad2")) 

DimPlot(subset(merFISH_integrated_inh, RNA_snn_res.2 %in% 20:24),
        label = T) &
        scale_color_manual(values = as.vector(pals::alphabet())) & NoAxes()

merFISH_integrated_inh = subset(merFISH_integrated_inh, idents = c(8, 17, 20, 21, 23, 24), invert = T) # 42907 samples

cluster_sample = as.matrix(table(merFISH_integrated_inh$orig.ident, merFISH_integrated_inh$RNA_snn_res.2))
pheatmap::pheatmap(log(cluster_sample[,colSums(cluster_sample)!=0]), scale = "row")

DimPlot(merFISH_integrated_inh, label = T, group.by = "RNA_snn_res.2",
        label.size = 5) &
        scale_color_manual(values = as.vector(pals::alphabet())) & NoAxes()

# re-label
Inh_subcluster <- c(1:length(levels(merFISH_integrated_inh)))
names(Inh_subcluster) <- levels(merFISH_integrated_inh)
merFISH_integrated_inh <- RenameIdents(merFISH_integrated_inh, Inh_subcluster)
merFISH_integrated_inh@meta.data$Inh_subcluster = Idents(merFISH_integrated_inh)


merFISH_integrated_inh <- RunUMAP(merFISH_integrated_inh, dims = 1:31, 
                                  reduction = 'harmony')

DimPlot(merFISH_integrated_inh, label = T, 
        label.size = 5) &
        scale_color_manual(values = as.vector(pals::alphabet2())) & NoAxes()

FeaturePlot(merFISH_integrated_inh, features = c("Pvalb", "Kcnc2", "Grin3a", "Lamp5", "Pdyn", "Th", "Penk", "Cnr1"),
            min.cutoff = 0.5, ncol = 4, pt.size = 0.1) &
        scale_color_gradientn(colours = c( "grey95", "pink", "red", "red4")) & NoAxes()


saveRDS(merFISH_integrated_inh, file = "./RDS/merFISH_integrated_inh.seurat_0712.rds")
merFISH_integrated_inh = readRDS(file = "./RDS/merFISH_integrated_inh.seurat.rds")



merFISH_integrated$Inh_subcluster = as.character(merFISH_integrated_inh$Inh_subcluster[colnames(merFISH_integrated)])
merFISH_integrated$Inh_subcluster[is.na(merFISH_integrated$Inh_subcluster)] <- "Other"
merFISH_integrated$Inh_subcluster = factor(merFISH_integrated$Inh_subcluster, levels = c(levels(merFISH_integrated_inh$Inh_subcluster), "Other"))
plot_cell_location(merFISH_integrated, "Inh_subcluster", color = c(as.vector(pals::alphabet2(20)), "gray90") , "../data/20211208_PFCL6/feature_metadata.csv", 1:500)
plot_cell_location(merFISH_integrated, "Inh_subcluster", color = c(as.vector(pals::alphabet2(20)), "gray90") , "../data/20220513_PFCL1_Results/feature_metadata.csv", 1:500)


SpatialDimPlot(merFISH_integrated_inh,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_inh), 
               images = "image.3",
               crop = T,
               stroke=0,
               pt.size.factor = 3,
               ncol = 10,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 2)

VlnPlot(merFISH_integrated, features = "nCount_RNA", pt.size = -1, group.by = "orig.ident")

#
merFISH_integrated_inh.markers = FindAllMarkers(merFISH_integrated_inh,
                                                logfc.threshold = 0.5,
                                                min.pct = 0.3,
                                                only.pos = T)
write.csv(merFISH_integrated_inh.markers, "./DEGs/merFISH_integrated_inh.markers_0712.csv", quote = F)


merFISH_integrated_inh.markers_top <- 
        merFISH_integrated_inh.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_inh, features = unique(merFISH_integrated_inh.markers_top$gene),
        scale.by = "radius", col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c(pals::coolwarm(10)[2:10]))


SpatialPlot(merFISH_integrated_inh,
            images = paste0("image.", c(1:3,17:18)),
            stroke=0.1,
            ncol = 3,
            pt.size.factor = 2) & theme(aspect.ratio = 1.5) &
        scale_fill_manual(values = as.vector(pals::alphabet2())) & NoAxes()&
        theme(legend.position = "none")




###  non-neuron
merFISH_integrated_non = subset(merFISH_integrated, idents = levels(merFISH_integrated)[3:4], invert = T) # 168399

merFISH_integrated_non = ScaleData(merFISH_integrated_non, vars.to.regress = c("abs_volume")) # regress out abs_volume
merFISH_integrated_non = RunPCA(merFISH_integrated_non, npcs = 50, verbose = FALSE)
ElbowPlot(merFISH_integrated_non, ndims=50)
determinePCA(merFISH_integrated_non) # 26

DimPlot(merFISH_integrated_non)

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_non <- RunHarmony(
        object = merFISH_integrated_non,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        dims.use = 1:26,
        plot_convergence = TRUE
)

merFISH_integrated_non <- RunUMAP(merFISH_integrated_non, dims = 1:26, 
                                  reduction = 'harmony')


KNN_K = determineK(merFISH_integrated_non, 2, 4, reduction="harmony", dims=1:26)
ggsave("../figure_new//merFISH_integrated_non_determineK.pdf", width = 4, height = 4)


# clust
merFISH_integrated_non <- FindNeighbors(merFISH_integrated_non, k.param = 15,
                                        reduction = "harmony", dims = 1:26)

merFISH_integrated_non <- FindClusters(merFISH_integrated_non, 
                                       resolution = 2,
                                       method = "igraph",
                                       algorithm = 4,
                                       group.singletons = F,
                                       verbose = FALSE)

DimPlot(merFISH_integrated_non)


FeaturePlot(merFISH_integrated_non, 
            features = c("nCount_RNA", "abs_volume"),
            min.cutoff = 0.5, pt.size = 0.1,raster=FALSE) &
        scale_color_gradientn(colours = c( "grey95", "pink", "red", "red4")) & NoAxes()

FeaturePlot(merFISH_integrated_non, features = c("Gjb6", "Cldn5", "C1qb", "Mog", "Pdgfra", "Serpinf1", "Gad1", "Grin2a"),
            min.cutoff = 0.5, ncol = 4, pt.size = 0.1) &
        scale_color_gradientn(colours = c( "grey95", "pink", "red", "red4")) & NoAxes()
FeaturePlot(merFISH_integrated_non, features = c("Enpp6"),
            min.cutoff = 0.5, pt.size = 0.1) &
        scale_color_gradientn(colours = c( "grey95", "pink", "red", "red4")) & NoAxes()


# remove cell cluster that contains less than 10 cells
cluster_select = table(merFISH_integrated_non$RNA_snn_res.2)
cluster_select = cluster_select[cluster_select > 10]
# remove singleton
cluster_select = cluster_select[names(cluster_select) != "singleton"]

merFISH_integrated_non = subset(merFISH_integrated_non, subset = RNA_snn_res.2 %in% names(cluster_select)) #168397

Idents(merFISH_integrated_non) = "RNA_snn_res.2"            

DimPlot(merFISH_integrated_non, label = T, group.by = "L1_cluster")

#source("stable_cluster.R")
#stable_c = stable_cluster(merFISH_integrated_non, 2, "RNA_snn_res.2", 1, reduction="harmony", dims = 1:23)
#stable_c = one_bootstrap(merFISH_integrated_non, 2, "RNA_snn_res.2", 1,reduction="harmony", dims=1:23)
#merFISH_integrated_non = subset(merFISH_integrated_non, idents = stable_c[stable_c$percent >= 0.4, 1] ) #
#DimPlot(merFISH_integrated_non, label = T)

# remove doublets
merFISH_integrated_non.markers = FindAllMarkers(merFISH_integrated_non,
                                                logfc.threshold = 0.5, 
                                                min.pct = 0.2,
                                                only.pos = T)
merFISH_integrated_non.markers_top <- 
        merFISH_integrated_non.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

merFISH_integrated_non = subset(merFISH_integrated_non, subset = RNA_snn_res.2 %in% c(24,25,20,13,6,26,19,21,22,17), invert = T) # 139823 samples


# re-label
Non_subcluster <- c("Oligo 1", "Endo 5", "Endo 2", "Astro 3", "Endo 3", "OPC 2", "Microglia 2", "Endo 3", "Astro 1", "Endo 4",
                    "Astro 2", "Endo 1", "Microglia 1", "OPC 1", "VLMC", "Oligo 2")
names(Non_subcluster) <- levels(merFISH_integrated_non)

SpatialDimPlot(merFISH_integrated_non,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_non), 
               images = "image.18",
               crop = T,
               stroke=0,
               pt.size.factor = 2,
               ncol = 8,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 2)

merFISH_integrated_non <- RenameIdents(merFISH_integrated_non, Non_subcluster)
merFISH_integrated_non@meta.data$Non_subcluster = Idents(merFISH_integrated_non)

table(merFISH_integrated_non$orig.ident, merFISH_integrated_non$Non_subcluster)

Idents(merFISH_integrated_non) = factor(Idents(merFISH_integrated_non), levels = sort(levels(merFISH_integrated_non)))
merFISH_integrated_non@meta.data$Non_subcluster = Idents(merFISH_integrated_non)


merFISH_integrated_non <- RunUMAP(merFISH_integrated_non, dims = 1:26, 
                                  reduction = 'harmony')

DimPlot(merFISH_integrated_non, raster=FALSE) &
        scale_color_manual(values = as.vector(pals::kelly()[-c(1:2)])) & NoAxes()

saveRDS(merFISH_integrated_non, file = "./RDS/merFISH_integrated_non.seurat_0715.rds")


merFISH_integrated_non.markers = FindAllMarkers(merFISH_integrated_non,
                                                logfc.threshold = 0.5,
                                                min.pct = 0.2,
                                                only.pos = T)

merFISH_integrated_non.markers_top <- 
        merFISH_integrated_non.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)
DotPlot(object = merFISH_integrated_non, features = unique(merFISH_integrated_non.markers_top$gene),
        scale.by = "radius", col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c(pals::coolwarm(10)[2:10]))



SpatialDimPlot(merFISH_integrated_non,
               images = "image.18", pt.size.factor = 1) &
        theme(aspect.ratio = 1.5) &
        scale_fill_manual(values = as.vector(pals::kelly()[-c(1:2)])) & NoAxes()

SpatialDimPlot(merFISH_integrated_non,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_non), 
               images = "image.3",
               crop = T,
               stroke=0,
               pt.size.factor = 3,
               ncol = 8,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 1.5)


merFISH_integrated$Non_subcluster = as.character(merFISH_integrated_non$Non_subcluster[colnames(merFISH_integrated)])
merFISH_integrated$Non_subcluster[is.na(merFISH_integrated$Non_subcluster)] <- "Other"
merFISH_integrated$Non_subcluster = factor(merFISH_integrated$Non_subcluster, levels = c(levels(merFISH_integrated_non$Non_subcluster), "Other"))
plot_cell_location(merFISH_integrated, "Non_subcluster", color = c(as.vector(pals::kelly()[3:17]), "gray90") , "../data/20211208_PFCL6/feature_metadata.csv", 1:500)
plot_cell_location(merFISH_integrated, "Non_subcluster", color = c(as.vector(pals::kelly()[3:17]), "gray90") , "../data/20220513_PFCL1_Results/feature_metadata.csv", 1:500)






# merge
merFISH_integrated_L2 = merFISH_integrated[, c(colnames(merFISH_integrated_ext),
                                               colnames(merFISH_integrated_inh),
                                               colnames(merFISH_integrated_non))] # 376470

merFISH_integrated_L2$L2_cluster = as.character(merFISH_integrated_L2$L1_cluster) 

merFISH_integrated_L2$L2_cluster[colnames(merFISH_integrated_ext)] <- paste0("Ext ", merFISH_integrated_ext$Ext_subcluster)
merFISH_integrated_L2$L2_cluster[colnames(merFISH_integrated_inh)] <- paste0("Inh ", merFISH_integrated_inh$Inh_subcluster)
merFISH_integrated_L2$L2_cluster[colnames(merFISH_integrated_non)] <- as.character(merFISH_integrated_non$Non_subcluster)


DimPlot(merFISH_integrated_L2, group.by = "L2_cluster", raster=FALSE)


merFISH_integrated_L2$L2_cluster = factor(merFISH_integrated_L2$L2_cluster, 
                                          levels = c(paste0("Ext ", levels(merFISH_integrated_ext$Ext_subcluster)),
                                                     paste0("Inh ", levels(merFISH_integrated_inh$Inh_subcluster)),
                                                     levels(merFISH_integrated_non$Non_subcluster)))
Idents(merFISH_integrated_L2) = "L2_cluster"

merFISH_integrated_L2$CellType = stringr::str_replace(merFISH_integrated_L2$L2_cluster, " .+", "")
merFISH_integrated_L2$CellType = factor(merFISH_integrated_L2$CellType, 
                                          levels = unique(merFISH_integrated_L2$CellType)[c(1,2,6,4,7,5,8,3)])

DimPlot(merFISH_integrated_L2, group.by = "CellType",label = T)

saveRDS(merFISH_integrated_L2, file = "./RDS/merFISH_integrated_L2.seurat_0715.rds")
merFISH_integrated_L2 = readRDS( file = "./RDS/merFISH_integrated_L2.seurat_0715.rds")

# Remove Ext 16, which is out of PFC
merFISH_integrated_L2 = subset(merFISH_integrated_L2, Ext_subcluster != 16)

merFISH_integrated_L2 = RunUMAP(merFISH_integrated_L2, dims = 1:39,
                                reduction = "harmony", #min.dist = 0.1,
                                umap.method = "umap-learn"
                                )

DimPlot(merFISH_integrated_L2,
        group.by = "L1_cluster",
        label.size = 5,
        pt.size = 0.1,
        shuffle = T,
        label = T,
        raster=FALSE) &
        scale_color_manual(values = L1_color) & NoAxes()


FeaturePlot(merFISH_integrated_L2, features = c("Gjb6", "Acsbg1", "Gad1", "Gad2"))


#
subcluster <- c("L6 CT 1", "L4/5 IT 1", "L5/6 NP", "L6 CT 2", "L2/3 IT 2",
                "L2/3 IT 1", "L5 IT 1", "L6 CT 3", "L6 IT 1", "L5 IT 2", 
                "L4/5 IT 2", "L6 CT 4", "L6 IT 2", "L5 IT 3", "L2/3 IT 3",
                "L2/3 IT 4", "L5 ET 1", "L5 ET 2", "Pvalb 1",
                "Lamp5 1", "Sncg 1", "Sncg 2", "Pvalb 2", "Sst 1",
                "Cdca7", "Pvalb 3", "Sst 2", "Vip 1", "Pvalb 4",
                "Sst 3", "Lamp5 2", "Pvalb 5", "Pvalb 6", "Sst 4",
                "Vip 2", "Sst 5", "Sst 6", "Lamp5 3", levels(merFISH_integrated_L2)[-c(1:38)])
names(subcluster) <- levels(merFISH_integrated_L2)


merFISH_integrated_L2 <- RenameIdents(merFISH_integrated_L2, subcluster)
Idents(merFISH_integrated_L2) = factor(Idents(merFISH_integrated_L2), 
                                       levels = c(sort(subcluster[1:18]), sort(subcluster[19:38]), levels(merFISH_integrated_L2)[-c(1:38)]))
merFISH_integrated_L2@meta.data$subcluster = Idents(merFISH_integrated_L2)


#
merFISH_integrated_L2@meta.data$L3_cluster = stringr::str_replace(Idents(merFISH_integrated_L2), " [0-9]", "")
merFISH_integrated_L2@meta.data$L3_cluster = factor(merFISH_integrated_L2@meta.data$L3_cluster,
                                                    levels = c("L2/3 IT", "L4/5 IT", "L5 IT", "L6 IT",
                                                               "L5 ET", "L5/6 NP", "L6 CT",
                                                               "Lamp5", "Pvalb", "Sncg", "Sst", "Vip", "Cdca7",
                                                               "Astro", "Endo", "Microglia", "Oligo", "OPC", "VLMC"))


L3_color = c( pals::parula(7), as.vector(pals::plasma(6)),
             "#AA0DFE","#993F00","#F99379","#604E97", "maroon2", "yellow3")

DimPlot(merFISH_integrated_L2,
        group.by = "L3_cluster",
        label.size = 4,
        pt.size = 0.1,
        label = T,
        shuffle = T,
        raster=FALSE) &
        scale_color_manual(values = L3_color) & NoAxes()


saveRDS(merFISH_integrated_L2, file = "./RDS/merFISH_integrated_L2.seurat_0716.rds")
merFISH_integrated_L2 = readRDS(file = "./RDS/merFISH_integrated_L2.seurat_0716.rds")

# remove Cdca7 cluster
merFISH_integrated_L3 = subset(merFISH_integrated_L2, L3_cluster != "Cdca7")

saveRDS(merFISH_integrated_L3, file = "./RDS/merFISH_integrated_L3.seurat_0716.rds")
merFISH_integrated_L3 = readRDS(file = "./RDS/merFISH_integrated_L3.seurat_0716.rds")


DimPlot(merFISH_integrated_L3,
        group.by = "L1_cluster",
        label.size = 4,
        pt.size = 0.01,
        shuffle = T,
        label = T,
        raster=FALSE) &
        scale_color_manual(values = L1_color) & NoAxes()


L3_color = c( pals::parula(7), as.vector(pals::plasma(5)),
              "#AA0DFE","#993F00","#F99379","#604E97", "maroon2", "yellow3")



DimPlot(merFISH_integrated_L3,
        group.by = "L3_cluster",
        label.size = 4,
        pt.size = 0.01,
        shuffle = T,
        label = T,
        raster=FALSE) &
        scale_color_manual(values = L3_color) & NoAxes()





L2_color = c( pals::parula(18), as.vector(pals::plasma(20)),
              as.vector(pals::alphabet2(20)))


DimPlot(merFISH_integrated_L3,
        group.by = "L2_cluster",
        label.size = 4,
        pt.size = 0.01,
        shuffle = T,
        label = F,
        raster=FALSE) &
        scale_color_manual(values = L2_color) & NoAxes()




# QC
DimPlot(object = merFISH_integrated_L3,
            group.by = c("orig.ident"),
            pt.size = .01,
        shuffle = T,
            raster=FALSE) &
        NoAxes() &
        scale_color_manual(values = as.character(pals::alphabet())) 

Idents(merFISH_integrated_L3) = "L3_cluster"
VlnPlot(merFISH_integrated_L3, features = "nFeature_RNA", pt.size = 0) &
        scale_fill_manual(values = L3_color) &
        theme(legend.position = "none")
ggsave("../figure/merFISH_QC_nFeature_RNA.pdf", width = 10, height = 3.5)

VlnPlot(merFISH_integrated_L3, features = "nCount_RNA", pt.size = 0) &
        scale_fill_manual(values = L3_color) &
        theme(legend.position = "none")
ggsave("../figure/merFISH_QC_nCount_RNA.pdf", width = 10, height = 3.5)

VlnPlot(merFISH_integrated_L3, features = "abs_volume", pt.size = 0) &
        scale_fill_manual(values = L3_color) &
        theme(legend.position = "none")




# train a glm
library(nnet)
cells_embed = merFISH_integrated_L2@reductions$harmony@cell.embeddings[, 1:30]
cells_embed = as.data.frame(cells_embed)
cells_embed$CT = as.factor(merFISH_integrated_L2$CellType)

glm.fit=multinom(CT~., data=cells_embed)
#summary(glm.fit)

#Prediction
model_pred = predict(glm.fit, cells_embed[,-31], "probs")

doublets = apply(model_pred, 1, function(x) ifelse(max(x) >= 0.8, "singlets", "doublets"))

merFISH_integrated_L2_1 = merFISH_integrated_L2
merFISH_integrated_L2_1$doublets = doublets

DimPlot(merFISH_integrated_L2_1,
        group.by = "CellType",
        split.by = "doublets",
        label.size = 5,
        raster=FALSE) &
        scale_color_manual(values = L1_color) & NoAxes()

# L1 marker
Idents(merFISH_integrated_L2) = "L1"
merFISH_integrated_L2.markers = FindAllMarkers(merFISH_integrated_L2,
                                               logfc.threshold = 1,
                                               min.pct = 0.3,
                                               only.pos = T)

merFISH_integrated_L2.markers_top <- 
        merFISH_integrated_L2.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_L2, features = unique(merFISH_integrated_L2.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_L1_markers.dotplot.pdf", width = 12, height = 5)

# L2 marker
Idents(merFISH_integrated_L2) = "L2"
merFISH_integrated_L2.markers = FindAllMarkers(merFISH_integrated_L2,
                                               logfc.threshold = 1,
                                               only.pos = T)

merFISH_integrated_L2.markers_top <- 
        merFISH_integrated_L2.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_L2, features = unique(merFISH_integrated_L2.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_L2_markers.dotplot.pdf", width = 12, height = 12)


# cell marker
DimPlot(merFISH_integrated_L2,
        group.by = "CellType",
        label.size = 5,
        shuffle = T,
        raster=FALSE) &
        scale_color_manual(values = L1_color[-c(2:7, 9)]) & NoAxes()

Idents(merFISH_integrated_L2) = "CellType"
merFISH_integrated_L2.markers = FindAllMarkers(merFISH_integrated_L2,
                                               logfc.threshold = 1,
                                               min.pct = 0.5,
                                               only.pos = T)

merFISH_integrated_L2.markers_top <- 
        merFISH_integrated_L2.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_L2, features = unique(merFISH_integrated_L2.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_CellType_markers.dotplot.pdf", width = 8, height = 3.5)


#
plot_cell_location(merFISH_integrated_L3, "L3_cluster", color = L3_color , "../data/20220514_PFCL1_Results/", 1:350)
plot_cell_location(merFISH_integrated_L3, "L3_cluster", color = L3_color , "../data/20220514_PFCL1_Results/", 380:1000)

plot_expr_spots(merFISH_integrated_L3, c("Syt6"), color = c("magenta"), "../data/20211208_PFCL6", 1:236)
plot_expr_spots(merFISH_integrated_L3, c("Gjb6", "Otof", "Cux2", "Rorb", "Syt6"), color = c("pink", "green", "blue", "yellow", "magenta"), "../data/20211208_PFCL6", 1:236)
plot_expr_spots(merFISH_integrated_L3, c("Gjb6", "Otof", "Cux2", "Rorb", "Syt6"), color = c("pink", "green", "blue", "yellow", "magenta"), "../data/20220514_PFCL1_Results/", 1:350)

gate_points(merFISH_integrated_L2, image = "image.3")

plot_expr_spots(merFISH_integrated_L2, c("Grin2a", "Gad1"), color = c("yellow", "magenta"), "../data/20211208_PFCL6", 1:236)
plot_expr_spots(subset(merFISH_integrated_L2, centroid_1 < 2000 & centroid_1 > 0 & centroid_2 > 2500 & centroid_2 < 4000), c("Grin2a", "Gad1"), color = c("yellow", "magenta"), "../data/20211208_PFCL6", 1:235)

#
plot_expr_spots(merFISH_integrated_L2, c("Gjb6", "Cux2", "Rorb", "Fezf2", "Syt6"), color = c("yellow",  "blue", "green", "red", "magenta"), "../data/20211208_PFCL6", 1:236, plot_cell = F)
plot_expr_spots(merFISH_integrated_L2, c("Gjb6", "Cux2", "Rorb", "Tshz2", "Syt6"), color = c("yellow",  "blue", "green", "red", "magenta"), "../data/20211208_PFCL6", 240:500, plot_cell = F)


#
plot_cell_location(merFISH_integrated_L2,
                   "L1", color = L1_color,
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:236)
plot_cell_location(merFISH_integrated_L2,
                   "L2", color = L2_color,
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:236)
plot_cell_location(merFISH_integrated_L2,
                   "L2", color = c(rep("gray90", 26), L2_color[27], rep("gray90", 10)),
                   "../data/20211208_PFCL6/feature_metadata.csv", 1:236)


Idents(merFISH_integrated_L2) = "L1"
SpatialDimPlot(merFISH_integrated_L2,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_L2), 
               images = "image.3",
               crop = T,
               stroke=0,
               pt.size.factor = 2,
               ncol = 6,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 1.5)
Idents(merFISH_integrated_L2) = "L2"
SpatialDimPlot(merFISH_integrated_L2,
               cells.highlight = CellsByIdentities(object = merFISH_integrated_L2), 
               images = "image.3",
               crop = T,
               stroke=0,
               pt.size.factor = 2,
               ncol = 10,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
        theme(aspect.ratio = 1.5)






#
# correlation with bulk RNAseq
bulk_expr = read.table("/nfs4/chaozhang/proj/Neuron/Aritra/bulk/stringtie/merge.TPM.txt", row.names = 1, header = T)

merFISH_avg_expr = AverageExpression(merFISH_integrated_L2, group.by = "orig.ident", slot = "counts")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 
ggplot(as.data.frame(merFISH_avg_expr), aes(log(mouse3.4 + 0.01), log(mouse4.1+ 0.01))) +
        geom_point(size=1, color="blue4") + 
        geom_abline(slope = 1, color = "gray60")+ 
        xlab("MERFISH rep1 (counts/cell)") + ylab("MERFISH rep2 (counts/cell)") +
        cowplot::theme_cowplot() + ggtitle("Spearman cor: 0.960")
ggsave("../figure/merFISH_rep_cor.scatter.pdf", width = 4, height = 4)
cor(log(merFISH_avg_expr$mouse3.4 + 0.01), log(merFISH_avg_expr$mouse4.1 + 0.01), method = "spearman")

merFISH_avg_expr.cor = cor(merFISH_avg_expr, method = "spearman")
pheatmap::pheatmap(merFISH_avg_expr.cor)

merFISH_avg_expr = as.data.frame(rowMeans(merFISH_avg_expr[,c(10,11)]))
genes = intersect(rownames(merFISH_avg_expr), rownames(bulk_expr))

merFISH_bulk = cbind(merFISH_avg_expr[genes, ], bulk_expr[genes, ])
colnames(merFISH_bulk) = c("merFISH", "bulk_rep1", "bulk_rep2")
merFISH_bulk$bulk = rowMeans(merFISH_bulk[, c(2:3)])

ggplot(merFISH_bulk, aes(log(bulk + 1), log(merFISH+ 0.01))) +
        geom_point(size=1, color="blue4") + 
        cowplot::theme_cowplot() + ggtitle("Spearman cor: 0.773")
ggsave("../figure/merFISH_bulk_cor.scatter.pdf", width = 4, height = 4)
cor(log(merFISH_bulk$merFISH + 0.01), log(merFISH_bulk$bulk + 1), method = "spearman")


# correlation with scRNAseq

# scRNA
load("../scRNA/all_cell_inDD_sobj.RData")
all_cell_inDD_sobj.expr = AverageExpression(all_cell_inDD_sobj, assays = "RNA", group.by = "Sample", slot = "counts")
all_cell_inDD_sobj.expr = as.data.frame(all_cell_inDD_sobj.expr$RNA)
all_cell_inDD_sobj.expr = as.data.frame(rowMeans(all_cell_inDD_sobj.expr[,c(10,11)]))

genes = intersect(rownames(merFISH_avg_expr), rownames(all_cell_inDD_sobj.expr))

merFISH_scRNA = cbind(merFISH_avg_expr[genes, ], all_cell_inDD_sobj.expr[genes, ])
merFISH_scRNA = as.data.frame(merFISH_scRNA)
colnames(merFISH_scRNA) = c("merFISH", "scRNA")

ggplot(merFISH_scRNA, aes(log(scRNA + 0.01), log(merFISH + 0.01))) +
        geom_point(size=1, color="blue4") + 
        cowplot::theme_cowplot() + ggtitle("Spearman cor: 0.646")
ggsave("../figure/merFISH_scRNA_cor.scatter.pdf", width = 4, height = 4)
cor(log(merFISH_scRNA$merFISH + 0.01), log(merFISH_scRNA$scRNA + 0.01), method = "spearman")


harmony_embed = aggregate(x= harmony_embed,     
                          by = list(merFISH_integrated_L2$L2_cluster),      
                          FUN = mean)

#
Idents(merFISH_integrated_L2) = "L2"
DimPlot(merFISH_integrated_L2) 

harmony_embed = merFISH_integrated_L2@reductions$harmony@cell.embeddings
harmony_embed = as.data.frame(harmony_embed)
harmony_embed = aggregate(x= harmony_embed,     
          by = list(merFISH_integrated_L2$L2),      
          FUN = mean)
rownames(harmony_embed) = harmony_embed$Group.1
harmony_embed = harmony_embed[, -1]
harmony_embed.cor = cor(t(harmony_embed))
pheatmap::pheatmap(harmony_embed.cor)



#
Idents(merFISH_integrated_L3) = "L3_cluster"
merFISH_integrated_L3.markers = FindAllMarkers(merFISH_integrated_L3,
                                               logfc.threshold = 0.5,
                                               min.pct = 0.3,
                                               only.pos = T)

merFISH_integrated_L3.markers_top <- 
        merFISH_integrated_L3.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_L3, features = unique(merFISH_integrated_L3.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))



g = c("Otof", "Cux2", "Fezf2", "Tshz2", "Syt6")
g = c("Pvalb", "Sncg", "Cnr1", "Lamp5", "Reln")
FeaturePlot(merFISH_integrated_L3, features = g,
            min.cutoff = "q02", max.cutoff = "q98", ncol = 5) &
        #scale_color_gradientn(colours = as.character(pals::viridis(20))) & NoAxes()
        scale_color_gradientn(colours = c("grey80", "grey95", "pink", "red", "red4")) & NoAxes()

SpatialPlot(merFISH_integrated_L3, features =  g,
            images = "image.18", 
            pt.size.factor = 1.5,
            stroke=0,
            min.cutoff = "q02", max.cutoff = "q98", ncol = 5) &
        theme(aspect.ratio = 2) & 
        scale_fill_gradientn(colours = c("grey80", "grey95", "pink", "red", "red4")) & NoAxes()
#

merFISH_avg_expr = AverageExpression(merFISH_integrated_L3, group.by = "subcluster", slot = "scale.data")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 
merFISH_avg_expr.cor = cor(merFISH_avg_expr)
pheatmap::pheatmap(merFISH_avg_expr.cor, clustering_method = "ward.D")


dist_matrix <- as.dist(1-merFISH_avg_expr.cor)
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
plot(dendrogram)

library(ggdendro)
dendrogram_data <- dendro_data(dendrogram)
dendrogram_segments <- dendrogram_data$segments # contains all dendrogram segment data


#install.packages("dendextend")
library(dendextend)

L2_color1 = L2_color
names(L2_color1) = levels(merFISH_integrated_L3$subcluster)
L2_color1 = L2_color1[dendrogram_data$labels$label]
dend110 <- dendrogram %>% color_labels(col = L2_color1) %>% color_branches(col = L2_color1)

pdf("../figure_new/cellType_hclust.pdf", 8,3)
dend110 %>%  plot(main = "Original tree") 
dev.off()


#

Idents(merFISH_integrated_L3) = "subcluster"
seurat_tmp = merFISH_integrated_L3
levels(seurat_tmp) <- (dendrogram_data$labels$label)


g = c( "Pdgfra", "Ctss",  "Cldn5", "Serpinf1", "Mog","Gjb6",
       "Pvalb", "Cnr1", "Lamp5","Sncg", "Reln",  "Lhx6",
       "Otof", "Cux2", "Rorb","Tshz2", "Syt6", "Fezf2",
       "Gad1","Slc32a1", "Grin2b", "Grin2a")

g = c( "Pdgfra", "Ctss",  "Cldn5", "Serpinf1", "Mog","Gjb6",
       "Reln", "Pvalb", "Cnr1", "Lamp5","Sncg", "Lhx6",
         "Otof", "Cux2","Nptx2", "Rorb","Tshz2", "Syt6", "Fezf2",
       "Gad1","Slc32a1", "Grin2b", "Grin2a")

DotPlot(object = seurat_tmp, features = g,
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        xlab("") + ylab("") + coord_flip()+
        scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) +
        theme(axis.title.x=element_blank(),
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
ggsave("../figure_new//Makers_L2_cluster.pdf", width = 10, height = 4)
rm(seurat_tmp) 

#
g = c( "Pdgfra", "Ctss",  "Cldn5", "Serpinf1", "Mog","Gjb6",
       "Lhx6", "Reln", "Pvalb", "Cnr1", "Lamp5","Sncg", 
       "Otof", "Cux2","Nptx2", "Rorb","Tshz2", "Fezf2", "Syt6",  "Pou3f1",
       "Gad1","Slc32a1", "Grin2b", "Grin2a")
temp = merFISH_integrated_L3
Idents(temp) = factor(as.character(temp$subcluster), levels = c("L5 ET 1", "L5 ET 2", "L6 CT 1", "L6 CT 2", "L6 CT 3", "L5/6 NP",
                                               "L6 CT 4", "L4/5 IT 2", "L5 IT 3", "L5 IT 1", "L6 IT 2", "L2/3 IT 2",
                                               "L2/3 IT 1", "L2/3 IT 3", "L2/3 IT 4", "L4/5 IT 1", "L5 IT 2", "L6 IT 1",
                                               "Lamp5 3", "Lamp5 1", "Lamp5 2", "Sncg 1", "Vip 1", "Vip 2", "Pvalb 5", "Pvalb 6",
                                               "Pvalb 1","Pvalb 2","Pvalb 3", "Pvalb 4", "Sst 1", "Sst 5","Sst 6","Sst 2","Sst 4","Sncg 2",
                                               "Sst 3", "Astro 3", "Astro 1", "Astro 2", "Oligo 1", "Oligo 2", "VLMC", "Endo 5", "Endo 3", 
                                               "Endo 4", "Endo 1", "Endo 2", "Microglia 1", "Microglia 2", "OPC 1", "OPC 2"))
DotPlot(object = temp, features = g,
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        xlab("") + ylab("") + coord_flip()+
        scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) +
        theme(axis.title.x=element_blank(),
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
ggsave("../figure_new//Makers_subcluster.pdf", width = 13, height = 5)
rm(temp)
gc()

Idents(merFISH_integrated_L3) = "subcluster"
SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
               images = "image.18",
               stroke=0.1,
               pt.size.factor = 1.5,
               facet.highlight = TRUE) & theme(aspect.ratio = 2) & 
        scale_fill_manual(values = L2_color) & NoAxes()

SpatialPlot(subset(merFISH_integrated_L3, CellType == "Inh"),
            images = "image.18",
            stroke=0.1,
            pt.size.factor = 2,
            facet.highlight = TRUE) & theme(aspect.ratio = 2) & 
        #coord_flip() &
        scale_fill_manual(values = as.character(pals::alphabet2())) & NoAxes()

SpatialPlot(subset(merFISH_integrated_L3, subset = CellType %in% c("Ext", "Inh"), invert = T),
            images = "image.18",
            stroke= 0.1,
            pt.size.factor = 1,
            facet.highlight = TRUE) & theme(aspect.ratio = 2) & 
        scale_fill_manual(values = L2_color[-c(1:38)]) & NoAxes()

SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            images = "image.18",
            stroke=0.1,
            pt.size.factor = 1.5,
            facet.highlight = TRUE) & theme(aspect.ratio = 2) & 
        scale_fill_manual(values = as.character(pals::alphabet2())) & NoAxes()
#


#### rerun umap for ext/inh
# ext
merFISH_integrated_ext = subset(merFISH_integrated_L3, L1_cluster == "Excitatory") # 182853 samples 

#merFISH_integrated_ext = NormalizeData(merFISH_integrated_ext)
merFISH_integrated_ext = ScaleData(merFISH_integrated_ext, vars.to.regress = c("abs_volume")) 
merFISH_integrated_ext = RunPCA(merFISH_integrated_ext, npcs = 50, verbose = FALSE)

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_ext <- RunHarmony(
        object = merFISH_integrated_ext,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        plot_convergence = TRUE
)

merFISH_integrated_ext <- RunUMAP(merFISH_integrated_ext, dims = 1:30,
                                  reduction = 'harmony')

DimPlot(merFISH_integrated_ext, group.by = "subcluster",
        label.size = 4,
        label = T,
        repel = T,
        pt.size = 0.01,
        raster=FALSE) &
        scale_color_manual(values = L2_color) & NoAxes()

FeaturePlot(merFISH_integrated_ext, features = "nCount_RNA") &
        scale_color_gradientn(colours = rainbow(10)[3:10]) 

saveRDS(merFISH_integrated_ext, file = "./RDS/merFISH_integrated_ext.seurat_0721.rds")


merFISH_integrated_ext.markers = FindAllMarkers(merFISH_integrated_ext,
                                               logfc.threshold = 1,
                                               min.pct = 0.5,
                                               only.pos = T)

merFISH_integrated_ext.markers_top <- 
        merFISH_integrated_ext.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_ext, features = unique(merFISH_integrated_ext.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_integrated_ext_markers.dotplot.pdf", width = 11, height = 4.5)


VlnPlot(merFISH_integrated_ext, features = unique(merFISH_integrated_ext.markers_top$gene),
        fill.by = "ident", adjust = 0.8,
        stack = T, pt.size = 0, flip = T) &
        scale_fill_manual(values = L2_color) 

# inh
merFISH_integrated_inh = subset(merFISH_integrated_L2, subset = CT == "Inh") # 27479 samples 

merFISH_integrated_inh = NormalizeData(merFISH_integrated_inh)
merFISH_integrated_inh = ScaleData(merFISH_integrated_inh, vars.to.regress = c("abs_volume")) 
merFISH_integrated_inh = RunPCA(merFISH_integrated_inh, npcs = 50, verbose = FALSE)

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_inh <- RunHarmony(
        object = merFISH_integrated_inh,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        plot_convergence = TRUE
)

merFISH_integrated_inh <- RunUMAP(merFISH_integrated_inh, dims = 1:29,
                                  reduction = 'harmony')

DimPlot(merFISH_integrated_inh, group.by = "L2",
        label.size = 5,
        label = T,
        raster=FALSE) &
        scale_color_manual(values = as.character(pals::alphabet2())) & NoAxes()

FeaturePlot(merFISH_integrated_inh, features = "centroid_3") &
        scale_color_gradientn(colours = rainbow(10)[3:10]) 

saveRDS(merFISH_integrated_inh, file = "./RDS/merFISH_integrated_inh.seurat_0421.rds")



merFISH_integrated_inh.markers = FindAllMarkers(merFISH_integrated_inh,
                                                logfc.threshold = 1,
                                                min.pct = 0.5,
                                                only.pos = T)

merFISH_integrated_inh.markers_top <- 
        merFISH_integrated_inh.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_inh, features = unique(merFISH_integrated_inh.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_integrated_inh_markers.dotplot.pdf", width = 10, height = 4)




# non neuronal 
merFISH_integrated_NN = subset(merFISH_integrated_L2, subset = CT == "Non-neuron") # 84675 samples 

merFISH_integrated_NN = NormalizeData(merFISH_integrated_NN)
merFISH_integrated_NN = ScaleData(merFISH_integrated_NN, vars.to.regress = c("abs_volume")) 
merFISH_integrated_NN = RunPCA(merFISH_integrated_NN, npcs = 50, verbose = FALSE)

# normalize batch effects
# harmony
require(harmony)
merFISH_integrated_NN <- RunHarmony(
        object = merFISH_integrated_NN,
        group.by.vars = 'orig.ident',
        assay.use = "RNA",
        plot_convergence = TRUE
)

merFISH_integrated_NN <- RunUMAP(merFISH_integrated_NN, dims = 1:20,
                                  reduction = 'harmony')

DimPlot(merFISH_integrated_NN, group.by = "L2",
        label.size = 5,
        label = T,
        raster=FALSE) &
        scale_color_manual(values = L2_color[-c(1:30)]) & NoAxes()


saveRDS(merFISH_integrated_NN, file = "./RDS/merFISH_integrated_NN.seurat_0421.rds")

merFISH_integrated_NN.markers = FindAllMarkers(merFISH_integrated_NN,
                                                logfc.threshold = 1,
                                                min.pct = 0.5,
                                                only.pos = T)

merFISH_integrated_NN.markers_top <- 
        merFISH_integrated_NN.markers %>%
        group_by(cluster) %>%
        slice_max(n = 3, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_NN, features = unique(merFISH_integrated_NN.markers_top$gene),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_integrated_NN_markers.dotplot.pdf", width = 9, height = 3.5)

#
merFISH_integrated_L2 = readRDS("./RDS/merFISH_integrated_final_0419.seurat.rds")
DimPlot(merFISH_integrated_L2, group.by = "L2_cluster")
Idents(merFISH_integrated_L2) = "L2"



DotPlot(object = merFISH_integrated_L2, features = unique(c(unique(merFISH_integrated_ext.markers_top$gene),
                                                     unique(merFISH_integrated_inh.markers_top$gene),
                                                     unique(merFISH_integrated_NN.markers_top$gene))),
        col.min = -2, col.max = 2) +
        theme(axis.text.x = element_text(angle = 30, hjust=1)) +
        scale_y_discrete(limits=rev) + xlab("") + ylab("") +
        scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4"))
ggsave("../figure/merFISH_integrated_allCell_markers.dotplot.pdf", width = 25, height = 10)



merFISH_integrated_L2_downsample = subset(merFISH_integrated_L2, downsample = 500)
merFISH_integrated_L2_downsample_stb = stabilize_expr(merFISH_integrated_L2_downsample, neighbor = 6, npcs = 30, weight.NN = 0.8)
DoHeatmap(merFISH_integrated_L2_downsample_stb, 
          features = unique(c(unique(merFISH_integrated_ext.markers_top$gene),
                              unique(merFISH_integrated_inh.markers_top$gene),
                              unique(merFISH_integrated_NN.markers_top$gene))),
          draw.lines=F, #disp.min = -2, disp.max = 2
          size = 3,
          group.colors = L2_color)&
        scale_fill_gradientn(colors = c("lightblue3", "lightblue", "white", "red", "red4")) &
        scale_color_manual(values = L2_color)

merFISH_integrated_L2_downsample_stb = ScaleData(merFISH_integrated_L2_downsample_stb)
VlnPlot(merFISH_integrated_L2_downsample_stb, features = unique(c(unique(merFISH_integrated_ext.markers_top$gene),
                                                                  unique(merFISH_integrated_inh.markers_top$gene),
                                                                  unique(merFISH_integrated_NN.markers_top$gene))),
        fill.by = "ident", adjust = 0.8, slot = "scale.data",
        stack = T, pt.size = 0, flip = T) &
        scale_fill_manual(values = L2_color) 


#





