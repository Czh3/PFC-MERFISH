# add new samples
setwd("/nfs4/chaozhang/proj/Neuron/Aritra/merFISH/script")
source("merFISH_pipeline.R")
source("functions.R")

merFISH_Behavior = merFISH_pipeline("../data/20221001_PFCL1_Results_BEHAVIOR-NOT PAIN/", "Behavior")

merFISH_pain10 = merFISH_pipeline("../data/20221012_PFCL1_Results/", "Pain10")
merFISH_pain11 = merFISH_pipeline("../data/20221022_PFCL1_Results/", "Pain11")

SpatialPlot(merFISH_pain10, features =  c("Cux2", "Syt6"),
            stroke = 0, pt.size.factor = 1)

merFISH_integrated = readRDS("./RDS/merFISH_final.0827.RDS")

# integrate
merFISH_integrated = merge(merFISH_integrated, c(merFISH_Behavior, merFISH_pain10, merFISH_pain11))

# scale data
all.genes = rownames(merFISH_integrated)
merFISH_integrated = FindVariableFeatures(merFISH_integrated, selection.method = "vst")
merFISH_integrated = ScaleData(merFISH_integrated, vars.to.regress = c("abs_volume"))
merFISH_integrated = RunPCA(merFISH_integrated, npcs = 50, verbose = FALSE)
ElbowPlot(merFISH_integrated, ndims=50)


# normalize batch effects
# harmony
require(harmony)
merFISH_integrated <- RunHarmony(
  object = merFISH_integrated,
  group.by.vars = 'orig.ident',
  assay.use = "RNA",
  dims.use = 1:40,
  plot_convergence = TRUE
)

merFISH_integrated <- RunUMAP(merFISH_integrated, dims = 1:40,
                              reduction = 'harmony')


merFISH_integrated = subset(merFISH_integrated, abs_volume <= max(merFISH_integrated_L3$abs_volume))

saveRDS(merFISH_integrated, "./RDS/merFISH_addNew.tmp1031.RDS")

library(pals)
DimPlot(object = merFISH_integrated,
        pt.size = .1,
        raster=FALSE, shuffle = T,
        group.by = "orig.ident") &
  scale_color_manual(values =  as.vector(alphabet()))

FeaturePlot(object = merFISH_integrated[, is.na(merFISH_integrated$L3_cluster) == T],
            pt.size = .1,
            raster=FALSE, 
            features = "abs_volume") &
  scale_color_gradientn(colours = rainbow(10)[4:10])



### Nearest neighbor transfer

new_data = colnames(merFISH_integrated)[is.na(merFISH_integrated$subcluster)]
old_data = colnames(merFISH_integrated)[!is.na(merFISH_integrated$subcluster)]

NN <- Seurat:::FindNN(
  object = merFISH_integrated,
  cells1 = new_data,
  cells2 = old_data,
  internal.neighbors = NULL,
  dims = 1:30,
  reduction = "harmony",
  nn.reduction = "harmony",
  k = 11,
  verbose = TRUE
)

neighbors <- Seurat::GetIntegrationData(object = NN, integration.name = "integrated", slot = 'neighbors')
neighbors = as.data.frame(t(neighbors[['nnab']]@nn.idx))
colnames(neighbors) = new_data

neighbors = as.list(neighbors)

neighbors = parallel::mclapply(neighbors, function(x){
  x = merFISH_integrated$subcluster[old_data[x[1:3]]]
  if(sum(x==min(x)) >= 2){
    max(x)
  }
}, mc.cores = 10)

neighbors[sapply(neighbors, is.null)] <- NA
neighbors = unlist(neighbors)


merFISH_integrated$predict_subcluster = merFISH_integrated$subcluster
merFISH_integrated$predict_subcluster[names(neighbors)] = neighbors


merFISH_integrated$subcluster = merFISH_integrated$predict_subcluster
merFISH_integrated$L3_cluster = stringr::str_replace(merFISH_integrated$subcluster, " [0-9]", "")
merFISH_integrated$L3_cluster = factor(merFISH_integrated$L3_cluster,
                                       levels = c("L2/3 IT", "L4/5 IT", "L5 IT", "L6 IT",
                                                  "L5 ET", "L5/6 NP", "L6 CT",
                                                  "Lamp5", "Pvalb", "Sncg", "Sst", "Vip",
                                                  "Astro", "Endo", "Microglia", "Oligo", "OPC", "VLMC"))


merFISH_integrated = merFISH_integrated[, !is.na(merFISH_integrated$subcluster)]


SpatialDimPlot(merFISH_integrated, group.by = "L3_cluster",
               images = "image.22")

# Add more infor
slice_n = 53
k = c(3, 2, 2)
names(k) = c("image.20", "image.21", "image.22")
slice = merFISH_integrated$slice[!is.na(merFISH_integrated$slice)]

for(i in c("image.20", "image.21", "image.22")){
  mat = merFISH_integrated@images[[i]]@coordinates
  
  cell_slice = cutree(hclust(dist(mat[,1:2])), k=k[i])
  plot(mat[,1:2], col = cell_slice)
  cell_slice = cell_slice + slice_n
  slice_n = slice_n + k[i]
  
  slice = c(slice, cell_slice)
}

merFISH_integrated$slice = slice[colnames(merFISH_integrated)]

# add spatial infor
slice_image = lapply(54:60, function(x){
  coord.df = merFISH_integrated@meta.data[merFISH_integrated$slice==x, c(7,8)]
  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
})

names(slice_image) = paste0("slice.", 54:60)

merFISH_integrated@images = c(merFISH_integrated@images[c(1:73,75:77)], slice_image)

saveRDS(merFISH_integrated, "./RDS/merFISH_final.1031.RDS")



L2_color = c( pals::parula(18), as.vector(pals::plasma(19)),
              as.vector(pals::alphabet2(20)))

DimPlot(object = merFISH_integrated[, is.na(merFISH_integrated$L2_cluster) == T],
        group.by = "L3_cluster",
        pt.size = .1,
        raster=FALSE, shuffle = T) &
  scale_color_manual(values = L2_color)

SpatialPlot(subset(merFISH_integrated, L3_cluster == "L5 ET"),
            group.by = "L3_cluster",
            images = "image.20")&
  scale_fill_manual(values = L3_color)


SpatialPlot(merFISH_integrated, features =  c("Cux2"),
            max.cutoff = "q99",
            images = c("image.20", "image.21", "image.22"),
            stroke = 0, pt.size.factor = 1)



