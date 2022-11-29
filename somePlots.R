# some plots
setwd("/nfs4/chaozhang/proj/Neuron/Aritra/merFISH/script")

source("functions.R")

L3_color = c( pals::parula(7), as.vector(pals::plasma(6))[1:5],
              "#AA0DFE","#993F00","#F99379","#604E97", "maroon2", "yellow3")

L3_color = as.character(pals::polychrome())[3:30]
# plot cell location
plot_cell_location(merFISH_integrated, "L3_cluster", color = L3_color , "../data/20211208_PFCL6", 1:500)

pdf("../figure_new/Plot.cell_loc.pdf", 10, 10)
plot_cell_location(merFISH_integrated, "L3_cluster", color = L3_color , "../data/20220514_PFCL1_Results", 1:360)
dev.off()

L2_color = c(as.character(pals::polychrome()), as.character(pals::kelly()[2:22]))
pdf("../figure_new/Plot.cell_loc_subcluster.pdf", 8, 10)
plot_cell_location(merFISH_integrated, "subcluster", color = L2_color , "../data/20220514_PFCL1_Results", 1:360)
dev.off()


SpatialPlot(merFISH_integrated, feature = "primary_fovID",
            images = "slice.44")

plot_expr_spots(merFISH_integrated, c("Gjb6", "Otof", "Cux2", "Tshz2", "Syt6"),
                color = c("pink", "green", "blue", "yellow", "magenta"),
                "../data/20220514_PFCL1_Results", 1:300, plot_cell = F)

plot_expr_spots(merFISH_integrated, c("Otof"),
                color = c("green"),
                "../data/20220514_PFCL1_Results", 1:360, plot_cell = T)
plot_expr_spots(merFISH_integrated, c("Cux2"),
                color = c("blue"),
                "../data/20220514_PFCL1_Results", 1:360, plot_cell = T)
plot_expr_spots(merFISH_integrated, c("Cnr1"),
                color = c("yellow"),
                "../data/20220514_PFCL1_Results", 1:360, plot_cell = T)



plot_expr_spots(merFISH_integrated, rownames(merFISH_integrated)[1:400],
                color = rainbow(400),
                "../data/20220514_PFCL1_Results", 120:150, plot_cell = T)


plot_expr_spots(merFISH_integrated, rownames(merFISH_integrated)[1:400],
                color = rep("gray90", 400),
                "../data/20220514_PFCL1_Results", 120:150, plot_cell = F)



merFISH_integrated_L3$PFC = "Out"
merFISH_integrated_L3@meta.data[colnames(merFISH_integrated)[is.na(merFISH_integrated$ABA_PFC) == F & merFISH_integrated$ABA_PFC == "In"],
                                "PFC"] <- "In"
merFISH_integrated_L3@meta.data = merFISH_integrated_L3@meta.data[!is.na(merFISH_integrated_L3$orig.ident), ]
p1 = DimPlot((merFISH_integrated_L3),
        group.by = "PFC",
        pt.size = 0.01,
        shuffle = T,
        raster=FALSE) & NoAxes() &
  scale_color_manual(values = c("red2", "gray90"))

p2 = DimPlot(merFISH_integrated_L3,
             group.by = "PFC",
             split.by = "PFC",
             pt.size = 0.01,
             shuffle = T,
             raster=FALSE) & NoAxes()
p1
ggsave("../figure_new/UMAP_in_outPFC.pdf", width = 8, height = 8)

merFISH_integrated_L3$ABA_metaRegion = merFISH_integrated@meta.data[colnames(merFISH_integrated_L3), "ABA_metaRegion"]
merFISH_integrated_L3_1 = subset(merFISH_integrated_L3, ABA_metaRegion %in% c("ACAd", "ACAv", "DP", "ILA", "PL", "ORBm", "MOp","MOs"))
merFISH_integrated_L3_1$ABA_metaRegion = factor(merFISH_integrated_L3_1$ABA_metaRegion,
                                                levels = c("ORBm","DP",  "ILA","PL", "ACAd", "ACAv", "MOs", "MOp"))
DimPlot(merFISH_integrated_L3_1,
             group.by = "ABA_metaRegion",
             pt.size = 0.01,
             shuffle = T,
             raster=FALSE) & NoAxes() &
  scale_color_manual(values = pals::tol(8))
ggsave("../figure_new/UMAP_PFC_subregions.pdf", width = 8, height = 8)

DimPlot(merFISH_integrated_L3_1,
        group.by = "ABA_metaRegion",
        split.by = "ABA_metaRegion",
        pt.size = 0.01,
        shuffle = T,
        raster=FALSE) & NoAxes() &
  scale_color_manual(values = pals::tol(8))
ggsave("../figure_new/UMAP_PFC_subregions_split.pdf", width = 30, height = 6)



DimPlot(merFISH_integrated_L3, group.by = "L3_cluster",
        pt.size = 0.01,
        shuffle = T,
        label = T,
        raster=FALSE) &
  scale_color_manual(values = L3_color) & NoAxes()
ggsave("../figure_new/UMAP_allcells_L3cluster.pdf", width = 8, height = 8)

SpatialPlot(merFISH_integrated, feature = c("Otof", "Cux2", "Fezf2"),
            pt.size.factor = 1,
            slot = "counts",
            images = "slice.44", stroke = 0, max.cutoff = "q97") &
  #scale_fill_gradientn(colours = c("gray90", "pink", "red2") ) &
  scale_fill_gradientn(colours =  as.vector(pals::gnuplot(20))[2:20]) &
  theme_void() &
  theme(panel.background = element_rect(fill = "black")) &
  NoAxes()
ggsave("../figure_new/SpatialPlot_markers", width = 8, height = 20)

SpatialPlot(merFISH_integrated, feature = c("Otof", "Cux2", "Fezf2"),
            pt.size.factor = 1,
            slot = "counts",
            images = "slice.44", stroke = 0, max.cutoff = "q97") &
  scale_fill_gradientn(colours = c("gray10","gray10", rev(rainbow(8)[1:4])) ) &
  #scale_fill_gradientn(colours =  as.vector(pals::gnuplot(20))[2:20]) &
  theme_void() &
  theme(panel.background = element_rect(fill = "black")) &
  NoAxes()




#A-P
mySpatialPlot_example(merFISH_integrated, L2_color, margin_x=2000, margin_y=2000, pt.size=1)

L2_color = c(as.character(pals::polychrome()), as.character(pals::kelly()[2:22]))
pdf("../figure_new/Plot.AP.cell_loc_subcluster.pdf", 20, 10)
plot_cell_location(merFISH_integrated, "subcluster", color = L2_color , "../data/20220128_PFCL1_Results", 1:1000,
                   legend = F)
dev.off()


#
merFISH_integrated_plot = merFISH_integrated@meta.data
my_col = pals::brewer.dark2(7)

cells = c("L2/3 IT")
cells = c("L4/5 IT", "L5 IT", "L5 ET")
cells = c("L5/6 NP", "L6 IT", "L6 CT")
colors = my_col[1:length(unique(merFISH_integrated_plot[merFISH_integrated_plot$L3_cluster %in% cells, "subcluster"]))]
names(colors) = unique(merFISH_integrated_plot[merFISH_integrated_plot$L3_cluster %in% cells, "subcluster"])

merFISH_integrated_plot$color = colors[merFISH_integrated_plot$subcluster]
merFISH_integrated_plot$color[is.na(merFISH_integrated_plot$color) & merFISH_integrated_plot$ABA_PFC == "Out"] <- "gray90"
merFISH_integrated_plot$color[is.na(merFISH_integrated_plot$color)] <- "gray80"

merFISH_integrated_plot1 = merFISH_integrated_plot[merFISH_integrated_plot$slice == 27, ]
ggplot(merFISH_integrated_plot1,
       aes(-reloc_x, reloc_y)) +
  geom_point(colour = merFISH_integrated_plot1$color, size = 0.3) +
  #scale_color_manual(values = c(rep("gray90",8), L2_color)) +
  theme_void()+coord_fixed() + scale_x_reverse()
ggsave("../figure_new/SpatialPlot_L45_subclusters.pdf", width = 8, height = 8)




plot(merFISH_integrated_plot$reloc_x, merFISH_integrated_plot$reloc_y,
     col = merFISH_integrated_plot$color)

SpatialPlot(subset(merFISH_integrated, L3_cluster %in% c("L5/6 NP", "L6 IT", "L6 CT")),
            pt.size.factor = 1,
            slot = "counts",
            images = "slice.27", stroke = 0, max.cutoff = "q97") &
  scale_fill_manual(values = as.character(pals::brewer.dark2(7))) &
  NoAxes()
ggsave("../figure_new/SpatialPlot_L56_subclusters111.pdf", width = 8, height = 8)


#
mySpatialDimPlot(merFISH_integrated,
               cells.highlight = CellsByIdentities(object = merFISH_integrated)[38:52],
               cells.highlight.size.factor = c(1, 1.8),
               images = "slice.44",
               crop = T,
               stroke = 0,
               ncol = 10,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) 
#ggsave("../figure_new/SpatialPlot_Exc_Cluster_exmaple.pdf", width = 20, height = 5)

mySpatialDimPlot(merFISH_integrated,
                 cells.highlight = CellsByIdentities(object = merFISH_integrated),
                 cells.highlight.size.factor = c(1, 1.8),
                 images = "slice.44",
                 crop = T,
                 stroke = 0,
                 ncol = 7,
                 cols.highlight = c("#DE2D26", "grey80"),
                 facet.highlight = TRUE) 
ggsave("../figure_new/SpatialPlot_Exc_Cluster_exmaple1.pdf", width = 15, height = 18)

pdf("../figure_new/SpatialDimPlot_L5ET1_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                       highlight = "L5 ET 1",
                       slice = 44,
                       max.cut = 2.5,
                       ABA_coronal = 37,
                       ABA_color = "left",
                       pt.size = 0.3,
                       plot.axis = F)
dev.off()
pdf("../figure_new/SpatialDimPlot_L45IT2_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                   highlight = "L4/5 IT 2",
                   slice = 44,
                   max.cut = 2.5,
                   ABA_coronal = 37,
                   ABA_color = "left",
                   pt.size = 0.3,
                   plot.axis = F)
dev.off()
pdf("../figure_new/SpatialDimPlot_L6CT2_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                   highlight = "L6 CT 2",
                   slice = 44,
                   max.cut = 2.5,
                   ABA_coronal = 37,
                   ABA_color = "left",
                   pt.size = 0.3,
                   plot.axis = F)
dev.off()
pdf("../figure_new/SpatialDimPlot_L5IT3_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                   highlight = "L5 IT 3",
                   slice = 44,
                   max.cut = 2.5,
                   ABA_coronal = 37,
                   ABA_color = "left",
                   pt.size = 0.3,
                   plot.axis = F)
dev.off()
pdf("../figure_new/SpatialDimPlot_L6IT1_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                   highlight = "L6 IT 1",
                   slice = 44,
                   max.cut = 2.5,
                   ABA_coronal = 37,
                   ABA_color = "left",
                   pt.size = 0.3,
                   plot.axis = F)
dev.off()
pdf("../figure_new/SpatialDimPlot_L45IT1_CCF.pdf", 15, 15)
SpatialDimPlot_CCF(subset(merFISH_integrated, ABA_metaRegion %in%  c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "MOp", "TTd") ),
                   highlight = "L4/5 IT 1",
                   slice = 44,
                   max.cut = 2.5,
                   ABA_coronal = 37,
                   ABA_color = "left",
                   pt.size = 0.3,
                   plot.axis = F)
dev.off()

# interaction
my_col = pals::brewer.set1(7)
cells = c("Pvalb 1", "L5 IT 3")
colors = my_col[1:length(unique(merFISH_integrated_plot[merFISH_integrated_plot$subcluster %in% cells, "subcluster"]))]
names(colors) = unique(merFISH_integrated_plot[merFISH_integrated_plot$subcluster %in% cells, "subcluster"])

merFISH_integrated_plot$color = colors[merFISH_integrated_plot$subcluster]
merFISH_integrated_plot$color[is.na(merFISH_integrated_plot$color) & merFISH_integrated_plot$ABA_PFC == "Out"] <- "gray90"
merFISH_integrated_plot$color[is.na(merFISH_integrated_plot$color)] <- "gray90"

merFISH_integrated_plot1 = merFISH_integrated_plot[merFISH_integrated_plot$slice == 44, ]
ggplot(merFISH_integrated_plot1,
       aes(-reloc_x, reloc_y)) +
  geom_point(colour = merFISH_integrated_plot1$color, size = ifelse(merFISH_integrated_plot1$subcluster %in% cells, 1, 0.3)) +
  #scale_color_manual(values = c(L2_color[1:4], rep("gray90",4))) +
  theme_void()+coord_fixed() + scale_x_reverse()
ggsave("../figure_new/Interaction_exmaple_L5IT3_Pv1.pdf", width = 8, height = 8)



# Channels expr 
channels_select = c("Chrm1", "Chrm2", "Chrna4", "Chrna7",
                    "Adora1", 
                    "Cacna1a", "Cacna1b", "Cacna1e", "Cacna1g", "Cacna1h", "Cacna2d1", "Cacna2d3", "Cacng2", "Cacng3", "Cacng7", "Cacng8",
                    "Clcn2", "Clcn3", "Clcn6", "Clic5", "Clic4",
                    "Drd1", "Drd5",
                    "Gabra5", "Gabrd", 
                    "Gria1", "Gria4", "Grm3", "Grm7",
                    "Hrh1", "Hrh3",
                    "Adra1d", "Adra2a", "Adrb1",
                    "Adcyap1r1", "Pdyn", "Crh", "Tacr1", "Cckbr", "Penk", "Oprk1", "Rxfp1", "Sstr2", "Sstr3", "Sstr4", "Vipr1",
                    "Htr1a", "Htr1b", "Htr2c", "Htr3a",
                    "Scn1a", "Scn1b", "Scn4b", "Scn8a",
                    "Kcna1", "Kcna2", "Kcna6", "Kcnab1", "Kcnab2", "Kcnb1", "Kcnc1", "Kcnc2", "Kcnc3", "Kcnd3", "Kcng1", "Kcnh1", "Kcnh5", "Kcnj4", "Kcnj9", "Kcnq5")

DotPlot(subset(merFISH_integrated, CellType %in% c("Ext", "Inh")), features = channels_select,
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") + 
  scale_color_gradientn(colours = pals::jet(20)) +
  coord_flip() 
ggsave("../figure_new/Channels_expr_scaled.pdf", width = 12, height = 16)


Sodium_select = read.table("genelist/Sodium_channels_genes.txt")
Sodium_select = Sodium_select$V1[c(3:12, 1:2)]

DotPlot(subset(merFISH_integrated, CellType %in% c("Ext", "Inh")), features = Sodium_select,
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") + 
  scale_color_gradientn(colours = pals::jet(20)) +
  coord_flip() 
ggsave("../figure_new/SodiumChannels_expr_scaled.pdf", width = 12, height = 5)



# pain slices


pain_slice = c(32, 35, 40, 42, 45, 57, 59)
pain_ctr_slice = c(33, 34, 41, 43, 44, 58, 60)
merFISH_integrated_plot = merFISH_integrated@meta.data

my_col = pals::alphabet(26)
cells = c("Pain")
colors = my_col[1:18]
names(colors) = levels(merFISH_integrated_plot$L3_cluster)

merFISH_integrated_plot$color = colors[merFISH_integrated_plot$L3_cluster]
merFISH_integrated_plot$color[merFISH_integrated_plot$ABA_PFC == "Out"] <- "gray90"
#merFISH_integrated_plot$color[is.na(merFISH_integrated_plot$color)] <- "gray90"

merFISH_integrated_plot1 = merFISH_integrated_plot[merFISH_integrated_plot$slice %in% pain_slice, ]
ggplot(merFISH_integrated_plot1,
       aes(-reloc_x, reloc_y)) +
  geom_point(colour = merFISH_integrated_plot1$color, size = 0.1) +
  #scale_color_manual(values = c(L2_color[1:4], rep("gray90",4))) +
  theme_void()+coord_fixed() + scale_x_reverse() + facet_wrap(~ slice, ncol = 7)
ggsave("../figure_new/Spatial_Pain_slices.pdf", width = 35, height = 5)


merFISH_integrated_plot1 = merFISH_integrated_plot[merFISH_integrated_plot$slice %in% pain_ctr_slice, ]
ggplot(merFISH_integrated_plot1,
       aes(-reloc_x, reloc_y)) +
  geom_point(colour = merFISH_integrated_plot1$color, size = 0.1) +
  #scale_color_manual(values = c(L2_color[1:4], rep("gray90",4))) +
  theme_void()+coord_fixed() + scale_x_reverse() + facet_wrap(~ slice, ncol = 7)
ggsave("../figure_new/Spatial_PainCtr_slices.pdf", width = 35, height = 5)



# Inh subtype

DotPlot(subset(merFISH_integrated, CellType %in% c("Inh")), features = c("Cck"),
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") + 
  scale_color_gradientn(colours = pals::jet(20)) +
  coord_flip() 

#
mySpatialDimPlot((merFISH_integrated),
                 cells.highlight = CellsByIdentities(object = merFISH_integrated)[38:52],
                 cells.highlight.size.factor = c(1, 1.5),
                 images = "slice.44",
                 crop = T,
                 stroke = 0,
                 ncol = 8,
                 cols.highlight = c("#DE2D26", "grey80"),
                 facet.highlight = TRUE) 
ggsave("../figure_new/SpatialPlot_Non_Cluster_exmaple.pdf", width = 20, height = 10)


# Cell cell communication
DotPlot(merFISH_integrated, features = c("Pdyn", "Oprk1", "Tac1", 'Tacr1', "Crh", "Crhr",
                                         "Vip", "Vipr", "Sst", "SStr"),
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") + 
  scale_color_gradientn(colours = pals::jet(20)) +
  coord_flip() 

