# subcluster

merFISH_integrated = readRDS(file = "./RDS/merFISH_final.0827.RDS")

#ET
merFISH_integrated_ET = subset(merFISH_integrated, L3_cluster == "L5 ET")

merFISH_integrated_ET.marker = FindAllMarkers(merFISH_integrated_ET)

merFISH_integrated_ET.marker = merFISH_integrated_ET.marker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)


DotPlot(merFISH_integrated_ET, features = c( unique(merFISH_integrated_ET.marker$gene)),
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_color_gradientn(colours = c("blue", "white", "red2"))

p1 = SpatialPlot(merFISH_integrated_ET, group.by = "subcluster",
                   images = "slice.44")

p2 = SpatialFeaturePlot(merFISH_integrated_ET, features = c("Nnat", "Nos1",  "Scn4b", "Scn1a"),
                        stroke = 0,
                        pt.size.factor = 2.5,
                        max.cutoff = "q97",
                   images = "slice.44") &
  scale_fill_gradientn(colours = pals::coolwarm(10) )

p1|p2


#CT
merFISH_integrated_CT = subset(merFISH_integrated, L3_cluster == "L6 CT")

merFISH_integrated_CT.marker = FindAllMarkers(merFISH_integrated_CT)

merFISH_integrated_CT.marker = merFISH_integrated_CT.marker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)


DotPlot(merFISH_integrated_CT, features = c( unique(merFISH_integrated_CT.marker$gene)),
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_color_gradientn(colours = c( "blue", "white", "red"))

p1 =SpatialPlot(merFISH_integrated_CT, group.by = "subcluster",
                 stroke = 0.1,
                 pt.size.factor = 1.3,
                 images = c( "slice.45"))

p2 = SpatialFeaturePlot(merFISH_integrated_CT, features = c("Scn4b", "Nnat",  "Nr4a1", "Kcnab1"),
                        stroke = 0,
                        pt.size.factor = 2.5,
                        max.cutoff = "q97",
                        images = "slice.45") &
  scale_fill_gradientn(colours = pals::coolwarm(10) )

p1|p2


#
#Pv
merFISH_integrated_Pv = subset(merFISH_integrated, L3_cluster == "Pvalb")

merFISH_integrated_Pv.marker = FindAllMarkers(merFISH_integrated_Pv)

merFISH_integrated_Pv.marker = merFISH_integrated_Pv.marker %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


p1 = DotPlot(merFISH_integrated_Pv, features = c( unique(merFISH_integrated_Pv.marker$gene)),
        scale = F) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_color_gradientn(colours = c( "blue", "white", "red"))

p1 = VlnPlot(merFISH_integrated_Pv, features = c( unique(merFISH_integrated_Pv.marker$gene)),
        stack = T, pt.size = -1, flip = T, fill.by = "ident") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        legend.position = "none") +
  ggsci::scale_fill_d3()

# Sst
merFISH_integrated_Sst = subset(merFISH_integrated, L3_cluster == "Sst")

merFISH_integrated_Sst.marker = FindAllMarkers(merFISH_integrated_Sst)

merFISH_integrated_Sst.marker = merFISH_integrated_Sst.marker %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DotPlot(merFISH_integrated_Sst, features = c( unique(merFISH_integrated_Sst.marker$gene)),
             scale = F) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_color_gradientn(colours = c( "blue", "white", "red"))

p2 = VlnPlot(merFISH_integrated_Sst, features = c( unique(merFISH_integrated_Sst.marker$gene)),
        stack = T, pt.size = -1, flip = T, fill.by = "ident") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        legend.position = "none") +
  ggsci::scale_fill_d3()

p1|p2
ggsave("../figure_new/Pv_Sst_subtype_markers_violinPlot.pdf", width = 8, height = 12)



#
merFISH_integrated_Inh = subset(merFISH_integrated, CellType == "Inh")
genes = read.table("./genelist/Inh_other_marker.shortList.txt")$V1
genes = genes[genes%in% rownames(merFISH_integrated_Inh)]
VlnPlot(merFISH_integrated_Inh, features = genes,
        stack = T, pt.size = -1, flip = T, fill.by = "ident") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        legend.position = "none")
ggsave("../figure_new/Inh_markers_shortList_violinPlot.pdf", width = 8, height = 6)

genes = read.table("./genelist/Inh_other_marker.longList.txt")$V1
genes = genes[genes%in% rownames(merFISH_integrated_Inh)]
VlnPlot(merFISH_integrated_Inh, features = genes,
        stack = T, pt.size = -1, flip = T, fill.by = "ident") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        legend.position = "none")
ggsave("../figure_new/Inh_markers_longList_violinPlot.pdf", width = 8, height = 12)

