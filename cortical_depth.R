# cortical depth
library(Seurat)

#
merFISH_integrated = readRDS(file = "./RDS/merFISH_final.0827.RDS")
corticalDepth = merFISH_integrated@meta.data
corticalDepth = corticalDepth[!is.na(corticalDepth$ABA_PFC),]
corticalDepth = corticalDepth[corticalDepth$ABA_PFC == "In", ]

ggplot(na.omit(corticalDepth[corticalDepth$CellType == "Ext", ]),
       aes(subcluster, abs(reloc_x)/7000, fill = subcluster, color = subcluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color) +
  scale_color_manual(values = L2_color) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1,0)) + 
  theme(legend.position = "none") + xlab("")+  ylab("Cortial depth") +
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))
ggsave("../figure_new//Cortical_depth_Ext.pdf", width = 12, height = 3)

ggplot(na.omit(corticalDepth[corticalDepth$L1_cluster == "Inhibitory", ]),
       aes(subcluster, abs(reloc_x)/7000, fill = subcluster, color = subcluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color[-c(1:18)]) +
  scale_color_manual(values = L2_color[-c(1:18)]) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1,0)) + 
  theme(legend.position = "none") + xlab("")+ylab("Cortial depth") +
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))

ggsave("../figure_new//Cortical_depth_Inh.pdf", width = 12, height = 3)


ggplot(na.omit(corticalDepth[!corticalDepth$L1_cluster %in% c("Excitatory", "Inhibitory"), ]), 
       aes(subcluster, abs(reloc_x)/7000, fill = subcluster, color = subcluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color[-c(1:37)]) +
  scale_color_manual(values = L2_color[-c(1:37)]) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1,0)) + 
  theme(legend.position = "none") + xlab("")+ylab("Cortial depth") +
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))
ggsave("../figure_new//Cortical_depth_Noneuronal.pdf", width = 10, height = 3)



library(ggridges)
corticalDepth$subcluster = factor(corticalDepth$subcluster, 
                                  levels = rev(levels(Idents(merFISH_integrated))))
ggplot(na.omit(corticalDepth), 
       aes( abs(reloc_x)/7000, subcluster, fill = subcluster))+
  geom_density_ridges(color = "black")+
  scale_fill_manual(values = rev(L2_color[1:52])) +
  theme_bw(base_size = 18) + xlim(0, 1)+
  theme(legend.position = "none") + ylab("")+xlab("Cortial depth") 









# old method
SpatialPlot(merFISH_integrated, images = "slice.25") &
  scale_fill_manual(values = L2_color)

FeaturePlot(merFISH_integrated, features = "Drd1")

plot(merFISH_integrated@images$slice.24@coordinates[, 1:2], pch=20, cex=0.1)
abline(v=9100, col="blue")


slice8.PFC = gate_points(merFISH_integrated, image = "slice.8")
s8 = merFISH_integrated@images$slice.8@coordinates[slice8.PFC, 1:2]
s8$cortical_depth = max(s8$centroid_1) - s8$centroid_1
#s8$cortical_depth = s8$cortical_depth/max(s8$cortical_depth)

slice9.PFC = gate_points(merFISH_integrated, image = "slice.9")
s9 = merFISH_integrated@images$slice.9@coordinates[slice9.PFC, 1:2]
s9$cortical_depth = max(s9$centroid_1) - s9$centroid_1
#s9$cortical_depth = s9$cortical_depth/max(s9$cortical_depth)


slice11.PFC = gate_points(merFISH_integrated, image = "slice.11")
s11 = merFISH_integrated@images$slice.11@coordinates[slice11.PFC, 1:2]
s11$cortical_depth = max(s11$centroid_1) - s11$centroid_1
#s11$cortical_depth = s11$cortical_depth/max(s11$cortical_depth)


slice12.PFC = gate_points(merFISH_integrated, image = "slice.12")
s12 = merFISH_integrated@images$slice.12@coordinates[slice12.PFC, 1:2]
s12$cortical_depth = max(s12$centroid_1) - s12$centroid_1
#s12$cortical_depth = s12$cortical_depth/max(s12$cortical_depth)


slice24.PFC = gate_points(merFISH_integrated, image = "slice.24")
s24 = merFISH_integrated@images$slice.24@coordinates[slice24.PFC, 1:2]
s24$cortical_depth = max(s24$centroid_2) - s24$centroid_2
#s24$cortical_depth = s24$cortical_depth/max(s24$cortical_depth)

slice25.PFC = gate_points(merFISH_integrated, image = "slice.25")
s25 = merFISH_integrated@images$slice.25@coordinates[slice25.PFC, 1:2]
s25$cortical_depth = max(s25$centroid_2) - s25$centroid_2
#s25$cortical_depth = s25$cortical_depth/max(s25$cortical_depth)


slice32.PFC = gate_points(merFISH_integrated, image = "slice.32")
s32 = merFISH_integrated@images$slice.32@coordinates[slice32.PFC, 1:2]
s32$cortical_depth = max(s32$centroid_1) - s32$centroid_1
#s32$cortical_depth = s32$cortical_depth/max(s32$cortical_depth)

slice33.PFC = gate_points(merFISH_integrated, image = "slice.33")
s33 = merFISH_integrated@images$slice.33@coordinates[slice33.PFC, 1:2]
s33$cortical_depth = max(s33$centroid_1) - s33$centroid_1
#s33$cortical_depth = s33$cortical_depth/max(s33$cortical_depth)

slice37.PFC = gate_points(merFISH_integrated, image = "slice.37")
s37 = merFISH_integrated@images$slice.37@coordinates[slice37.PFC, 1:2]
s37$cortical_depth = max(s37$centroid_1) - s37$centroid_1
#s37$cortical_depth = s37$cortical_depth/max(s37$cortical_depth)

slice40.PFC = gate_points(merFISH_integrated, image = "slice.40")
s40 = merFISH_integrated@images$slice.40@coordinates[slice40.PFC, 1:2]
s40$cortical_depth = max(s40$centroid_1) - s40$centroid_1
#s40$cortical_depth = s40$cortical_depth/max(s40$cortical_depth)

slice40.PFC.1 = gate_points(merFISH_integrated, image = "slice.40")
s40.1 = merFISH_integrated@images$slice.40@coordinates[slice40.PFC.1, 1:2]
s40.1$cortical_depth = s40.1$centroid_1 - min(s40.1$centroid_1)

slice42.PFC = gate_points(merFISH_integrated, image = "slice.42")
s42 = merFISH_integrated@images$slice.42@coordinates[slice42.PFC, 1:2]
s42$cortical_depth = max(s42$centroid_1) - s42$centroid_1

slice43.PFC = gate_points(merFISH_integrated, image = "slice.43")
s43 = merFISH_integrated@images$slice.43@coordinates[slice43.PFC, 1:2]
s43$cortical_depth = s43$centroid_2 - min(s43$centroid_2)

slice43.PFC.1 = gate_points(merFISH_integrated, image = "slice.43")
s43.1 = merFISH_integrated@images$slice.43@coordinates[slice43.PFC.1, 1:2]
s43.1$cortical_depth = max(s43.1$centroid_2) - s43.1$centroid_2


slice45.PFC = gate_points(merFISH_integrated, image = "slice.45")
s45 = merFISH_integrated@images$slice.45@coordinates[slice45.PFC, 1:2]
s45$cortical_depth = max(s45$centroid_1) - s45$centroid_1



corticalDepth = rbind(s8, s11, s12, s24, s25, s33, s37, s40, s40.1, s42, s43, s43.1, s45)
corticalDepth$cluster = merFISH_integrated$subcluster[rownames(corticalDepth)]
corticalDepth$L1 = merFISH_integrated$L1_cluster[rownames(corticalDepth)]
corticalDepth$slice = merFISH_integrated$slice[rownames(corticalDepth)]
corticalDepth$CellType = merFISH_integrated$CellType[rownames(corticalDepth)]
corticalDepth$CellType[!corticalDepth$CellType %in% c("Ext", "Inh")] = "Non-neuron"

saveRDS(corticalDepth, "./RDS/corticalDepth.rds")

#corticalDepth$cortical_depth[corticalDepth$cortical_depth > 1000] = 1000
corticalDepth = corticalDepth[corticalDepth$slice %in% c(11, 12, 24, 25, 30, 37,40,42,43),]

ggplot(na.omit(corticalDepth[corticalDepth$CellType == "Ext", ]), aes(cluster, cortical_depth, fill = cluster, color = cluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color) +
  scale_color_manual(values = L2_color) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1000,0)) + 
  theme(legend.position = "none") + xlab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))
  facet_grid(CellType ~ .)
 
ggsave("../figure_new//Cortical_depth_Ext.pdf", width = 12, height = 3)

ggplot(na.omit(corticalDepth[corticalDepth$L1 == "Inhibitory", ]), aes(cluster, cortical_depth, fill = cluster, color = cluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color[-c(1:18)]) +
  scale_color_manual(values = L2_color[-c(1:18)]) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1000,0)) + 
  theme(legend.position = "none") + xlab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))
facet_grid(CellType ~ .)

ggsave("../figure_new//Cortical_depth_Inh.pdf", width = 12, height = 3)

corticalDepth$CellType = NULL
ggplot(na.omit(corticalDepth[!corticalDepth$L1 %in% c("Excitatory", "Inhibitory"), ]), aes(cluster, cortical_depth, fill = cluster, color = cluster))+
  geom_violin(draw_quantiles = T, scale = "width")+
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .5)+
  scale_fill_manual(values = L2_color[-c(1:38)]) +
  scale_color_manual(values = L2_color[-c(1:38)]) +
  theme_bw(base_size = 18) + scale_y_reverse(limits = c(1000,0)) + 
  theme(legend.position = "none") + xlab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1))
ggsave("../figure_new//Cortical_depth_Noneuronal.pdf", width = 10, height = 3)

