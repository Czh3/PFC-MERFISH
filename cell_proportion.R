# cell proportion

# all cell
merFISH_integrated = readRDS(file = "./RDS/merFISH_final.0827.RDS")


merFISH_integrated$CT = as.character(merFISH_integrated$CellType)
merFISH_integrated$CT[!merFISH_integrated$CT %in% c("Ext", "Inh")] <- "Non-neuron"

plot_cell_proportion = function(obj, label){
  cell_prop = as.data.frame(table(obj@meta.data[, label]))
  cell_prop = cell_prop[cell_prop$Freq != 0, ]
  cell_prop$fraction = cell_prop$Freq / sum(cell_prop$Freq)
  
  
  cell_prop$ymax <- cumsum(cell_prop$fraction)
  cell_prop$ymin <- c(0, head(cell_prop$ymax, n=-1))
  
  # Compute label position
  cell_prop$labelPosition <- (cell_prop$ymax + cell_prop$ymin) / 2
  
  # Compute a good label
  cell_prop$label <-  sprintf("%.1f%%", cell_prop$fraction * 100)
  
  cell_prop
}


cell_prop = plot_cell_proportion(merFISH_integrated, "CT")

L1_color = c("#AA0DFE","#993F00",
             pals::kelly()[c(5:8, 10:14)],
             "maroon2", "yellow3")

p1 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = c("#0D77FA", "#CC2096", "yellow3")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()

cell_prop = plot_cell_proportion(subset(merFISH_integrated, CT == "Ext"), "L3_cluster")

p2 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L1_color[1:7]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()


cell_prop = plot_cell_proportion(subset(merFISH_integrated, CT == "Inh"), "L3_cluster")
p3 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L1_color[8:15]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()

p1|p2|p3
ggsave("../figure_new/Cell_proprotion_allRegion.pdf", width = 12, height = 3)





# only PFC region
merFISH_integrated_PFC = subset(merFISH_integrated, ABA_PFC == "In")

cell_prop = plot_cell_proportion(merFISH_integrated_PFC, "CT")
p1 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = c("#0D77FA", "#CC2096", "yellow3")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()


cell_prop = plot_cell_proportion(subset(merFISH_integrated_PFC, CT == "Ext"), "L3_cluster")
p2 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L1_color[1:7]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()


cell_prop = plot_cell_proportion(subset(merFISH_integrated_PFC, CT == "Inh"), "L3_cluster")
p3 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L1_color[8:15]) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()

p1|p2|p3
ggsave("../figure_new//Cell_proprotion_PFC.pdf", width = 12, height = 3)


cell_prop = plot_cell_proportion(merFISH_integrated_PFC, "CT")

# scRNA
load("../scRNA/all_cell_inDD_sobj.RData")
all_cell_inDD_sobj$CT = all_cell_inDD_sobj$CellType
all_cell_inDD_sobj$CT[!all_cell_inDD_sobj$CellType %in% c("Excitatory", "Inhibitory")] <- "Non-neuron"
cell_prop_scRNA = plot_cell_proportion(all_cell_inDD_sobj, "CT")


cell_prop$data = "MERFISH"
cell_prop_scRNA$data = "scRNA"

cell_prop_compare = rbind(cell_prop, cell_prop_scRNA)
cell_prop_compare$Var1[cell_prop_compare$Var1 == "Ext"] <- "Excitatory"
cell_prop_compare$Var1[cell_prop_compare$Var1 == "Inh"] <- "Inhibitory"
cell_prop_compare$Var1 = factor(cell_prop_compare$Var1, levels = c("Excitatory", "Inhibitory", "Non-neuron"))

ggplot(cell_prop_compare, aes( ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = c("#0D77FA", "#CC2096", "yellow3")) +
  xlim(c(2, 4)) +
  scale_y_reverse()+
  theme_void() +
  facet_wrap(~ data)
ggsave("../figure_new//Cell_proprotion_PFC_compare_scRNA.pdf", width = 4, height = 5)


# ABA subregions
merFISH_integrated_plot = merFISH_integrated@meta.data
merFISH_integrated_plot = merFISH_integrated_plot[!is.na(merFISH_integrated_plot$ABA_metaRegion), ]

p1 = ggplot(merFISH_integrated_plot, aes( L3_cluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(c(pals::alphabet(), pals::alphabet2()))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

merFISH_integrated_plot$subcluster = factor(merFISH_integrated_plot$subcluster,
                                                levels = levels(Idents(merFISH_integrated)))
p2 = ggplot(merFISH_integrated_plot, aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(c(pals::alphabet(), pals::alphabet2()))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))


library(cowplot)
ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = .5)

ggsave("../figure_new//Cell_proprotion_All_subregions.pdf", width = 20, height = 8)

L3_color = c( pals::parula(9)[1:7], as.vector(pals::plasma(6)),
              "#AA0DFE","#993F00","#F99379","#604E97", "maroon2", "yellow3")
merFISH_integrated_plot = merFISH_integrated_plot[merFISH_integrated_plot$ABA_metaRegion 
                                                  %in% c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm",
                                                         "MOp", "MOs"),]
  
p1 = ggplot(merFISH_integrated_plot, aes( ABA_metaRegion, fill = L3_cluster)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::polychrome())) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

p2 = ggplot(merFISH_integrated_plot, aes(ABA_metaRegion, fill = subcluster)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values =  c( pals::parula(18), as.vector(pals::plasma(20)),
                                 as.vector(pals::alphabet2(20)))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = .5)

ggsave("../figure_new//Cell_proprotion_All_subregions.pdf", width = 20, height = 8)




# PFC subregions
merFISH_integrated_PFC_plot = merFISH_integrated_PFC@meta.data

p1 = ggplot(merFISH_integrated_PFC_plot, aes( L3_cluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::okabe(8))[2:8]) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

merFISH_integrated_PFC_plot$subcluster = factor(merFISH_integrated_PFC_plot$subcluster,
                                                levels = levels(Idents(merFISH_integrated_PFC)))
p2 = ggplot(merFISH_integrated_PFC_plot, aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::okabe(8))[2:8]) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

p1/p2

library(cowplot)
ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = .5)

ggsave("../figure_new//Cell_proprotion_PFC_subregions.pdf", width = 13, height = 6)



#

SpatialPlot(merFISH_integrated,
            pt.size.factor = 2,
            stroke=0, images = paste0("slice.", c(29,28,30,31,27))) &
  theme(legend.position = "none")&
  scale_fill_manual(values = L2_color) & NoAxes()

mySpatialPlot_example(subset(merFISH_integrated, CellType == "Ext"), L2_color, margin_x=2000, margin_y=2000)
mySpatialPlot_example(subset(merFISH_integrated, CellType == "Inh"), as.character(pals::alphabet()))
mySpatialPlot_example(merFISH_integrated, L2_color, margin_x=2000, margin_y=2000, pt.size = 0.5)

### A - P
cell_prop_ap = as.data.frame(table(merFISH_integrated_PFC$slice, merFISH_integrated_PFC$subcluster))
cell_prop_ap = cell_prop_ap[cell_prop_ap$Var1 %in% c(29,28,30,31,27), ]
cell_prop_ap$Var1 = factor((cell_prop_ap$Var1), levels = c(29,28,30,31,27))
cell_prop_ap$Var1 = as.numeric(cell_prop_ap$Var1)

cell_prop_ap = cell_prop_ap[cell_prop_ap$Var2 %in% levels(Idents(merFISH_integrated))[1:37], ]
cell_prop_ap$Var2 = factor(as.character(cell_prop_ap$Var2), levels = levels(Idents(merFISH_integrated))[1:37])

colnames(cell_prop_ap) = c("Var1", "Ext subtype", "Freq")
cell_prop_ap1_sum <- cell_prop_ap %>% group_by(Var1) %>% dplyr::summarise(Sum=sum(Freq))
cell_prop_ap$percent = cell_prop_ap$Freq / cell_prop_ap1_sum$Sum[cell_prop_ap$Var1]


cell_prop_ap_hp = dcast(cell_prop_ap, `Ext subtype` ~ Var1, value.var = "percent")
rownames(cell_prop_ap_hp) = cell_prop_ap_hp$`Ext subtype`
cell_prop_ap_hp = cell_prop_ap_hp[,-1]
cell_prop_ap_hp[cell_prop_ap_hp>0.18] <- 0.18

pdf("../figure_new/Plot.AP.cell_percent.heatmap.pdf", 14, 3)
pheatmap::pheatmap(t(cell_prop_ap_hp), cluster_rows = F, cluster_cols = F,
                   border_color = NA, display_numbers = F,
                   color = pals::viridis(30))
dev.off()


pdf("../figure_new/Plot.Exc_AP.cell_percent.heatmap.pdf", 14, 3)
pheatmap::pheatmap(t(cell_prop_ap_hp)[,1:18], cluster_rows = F, cluster_cols = F,
                   border_color = NA, display_numbers = F,
                   color = colorRampPalette(c("white", "red" ))(50))
dev.off()
pdf("../figure_new/Plot.Inh_AP.cell_percent.heatmap.pdf", 14, 3)
pheatmap::pheatmap(t(cell_prop_ap_hp)[,19:37], cluster_rows = F, cluster_cols = F,
                   border_color = NA, display_numbers = F,
                   color = colorRampPalette(c("white", "red", "red4" ))(50))
dev.off()

ggplot(cell_prop_ap[cell_prop_ap$`Ext subtype` %in% rownames(cell_prop_ap_hp)[1:18],], aes(Var1, percent,color = `Ext subtype`))+
  geom_point() +
  geom_smooth(method = "loess") +
  scale_color_manual(values = L2_color) +
  theme_cowplot() + xlab("Slice: Anterior to Posterior")


cell_prop_ap = as.data.frame(table(merFISH_integrated$slice, merFISH_integrated$subcluster))
cell_prop_ap = cell_prop_ap[cell_prop_ap$Var1 %in% c(29,28,30,31,27), ]
cell_prop_ap$Var1 = factor(as.numeric(cell_prop_ap$Var1), levels = c(29,28,30,31,27))
cell_prop_ap$Var1 = as.numeric(cell_prop_ap$Var1)

cell_prop_ap1 = cell_prop_ap[cell_prop_ap$Var2 %in% levels(Idents(merFISH_integrated))[1:18], ]
cell_prop_ap1$Var2 = factor(as.character(cell_prop_ap1$Var2), levels = levels(Idents(merFISH_integrated))[1:18])

colnames(cell_prop_ap1) = c("Var1", "Ext subtype", "Freq")
L2_color = c(as.character(pals::polychrome()), as.character(pals::kelly()[2:22]))
ggplot(cell_prop_ap1, aes(Var1, Freq, fill = `Ext subtype`)) +
  geom_bar(stat="identity", position="fill", width = 0.7) +
  scale_fill_manual(values = L2_color[9:30]) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ggsave("../figure_new//Cell_proprotion_AtoP.pdf", width = 4, height = 5)

cell_prop_ap1_sum <- cell_prop_ap1 %>% group_by(Var1) %>% dplyr::summarise(Sum=sum(Freq))
cell_prop_ap1$percent = cell_prop_ap1$Freq / cell_prop_ap1_sum$Sum[cell_prop_ap1$Var1]

ggplot(cell_prop_ap1, aes(Var1, percent, color = `Ext subtype`))+
  geom_point() +
  geom_smooth(method = "loess") +
  scale_color_manual(values = L2_color[9:30]) +
  theme_cowplot() + xlab("Slice: Anterior to Posterior")
ggsave("../figure_new//Cell_proprotion_AtoP.linePlot.pdf", width = 6, height = 5)


cell_prop_ap2 = cell_prop_ap[cell_prop_ap$Var2 %in% levels(Idents(merFISH_integrated))[19:37], ]
cell_prop_ap2$Var2 = factor(as.character(cell_prop_ap2$Var2), levels = levels(Idents(merFISH_integrated))[19:37])

colnames(cell_prop_ap2) = c("Var1", "Inh subtype", "Freq")
ggplot(cell_prop_ap2, aes(Var1, Freq, fill = `Inh subtype`)) +
  geom_bar(stat="identity", position="fill", width = 0.7) +
  scale_fill_manual(values = L2_color[20:38]) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("../figure_new//Cell_proprotion_Inh_AtoP.pdf", width = 4, height = 5)

cell_prop_ap2_sum <- cell_prop_ap1 %>% group_by(Var1) %>% dplyr::summarise(Sum=sum(Freq))
cell_prop_ap2$percent = cell_prop_ap2$Freq / cell_prop_ap2_sum$Sum[cell_prop_ap2$Var1]

ggplot(cell_prop_ap2, aes(Var1, percent,color = `Inh subtype`))+
  geom_point() +
  geom_smooth(method = "loess") +
  scale_color_manual(values = L2_color[20:38]) +
  theme_cowplot() + xlab("Slice: Anterior to Posterior")
ggsave("../figure/Cell_proprotion_Inh_AtoP.linePlot.pdf", width = 6, height = 5)



### A - P # in PFC

cell_prop_ap = as.data.frame(table(merFISH_integrated_PFC$slice, merFISH_integrated_PFC$subcluster))
cell_prop_ap = cell_prop_ap[cell_prop_ap$Var1 %in% c(29,28,30,31,27), ]
cell_prop_ap$Var1 = factor(as.numeric(as.character(cell_prop_ap$Var1)), levels = c(29,28,30,31,27))
cell_prop_ap$Var1 = as.numeric(cell_prop_ap$Var1)

cell_prop_ap1 = cell_prop_ap[cell_prop_ap$Var2 %in% levels(Idents(merFISH_integrated))[1:18], ]
cell_prop_ap1$Var2 = factor(as.character(cell_prop_ap1$Var2), levels = levels(Idents(merFISH_integrated))[1:18])

colnames(cell_prop_ap1) = c("Var1", "Ext subtype", "Freq")
ggplot(cell_prop_ap1, aes(Var1, Freq, fill = `Ext subtype`)) +
  geom_bar(stat="identity", position="fill", width = 0.7) +
  scale_fill_manual(values = L2_color) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

cell_prop_ap1_sum <- cell_prop_ap1 %>% group_by(Var1) %>% dplyr::summarise(Sum=sum(Freq))
cell_prop_ap1$percent = cell_prop_ap1$Freq / cell_prop_ap1_sum$Sum[cell_prop_ap1$Var1]

ggplot(cell_prop_ap1, aes(Var1, percent,color = `Ext subtype`))+
  geom_point() +
  geom_smooth(method = "loess") +
  scale_color_manual(values = L2_color) +
  theme_cowplot() + xlab("Slice: Anterior to Posterior")
ggsave("../figure_new/PFC_Cell_proprotion_AtoP.linePlot.pdf", width = 6, height = 5)


cell_prop_ap2 = cell_prop_ap[cell_prop_ap$Var2 %in% levels(Idents(merFISH_integrated))[19:37], ]
cell_prop_ap2$Var2 = factor(as.character(cell_prop_ap2$Var2), levels = levels(Idents(merFISH_integrated))[19:37])

colnames(cell_prop_ap2) = c("Var1", "Inh subtype", "Freq")
ggplot(cell_prop_ap2, aes(Var1, Freq, fill = `Inh subtype`)) +
  geom_bar(stat="identity", position="fill", width = 0.7) +
  scale_fill_manual(values = L2_color[20:38]) +
  theme_void() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

cell_prop_ap2_sum <- cell_prop_ap1 %>% group_by(Var1) %>% dplyr::summarise(Sum=sum(Freq))
cell_prop_ap2$percent = cell_prop_ap2$Freq / cell_prop_ap2_sum$Sum[cell_prop_ap2$Var1]

ggplot(cell_prop_ap2, aes(Var1, percent,color = `Inh subtype`))+
  geom_point() +
  geom_smooth(method = "loess") +
  scale_color_manual(values = L2_color[19:37]) +
  theme_cowplot() + xlab("Slice: Anterior to Posterior")
ggsave("../figure/PFC_Cell_proprotion_Inh_AtoP.linePlot.pdf", width = 6, height = 5)






### subcluster
## # ext
# cells out PFC
cell_prop = plot_cell_proportion(subset(merFISH_integrated, merFISH_integrated$ABA_PFC == "Out" & CT == "Ext"), "subcluster")

p1 = ggplot(cell_prop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L2_color) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()

# cells in PFC
cell_prop_PFC = plot_cell_proportion(subset(merFISH_integrated, merFISH_integrated$ABA_PFC == "In" & CT == "Ext"), "subcluster")

p2 = ggplot(cell_prop_PFC, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
  scale_fill_manual(values = L2_color) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()
p1|p2

cell_prop_enrichPFC = cbind(cell_prop, cell_prop_PFC)
cell_prop_enrichPFC$ratio = cell_prop_enrichPFC[,10] / cell_prop_enrichPFC[,3]
cell_prop_enrichPFC = cell_prop_enrichPFC[,c(1,15)]
ggplot(cell_prop_enrichPFC, aes(Var1, log2(ratio), fill = Var1))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = L2_color) +
  theme_cowplot() + xlab("") + ylab("#Cells: log2(inPFC/outPFC)")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = "None")
ggsave("../figure_new/Ext_cellType_enriched_PFC.pdf", width = 10, height = 4)

## all cell types
# cells out PFC
cell_prop = plot_cell_proportion(subset(merFISH_integrated, ABA_PFC == "Out"), "subcluster")

# cells in PFC
cell_prop_PFC = plot_cell_proportion(subset(merFISH_integrated, ABA_PFC == "In"), "subcluster")


cell_prop_enrichPFC = cbind(cell_prop, cell_prop_PFC)
cell_prop_enrichPFC$ratio = cell_prop_enrichPFC[,10] / cell_prop_enrichPFC[,3]
cell_prop_enrichPFC = cell_prop_enrichPFC[,c(1,15)]

L2_color = c( pals::parula(18), as.vector(pals::plasma(19)),
              as.vector(pals::alphabet2(20)))
cell_prop_enrichPFC$Var1 = factor(cell_prop_enrichPFC$Var1, levels = levels(Idents(merFISH_integrated)))
ggplot(cell_prop_enrichPFC, aes(Var1, log2(ratio), fill = Var1))+
  geom_bar(stat="identity") +
  scale_fill_manual(values = L2_color) +
  theme_cowplot() + xlab("") + ylab("#Cells: log2(inPFC/outPFC)")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = "None")
ggsave("../figure_new/cellType_enriched_PFC.pdf", width = 15, height = 3.5)





# cluster enriched in PFC

SpatialDimPlot(merFISH_integrated,
            pt.size.factor = 2,
            stroke=0, images = paste0("slice.", c(20:27))) &
  theme(legend.position = "none")& NoAxes()



SpatialDimPlot(subset(merFISH_integrated, CellType == "Ext"),
               cells.highlight = CellsByIdentities(object = subset(merFISH_integrated, CellType == "Ext")), 
               images = "slice.44",
               stroke=0,
               pt.size.factor = 1,
               ncol = 6,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) & theme(aspect.ratio = 0.5)

SpatialDimPlot(subset(merFISH_integrated, CellType == "Ext"),
               cells.highlight = CellsByIdentities(object = subset(merFISH_integrated, CellType == "Ext")), 
               images = "slice.45",
               stroke=0,
               pt.size.factor = 1,
               ncol = 6,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) & theme(aspect.ratio = 0.5) & coord_flip() & scale_x_reverse()


SpatialDimPlot(subset(merFISH_integrated, CellType == "Inh"),
               cells.highlight = CellsByIdentities(object = subset(merFISH_integrated, CellType == "Inh")), 
               images = "slice.44",
               stroke=0,
               pt.size.factor = 1.6,
               ncol = 6,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) & theme(aspect.ratio = 0.5)


# cell types enrich in PFC
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L2/3 IT 1", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L5 ET 1", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L5 IT 1", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L6 CT 2", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L6 CT 3", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L4/5 IT 2", margin_x=2000, margin_y=2000, pt.size = 0.3)

# cell types deplete in PFC
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L6 IT 1", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L4/5 IT 1", margin_x=2000, margin_y=2000, pt.size = 0.3)
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Ext"), "L6 CT 1", margin_x=2000, margin_y=2000, pt.size = 0.3)

#
mySpatialPlot_example2(subset(merFISH_integrated, CellType == "Inh"), "Pvalb3", margin_x=2000, margin_y=2000, pt.size = 0.3)


#
merFISH_integrated_ext = subset(merFISH_integrated, CellType == "Ext")
PFC_ct_plot = merFISH_integrated_ext@meta.data[merFISH_integrated_ext$slice == 44, ]
PFC_ct_plot$subcluster[!PFC_ct_plot$subcluster %in% c("L2/3 IT 1", "L5 ET 1", "L5 IT 1", "L6 CT 2", "L6 CT 3")] <- "Other"
ggplot(PFC_ct_plot, aes(centroid_2, centroid_1)) +
  geom_point(aes(color = subcluster), size = 0.8) +
  theme_void()+
  scale_color_manual(values = c(as.character(pals::kelly(7)[3:7])), "gray95") +
  scale_y_reverse()


merFISH_integrated_ext$PFC_enrich = as.character(merFISH_integrated_ext$subcluster)
merFISH_integrated_ext$PFC_enrich[!merFISH_integrated_ext$PFC_enrich %in% c("L2/3 IT 1", "L5 ET 1", "L5 IT 1", "L6 CT 2", "L6 CT 3")] <- "Other"
merFISH_integrated_ext$PFC_enrich[merFISH_integrated_ext$PFC_enrich == "Other" & merFISH_integrated_ext$ABA_PFC == "In"] <- "Other_PFC"
pdf("../figure_new/PFC_enrich_CellTypes.pdf", 8, 8)
plot_cell_location(merFISH_integrated_ext, "PFC_enrich", color = c(as.vector(pals::kelly()[3:7]), "gray90", "gray80") , "../data/20220514_PFCL1_Results", 1:360, legend = T)
dev.off()

pdf("../figure_new/PFC_regions.pdf", 8, 8)
plot_cell_location(merFISH_integrated_ext, "ABA_PFC", color = c("red2", "gray80") , "../data/20220514_PFCL1_Results", 1:360, legend = T)
dev.off()

merFISH_integrated_ext$PFC_deplete = as.character(merFISH_integrated_ext$subcluster)
merFISH_integrated_ext$PFC_deplete[!merFISH_integrated_ext$PFC_deplete %in% c("L2/3 IT 3","L2/3 IT 4","L4/5 IT 1", "L5 ET 2", "L5 IT 3", "L6 CT 1", "L6 IT 1")] <- "Other"
merFISH_integrated_ext$PFC_deplete[merFISH_integrated_ext$PFC_deplete == "Other" & merFISH_integrated_ext$ABA_PFC == "In"] <- "Other_PFC"
pdf("../figure_new/PFC_deplete_CellTypes.pdf", 8, 8)
plot_cell_location(merFISH_integrated_ext, "PFC_deplete", color = c(as.vector(pals::kelly()[10:16]), "gray90", "gray80") , "../data/20220514_PFCL1_Results", 1:360, legend = T)
dev.off()

merFISH_integrated_ext$unbiased = as.character(merFISH_integrated_ext$subcluster)
merFISH_integrated_ext$unbiased[merFISH_integrated_ext$unbiased %in% c("L2/3 IT 1", "L5 ET 1", "L5 IT 1", "L6 CT 2", "L6 CT 3",
                                                                        "L2/3 IT 3","L2/3 IT 4","L4/5 IT 1", "L5 ET 2", "L5 IT 3", "L6 CT 1", "L6 IT 1")] <- "Other"
merFISH_integrated_ext$unbiased[merFISH_integrated_ext$unbiased == "Other" & merFISH_integrated_ext$ABA_PFC == "In"] <- "Other_PFC"
pdf("../figure_new/PFC_unbiased_CellTypes.pdf", 8, 8)
plot_cell_location(merFISH_integrated_ext, "unbiased", color = c(as.vector(pals::polychrome()[3:8]), "gray90", "gray80") , "../data/20220514_PFCL1_Results", 1:360, legend = T)
dev.off()



# PFC enrich gene
merFISH_integrated_PFC_enrich = subset(merFISH_integrated_L3, CellType == "Ext" & slice %in% c(unique(corticalDepth$slice)))

merFISH_integrated_PFC_enrich$PFC = "Out"
merFISH_integrated_PFC_enrich$PFC[colnames(merFISH_integrated_PFC_enrich) %in% rownames(corticalDepth)] <- "In"

PFC_enrich_gene = FindMarkers(merFISH_integrated_PFC_enrich,
                              pseudocount.use = 0.1,
                              ident.1 = "In",
                              ident.2 = "Out",
                              group.by = "PFC",
                              logfc.threshold = 0)
PFC_enrich_gene$DEG = "No"
PFC_enrich_gene$DEG[PFC_enrich_gene$p_val_adj < 0.01 & PFC_enrich_gene$avg_log2FC > log2(1.2)] <- "Enrich"
PFC_enrich_gene$DEG[PFC_enrich_gene$p_val_adj < 0.01 & PFC_enrich_gene$avg_log2FC < -log2(1.2)] <- "Deplete"
PFC_enrich_gene$DEG = factor(PFC_enrich_gene$DEG, levels = c("Enrich", "Deplete", "No"))

PFC_enrich_gene$gene = rownames(PFC_enrich_gene)
write.csv(PFC_enrich_gene, "./DEGs/PFC_enrich_gene.DEG.csv")


ggplot(PFC_enrich_gene, aes(avg_log2FC, -log10(p_val)))+
  geom_point(aes(color = DEG), shape = 19) +
  theme_cowplot() + xlab("log2FC (inPFC/outPFC)") + ylab("-log10 (P value)") +
  geom_text_repel(
    data = subset(PFC_enrich_gene, abs(avg_log2FC) > 0.5 & -log10(p_val) > 100),
    aes(label = gene),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  scale_color_manual(values = c("red2", "blue2", "gray60")) 
ggsave("../figure_new/PFC_enrich_gene.pdf", width = 7, height = 6)

update_geom_defaults("point", list(shape = 16))

SpatialPlot(merFISH_integrated_ext, features = c("Nnat", "Fezf2", "Scn4b", "Scn1a"),
            stroke = 0, max.cutoff = "q98", #slot = "counts",
            pt.size.factor = 1.6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
ggsave("../figure_new//PFC_enrich_gene.examples.png", width = 12, height = 12, dpi = 600)

SpatialPlot(merFISH_integrated_ext, features = c("Nnat", "Scn4b"),
            stroke = 0, max.cutoff = "q98", #slot = "counts",
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.44")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) 
ggsave("../figure_new//PFC_enrich_gene.examples.pdf", width = 6, height = 12)


pdf("../figure_new//PFC_enrich_gene.spatialPlot.pdf", 5, 5)
for(i in c(PFC_enrich_gene[PFC_enrich_gene$DEG == "Enrich", "gene"],  PFC_enrich_gene[PFC_enrich_gene$DEG == "Deplete", "gene"])){
  p = SpatialPlot(merFISH_integrated_ext, features = i,
              stroke = 0, max.cutoff = "q99",
              pt.size.factor = 1.5,
              images = c("slice.44")) &
    scale_fill_gradientn(colours = pals::gnuplot(10)) &
    ggtitle(paste0("log2FC: ", sprintf("%.2f", PFC_enrich_gene[PFC_enrich_gene$gene==i, "avg_log2FC"]) ))
  print(p)
}
dev.off()

#
DotPlot(merFISH_integrated_ext, 
        features = c(PFC_enrich_gene[PFC_enrich_gene$DEG == "Enrich", "gene"],  PFC_enrich_gene[PFC_enrich_gene$DEG == "Deplete", "gene"]),
        scale = T) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") + 
  scale_color_gradientn(colours = pals::jet(20)) +
  coord_flip() 
ggsave("../figure_new/PFC_enrich_gene_cellTypeExpr.pdf", width = 8, height = 20)


#
PFC_enrich_gene = read.csv("./DEGs/PFC_enrich_gene.DEG.csv", row.names = 1)
PFC_signature_pos = PFC_enrich_gene$gene[PFC_enrich_gene$avg_log2FC > 0][1:10]
PFC_signature_neg = PFC_enrich_gene$gene[PFC_enrich_gene$avg_log2FC < 0][1:10]

PFC_signature = colSums(merFISH_integrated_L3@assays$RNA@data[PFC_signature_pos,]) - 
  colSums(merFISH_integrated_L3@assays$RNA@data[PFC_signature_neg,])
merFISH_integrated_L3$PFC_signature = PFC_signature/max(PFC_signature)

FeaturePlot(merFISH_integrated_L3,
            feature = "PFC_signature",
            #repel = T,
            #max.cutoff = 0.5,
            #min.cutoff = -0.5,
            raster=FALSE) &
  #scale_color_gradientn(colours = c("gray90", "pink", "red2") ) &
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) &
  NoAxes() 
ggsave("../figure_new//PFC_signature.UMAP.pdf", width = 8, height = 8)

SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            feature = "PFC_signature",
            stroke = 0, 
            #max.cutoff = 0.5,
            #min.cutoff = -0.5,
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.44")) &
  scale_fill_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) 
ggsave("../figure_new//PFC_signature.Spatial.pdf", width = 6, height = 6)


# class
PFC_enrich_gene_class = PFC_enrich_gene
pdf("../figure_new/PFC_enrich_genes_classes.pdf", 5, 5)
for(f in list.files("./genelist/", "_genes.txt", full.names = T)){
  genes = read.table(f)$V1
  class = stringr::str_replace(basename(f), "_genes.txt", "")
  PFC_enrich_gene_class$highlight = "No"
  PFC_enrich_gene_class$highlight[PFC_enrich_gene_class$gene %in% genes &
                                    PFC_enrich_gene_class$DEG != "No" ] <- "Yes"

  p = ggplot(PFC_enrich_gene_class, aes(avg_log2FC, -log10(p_val)))+
    geom_point(data = PFC_enrich_gene_class[PFC_enrich_gene_class$highlight == "No",], shape = 19, size = 1, color = "gray80") +
    geom_point(data = PFC_enrich_gene_class[PFC_enrich_gene_class$highlight == "Yes",], shape = 19, size = 2, color = "red2") +
    theme_cowplot() + xlab("log2FC (inPFC/outPFC)") + ylab("-log10 (P value)") +
    geom_text_repel(
      data = subset(PFC_enrich_gene_class, highlight == "Yes"),
      aes(label = gene),
      size = 4,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    ) +
    #scale_color_manual(values = c("gray90", "red2")) +
    ggtitle(class)
  print(p)
}
dev.off()

pdf("../figure_new/PFC_enrich_genes_classes_spatialPlot.pdf", 5, 5)
for(f in list.files("./genelist/", "_genes.txt", full.names = T)){
  genes = read.table(f)$V1
  class = stringr::str_replace(basename(f), "_genes.txt", "")
  PFC_enrich_gene_class$highlight = "No"
  PFC_enrich_gene_class$highlight[PFC_enrich_gene_class$gene %in% genes &
                                    PFC_enrich_gene_class$DEG != "No" ] <- "Yes"
  
  for(g in PFC_enrich_gene_class[PFC_enrich_gene_class$highlight == "Yes", "gene"]){
    p = SpatialPlot(merFISH_integrated_ext, features = g,
                stroke = 0, max.cutoff = "q98", #slot = "counts",
                pt.size.factor = 1.2,
                ncol = 1,
                images = c("slice.44")) &
      scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) 
      p = p + ggtitle(paste0(class, ", ", g, ", DEG:", PFC_enrich_gene_class[g, "DEG"]))
      print(p)
  }
}
dev.off()




#
#
Layer_gene = FindMarkers(merFISH_integrated_L3, ident.1 = "L2/3 IT", ident.2 = "L6 IT", group.by = "L3_cluster")
Layer_gene$gene = rownames(Layer_gene)

Layer_gene = Layer_gene[order(abs(Layer_gene$avg_log2FC), decreasing = T),]
Layer_gene_pos = Layer_gene$gene[Layer_gene$avg_log2FC > 0][1:30]
Layer_gene_neg = Layer_gene$gene[Layer_gene$avg_log2FC < 0][1:30]

Layer_signature = colSums(merFISH_integrated_L3@assays$RNA@scale.data[Layer_gene_pos,]) -
  colSums(merFISH_integrated_L3@assays$RNA@scale.data[Layer_gene_neg,])
merFISH_integrated_L3$Layer_signature = Layer_signature/max(Layer_signature)

FeaturePlot(merFISH_integrated_L3,
            feature = "Layer_signature",
            #repel = T,
            #max.cutoff = 0.5,
            #min.cutoff = -0.5,
            raster=FALSE) &
  #scale_color_gradientn(colours = pals::coolwarm(10)) &
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) &
  NoAxes() 
ggsave("../figure_new//Layer_signature.UMAP.pdf", width = 8, height = 8)

SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
  #subset(merFISH_integrated_L3, L3_cluster %in% c("L2/3 IT", "L4/5 IT", "L5 IT", "L6 IT")),
            feature = "Layer_signature",
            stroke = 0, 
            #max.cutoff = 0.5,
            #min.cutoff = -0.5,
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.44")) &
  #scale_fill_gradientn(colours = rainbow(10)) 
  scale_fill_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) 
ggsave("../figure_new//Layer_signature.Spatial.pdf", width = 6, height = 6)


PFC_enrich_gene_class = PFC_enrich_gene
PFC_enrich_gene_class_count = matrix(nrow = 0, ncol= 2)
cn = NULL
for(f in list.files("./genelist/", "_genes.txt", full.names = T)){
  genes = read.table(f)$V1
  class = stringr::str_replace(basename(f), "_genes.txt", "")
  cn = c(cn, class)
  tmp = PFC_enrich_gene_class[PFC_enrich_gene_class$DEG != "No", ]
  tmp = tmp[tmp$gene %in% genes, ]
  PFC_enrich_gene_class_count = rbind(PFC_enrich_gene_class_count, c(sum(tmp$DEG == "Enrich"), sum(tmp$DEG == "Deplete")))
}

rownames(PFC_enrich_gene_class_count) = cn
colnames(PFC_enrich_gene_class_count) = c("Enrich", "Deplete")
PFC_enrich_gene_class_count1 = PFC_enrich_gene_class_count[rowSums(PFC_enrich_gene_class_count) != 0, ]

pdf("../figure_new/PFC_enrich_genes_classes_quantify_heatmap.pdf", 3.5, 5)
pheatmap::pheatmap(PFC_enrich_gene_class_count1[order(rowSums(PFC_enrich_gene_class_count1), decreasing = T), ], 
                   cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   display_numbers = T,
                   number_format = "%d",
                   number_color = "black",
                   fontsize_number = 14,
                   color = colorRampPalette(c( "white", "orange", "red2"))(50))
dev.off()
#
  



load("../scRNA/all_cell_inDD_sobj.RData")
all_cell_inDD_sobj = subset(all_cell_inDD_sobj, CellType == "Excitatory")
all_cell_inDD_sobj = NormalizeData(all_cell_inDD_sobj)


merFISH_integrated_PFC_enrich = NormalizeData(merFISH_integrated_PFC_enrich)

options(future.globals.maxSize = 8000 * 1024^2)
merFISH_integrated_PFC_enrich_iSpatial = iSpatial::infer(merFISH_integrated_PFC_enrich,
                                                            all_cell_inDD_sobj,
                                                            k.neighbor = 50)
merFISH_integrated_PFC_enrich_iSpatial@images = merFISH_integrated_PFC_enrich[,colnames(merFISH_integrated_PFC_enrich_iSpatial)]@images
saveRDS(merFISH_integrated_PFC_enrich_iSpatial, "./RDS/merFISH_integrated_PFC_enrich_iSpatial.rds")
merFISH_integrated_PFC_enrich_iSpatial = readRDS("./RDS/merFISH_integrated_PFC_enrich_iSpatial.rds")

PFC_enrich_gene_iSpatial = FindMarkers(merFISH_integrated_PFC_enrich_iSpatial,
                              pseudocount.use = 0.1,
                              ident.1 = "In",
                              ident.2 = "Out",
                              group.by = "PFC",
                              logfc.threshold = 0)
PFC_enrich_gene_iSpatial$DEG = "No"
PFC_enrich_gene_iSpatial$DEG[PFC_enrich_gene_iSpatial$p_val_adj < 0.01 & PFC_enrich_gene_iSpatial$avg_log2FC > log2(1.2)] <- "Enrich"
PFC_enrich_gene_iSpatial$DEG[PFC_enrich_gene_iSpatial$p_val_adj < 0.01 & PFC_enrich_gene_iSpatial$avg_log2FC < -log2(1.2)] <- "Deplete"
PFC_enrich_gene_iSpatial$DEG = factor(PFC_enrich_gene_iSpatial$DEG, levels = c("Enrich", "Deplete", "No"))

PFC_enrich_gene_iSpatial$gene = rownames(PFC_enrich_gene_iSpatial)
write.csv(PFC_enrich_gene_iSpatial, "./DEGs/PFC_enrich_gene_iSpatial.DEG.csv")
PFC_enrich_gene_iSpatial = read.csv("./DEGs/PFC_enrich_gene_iSpatial.DEG.csv")

PFC_enrich_gene_iSpatial$measure = "Out"
PFC_enrich_gene_iSpatial$measure[PFC_enrich_gene_iSpatial$gene %in% rownames(merFISH_integrated)] = "In"
library(ggrepel)
ggplot(PFC_enrich_gene_iSpatial, aes(avg_log2FC, -log10(p_val)))+
  geom_point(aes(color = DEG)) +
  theme_cowplot() + xlab("log2FC (inPFC/outPFC)") + ylab("-log10 (P value)") +
  geom_text_repel(
    data = subset(PFC_enrich_gene_iSpatial, abs(avg_log2FC) > 0.5 & -log10(p_val) > 100),
    aes(label = gene, color = measure),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  scale_color_manual(values = c( "blue2", "red2", "black","gray60", "orange")) 
ggsave("../figure_new/PFC_enrich_gene_iSpatial.pdf", width = 7, height = 6)


SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial, features = c("Cdh13", "Cpne7", "Tmem215", "Abcd2"),
            stroke = 0, max.cutoff = "q98",
            pt.size.factor = 1.6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
ggsave("../figure_new//PFC_enrich_gene_iSpatial.examples.png", width = 12, height = 12, dpi = 600)


SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial,
            features = c("Cacna1h", "Cxcl12", "Cdh13", "Abcd2", "Nnat", "Fezf2", "Nr4a1", "Scn4b"),
            stroke = 0, max.cutoff = "q99",
            ncol = 4,
            pt.size.factor = 1.6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("gray10", "gray50", "green", "yellow", "red", "red3")) &
  theme_void() &
  theme(panel.background = element_rect(fill = "black")) &
  # scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
ggsave("../figure_new//PFC_enrich_gene_iSpatial_ColorBlack.examples.png", width = 12, height = 24, dpi = 600)

SpatialPlot(subset(merFISH_integrated, CellType == "Ext"),
            features = c("Cacna1h", "Cxcl12"),
            stroke = 0, max.cutoff = "q98",
            pt.size.factor = 1.6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
ggsave("../figure_new//PFC_enrich_gene.examples1.png", width = 12, height = 6, dpi = 600)

SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial, features = c("Cdh13", "Abcd2"),
            stroke = 0, max.cutoff = "q98",
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
ggsave("../figure_new//PFC_enrich_gene_iSpatial.examples.pdf", width = 12, height = 12, dpi = 600)


SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial, features = c("Kcnc2", "Kcnh7"),
            stroke = 0, max.cutoff = "q98",
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()

#
merFISH_integrated_PFC_enrich_iSpatial[['umap']] = merFISH_integrated_L3[,colnames(merFISH_integrated_PFC_enrich_iSpatial)][['umap']]
FeaturePlot(merFISH_integrated_PFC_enrich_iSpatial, features = c("Figf", "Npr3", "Cxcr7", "Syt6"))
FeaturePlot(merFISH_integrated_PFC_enrich_iSpatial, features = c("Pou3f1", "Foxp2"))

g = c("Otof", "Cux2", "Rorb", "Rspo1", "Osr1", "Syt6")
p1 = SpatialPlot(subset(merFISH_integrated, CellType == "Ext"),
            features = g,
            stroke = 0, max.cutoff = "q98",
            ncol = 6,
            pt.size.factor = 1.6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()

p2 = SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial, features = g,
            stroke = 0, max.cutoff = "q98",
            pt.size.factor = 1.6,
            ncol = 6,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()

p1/p2

#GO
library(enrichR)
dbs = listEnrichrDbs()
dbs = c("GO_Biological_Process_2018", "Allen_Brain_Atlas_up", "KEGG_2016")

enriched <- enrichr(PFC_enrich_gene_iSpatial[PFC_enrich_gene_iSpatial$DEG == "Enrich", "gene"], dbs)
p1 = plotEnrich(enriched[[1]], showTerms = 20, numChar = 80, y = "Count")


enriched <- enrichr(PFC_enrich_gene_iSpatial[PFC_enrich_gene_iSpatial$DEG == "Deplete", "gene"], dbs)
p2 = plotEnrich(enriched[[1]], showTerms = 20, numChar = 80, y = "Count")

p1|p2
ggsave("../figure_new//PFC_enrich_gene_iSpatial.GO.pdf", width = 20, height = 8)


library(org.Mm.eg.db)
library(clusterProfiler)

ego2 <- enrichGO(gene         = PFC_enrich_gene_iSpatial[PFC_enrich_gene_iSpatial$DEG == "Enrich", "gene"],
                 universe = PFC_enrich_gene_iSpatial$gene,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL')

barplot(ego2, showCategory = 10)



ego <- compareCluster(gene~DEG, data = PFC_enrich_gene_iSpatial,
                      fun = "enrichGO",
                 universe = PFC_enrich_gene_iSpatial$gene,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 pvalueCutoff = 0.6,
                 qvalueCutoff = 0.5)

write.csv(as.data.frame(ego), "./DEGs/PFC_enrich_gene_iSpatial.compareGO.csv")

dotplot(ego, showCategory=10) +
  scale_color_gradientn(colours = terrain.colors(7))

ggsave("../figure_new//PFC_enrich_gene_iSpatial.compareGO.pdf", width = 10, height = 6)


SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            feature = stringr::str_split(ego[1,"geneID"], "/")[[1]],
            stroke = 0, 
            pt.size.factor = 1.6,
            ncol = 10,
            max.cutoff = "q98",
            images = c("slice.44")) &
  scale_fill_gradientn(colours = pals::jet(20)) 
ggsave("../figure_new//PFC_enrich_gene_iSpatial.compareGO.pdf", width = 20, height = 6)


SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            feature = stringr::str_split(ego[1,"geneID"], "/")[[1]],
            stroke = 0, 
            pt.size.factor = 1.6,
            ncol = 10,
            max.cutoff = "q98",
            images = c("slice.44")) &
  scale_fill_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) 


SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            feature = stringr::str_split(ego[1,"geneID"], "/")[[1]],
            stroke = 0, 
            pt.size.factor = 1.6,
            ncol = 10,
            max.cutoff = "q98",
            images = c("slice.44")) &
  scale_fill_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) 


#
 
Sodium_channels =  read.table("genelist/Sodium_channels_genes.txt")$V1
Potassium_chennels =  read.table("genelist/Potassium_chennels_genes.txt")$V1
voltage_select = c(Sodium_channels, Potassium_chennels)
voltage_select = voltage_select[voltage_select %in% rownames(merFISH_integrated_L3)]

voltage_signature = colMeans(merFISH_integrated_L3@assays$RNA@scale.data[voltage_select,])
merFISH_integrated_L3$voltage_signature = voltage_signature

FeaturePlot(merFISH_integrated_L3,
            feature = "voltage_signature",
            #repel = T,
            #max.cutoff = 0.5,
            #min.cutoff = -0.5,
            raster=FALSE) &
  #scale_color_gradientn(colours = c("gray90", "pink", "red2") ) &
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) &
  NoAxes() 
ggsave("../figure_new//voltage_signature.UMAP.pdf", width = 8, height = 8)

SpatialPlot(subset(merFISH_integrated_L3, CellType == "Ext"),
            feature = "voltage_signature",
            stroke = 0, 
            max.cutoff = 0.4,
            min.cutoff = -0.2,
            pt.size.factor = 1.6,
            ncol = 1,
            images = c("slice.44")) &
  scale_fill_gradientn(colours = c("blue4", "blue", "gray90", "red", "red3")) 
ggsave("../figure_new//Voltage_signature.Spatial.pdf", width = 8, height = 8)


# GWAS

dbs = listEnrichrDbs()
dbs = c("GWAS_Catalog_2019", "UK_Biobank_GWAS_v1", "ClinVar_2019")

enriched <- enrichr(PFC_enrich_gene_iSpatial[PFC_enrich_gene_iSpatial$DEG == "Enrich", "gene"], dbs)
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count")

enriched <- enrichr(PFC_enrich_gene_iSpatial[PFC_enrich_gene_iSpatial$DEG == "Deplete", "gene"], dbs)
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count")




### PFC Subregion

merFISH_integrated_L3$ABA_metaRegion = merFISH_integrated@meta.data[colnames(merFISH_integrated_L3), "ABA_metaRegion"]
merFISH_integrated_L3_1 = subset(merFISH_integrated_L3, ABA_metaRegion %in% c("ACAd", "DP", "ILA", "PL", "ORBm", "MOp","MOs"))
merFISH_integrated_L3_1$ABA_metaRegion = factor(merFISH_integrated_L3_1$ABA_metaRegion,
                                                levels = c("ORBm","DP",  "ILA","PL", "ACAd", "MOs", "MOp"))

# cells out PFC
merFISH_integrated_L3_Plot = merFISH_integrated_L3_1@meta.data

ggplot(merFISH_integrated_L3_Plot, aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::tol(8))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("../figure_new//Cell_proprotion_PFC_subregions_inculdeMO.pdf", width = 13, height = 4)

ggplot(merFISH_integrated_L3_Plot[merFISH_integrated_L3_Plot$CellType == "Ext", ],
       aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::tol(8))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("../figure_new//Cell_proprotion_PFC_Ext_subregions_inculdeMO.pdf", width = 8, height = 4)


ggplot(merFISH_integrated_L3_Plot[merFISH_integrated_L3_Plot$ABA_metaRegion %in% c("ORBm","DP",  "ILA","PL", "ACAd", "ACAv"),], 
       aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::tol(8))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("../figure_new//Cell_proprotion_PFC_subregions_withoutMO.pdf", width = 13, height = 4)

ggplot(merFISH_integrated_L3_Plot[!merFISH_integrated_L3_Plot$subcluster %in% c("L2/3 IT 3","L2/3 IT 4","L4/5 IT 1", "L5 ET 2", "L5 IT 3", "L6 CT 1", "L6 IT 1") &
                                  merFISH_integrated_L3_Plot$CellType == "Ext" &
                                  merFISH_integrated_L3_Plot$ABA_metaRegion %in% c("ORBm","DP",  "ILA","PL", "ACAd", "ACAv"),], 
       aes(subcluster, fill = ABA_metaRegion)) +
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(values = as.character(pals::tol(8))) +
  cowplot::theme_minimal_vgrid() +
  ylab("Cell proprotion") + xlab("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
ggsave("../figure_new//Cell_proprotion_PFC_subregions_withoutMO_withoutDepletedClusters.pdf", width = 10, height = 4)

#heatmap
merFISH_integrated_L3_Plot_Heatmap = merFISH_integrated_L3_Plot[merFISH_integrated_L3_Plot$CellType == "Ext", ]
merFISH_integrated_L3_Plot_Heatmap = as.data.frame(table(merFISH_integrated_L3_Plot_Heatmap$subcluster, merFISH_integrated_L3_Plot_Heatmap$ABA_metaRegion))
merFISH_integrated_L3_Plot_Heatmap = dcast(merFISH_integrated_L3_Plot_Heatmap, Var1 ~ Var2)
rownames(merFISH_integrated_L3_Plot_Heatmap) = merFISH_integrated_L3_Plot_Heatmap$Var1
merFISH_integrated_L3_Plot_Heatmap = merFISH_integrated_L3_Plot_Heatmap[,-1]
merFISH_integrated_L3_Plot_Heatmap = merFISH_integrated_L3_Plot_Heatmap[levels(Idents(merFISH_integrated_pain_ext)), ]

#random sample
merFISH_integrated_L3_Plot_Heatmap_Random = merFISH_integrated_L3_Plot[merFISH_integrated_L3_Plot$CellType == "Ext", ]
set.seed(1)
merFISH_integrated_L3_Plot_Heatmap_Random$ABA_metaRegion = sample(merFISH_integrated_L3_Plot_Heatmap_Random$ABA_metaRegion)
merFISH_integrated_L3_Plot_Heatmap_Random = as.data.frame(table(merFISH_integrated_L3_Plot_Heatmap_Random$subcluster, merFISH_integrated_L3_Plot_Heatmap_Random$ABA_metaRegion))
merFISH_integrated_L3_Plot_Heatmap_Random = dcast(merFISH_integrated_L3_Plot_Heatmap_Random, Var1 ~ Var2)
rownames(merFISH_integrated_L3_Plot_Heatmap_Random) = merFISH_integrated_L3_Plot_Heatmap_Random$Var1
merFISH_integrated_L3_Plot_Heatmap_Random = merFISH_integrated_L3_Plot_Heatmap_Random[,-1]
merFISH_integrated_L3_Plot_Heatmap_Random = merFISH_integrated_L3_Plot_Heatmap_Random[levels(Idents(merFISH_integrated_pain_ext)), ]


merFISH_integrated_L3_Plot_Heatmap_enrich = log2((merFISH_integrated_L3_Plot_Heatmap+1)/(merFISH_integrated_L3_Plot_Heatmap_Random+1))
merFISH_integrated_L3_Plot_Heatmap_enrich[merFISH_integrated_L3_Plot_Heatmap_enrich < -2.5] <- -2.5
merFISH_integrated_L3_Plot_Heatmap_enrich[merFISH_integrated_L3_Plot_Heatmap_enrich > 2.5] <- 2.5

pdf("../figure_new//Cell_proprotion_PFC_enrich_clusters_subregions.pdf", 6, 8)
pheatmap::pheatmap(merFISH_integrated_L3_Plot_Heatmap_enrich, cluster_cols = F,
                   scale = "none", #clustering_method = "ward.D2",
                   border_color = NA,
                   color = colorRampPalette(c("blue", "white", "red2"))(50))

dev.off()
pdf("../figure_new//Cell_proprotion_PFC_enrich_clusters_subregions_withoutMO.pdf", 5, 8)
pheatmap::pheatmap(merFISH_integrated_L3_Plot_Heatmap_enrich[,-c(6,7)], cluster_cols = F,
                   scale = "none", #clustering_method = "ward.D",
                   border_color = NA,
                   color = colorRampPalette(c("blue", "white", "red2"))(50))
dev.off()

### gene enriched in subregions
merFISH_integrated_L3_1 = subset(merFISH_integrated_L3, ABA_metaRegion %in% c("ACAd", "DP", "ILA", "PL","MOs", "MOp"))

selectSlice = table(merFISH_integrated_L3_1$slice, merFISH_integrated_L3_1$ABA_metaRegion)
selectSlice = apply(selectSlice, 1, function(x) sum(x==0))
selectSlice = names(selectSlice[selectSlice == 0])
merFISH_integrated_L3_1 = subset(merFISH_integrated_L3_1, slice %in% selectSlice)

merFISH_integrated_L3_1$ABA_metaRegion = factor(merFISH_integrated_L3_1$ABA_metaRegion,
                                                levels = c( "DP", "ILA", "PL", "ACAd","MOs", "MOp"))

Idents(merFISH_integrated_L3_1) = "ABA_metaRegion"
PFCsubregion_enrich_gene = FindAllMarkers(merFISH_integrated_L3_1,
                                          test.use = "LR",
                                          latent.vars = "orig.ident",
                                          pseudocount.use = 1,
                                          logfc.threshold = 0.5,
                                          group.by = "ABA_metaRegion")

PFCsubregion_enrich_gene_sig = PFCsubregion_enrich_gene[PFCsubregion_enrich_gene$p_val_adj < 0.05, ]
write.csv(PFCsubregion_enrich_gene, "./DEGs/PFCsubregion_enrich_gene.DEG.csv")

PFCsubregion_enrich_gene_sig = PFCsubregion_enrich_gene_sig[PFCsubregion_enrich_gene_sig$avg_log2FC > 0, ]
write.csv(PFCsubregion_enrich_gene_sig, "./DEGs/PFCsubregion_enrich_gene_sig.DEG.csv")

PFCsubregion_enrich_gene_sig = PFCsubregion_enrich_gene_sig[PFCsubregion_enrich_gene_sig$pct.1 > 0.4, ]

PFCsubregion_enrich_gene_sig_top <- 
  PFCsubregion_enrich_gene_sig %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_L3_1, scale.max = 60,
        features = unique(PFCsubregion_enrich_gene_sig_top$gene)) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  coord_flip() +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") +
  scale_color_gradientn(colours = c(pals::coolwarm(10)[3:10]))
ggsave("../figure_new//GeneExpr_enriched_subregions_dotplot.pdf", width = 5, height = 4)


SpatialPlot(subset(merFISH_integrated, CellType == "Ext"),
            features = c( "Nnat", "Fezf2", "Nr4a1", "Scn4b"),
            stroke = 0, max.cutoff = "q99",
            pt.size.factor = 1.6,
            ncol = 4,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) &
  scale_y_reverse()
  #scale_fill_gradientn(colours = pals::gnuplot(10)) 
ggsave("../figure_new//GeneExpr_enriched_subregions_example.pdf", width = 12, height = 4)

# iSpatial infer
merFISH_integrated_PFC_enrich_iSpatial = readRDS("./RDS/merFISH_integrated_L3_PFC_enrich_iSpatial.rds")

merFISH_integrated_PFC_enrich_iSpatial$ABA_metaRegion = merFISH_integrated@meta.data[colnames(merFISH_integrated_PFC_enrich_iSpatial), "ABA_metaRegion"]
merFISH_integrated_PFC_enrich_iSpatial = subset(merFISH_integrated_PFC_enrich_iSpatial, ABA_metaRegion %in% c("ACAd", "DP", "ILA", "PL", "MOp","MOs"))
merFISH_integrated_PFC_enrich_iSpatial$ABA_metaRegion = factor(merFISH_integrated_PFC_enrich_iSpatial$ABA_metaRegion,
                                                levels = c("DP",  "ILA","PL", "ACAd", "MOs", "MOp"))


Idents(merFISH_integrated_PFC_enrich_iSpatial) = "ABA_metaRegion"
merFISH_integrated_PFC_enrich_iSpatial_gene = FindAllMarkers(merFISH_integrated_PFC_enrich_iSpatial,
                                          pseudocount.use = 1,
                                          group.by = "ABA_metaRegion")
write.csv(merFISH_integrated_PFC_enrich_iSpatial_gene, "./DEGs/PFCsubregion_enrich_iSpatial_gene.DEG.csv")

merFISH_integrated_PFC_enrich_iSpatial_gene_sig = merFISH_integrated_PFC_enrich_iSpatial_gene[merFISH_integrated_PFC_enrich_iSpatial_gene$p_val_adj < 0.05, ]
merFISH_integrated_PFC_enrich_iSpatial_gene_sig = merFISH_integrated_PFC_enrich_iSpatial_gene_sig[merFISH_integrated_PFC_enrich_iSpatial_gene_sig$avg_log2FC > 0, ]
write.csv(merFISH_integrated_PFC_enrich_iSpatial_gene_sig, "./DEGs/PFCsubregion_enrich_iSpatial_gene_sig.DEG.csv")

#merFISH_integrated_PFC_enrich_iSpatial_gene_sig = merFISH_integrated_PFC_enrich_iSpatial_gene_sig[merFISH_integrated_PFC_enrich_iSpatial_gene_sig$pct.1 > 0.5, ]

merFISH_integrated_PFC_enrich_iSpatial_gene_sig_top <- 
  #merFISH_integrated_PFC_enrich_iSpatial_gene_sig[!merFISH_integrated_PFC_enrich_iSpatial_gene_sig$gene %in% rownames(merFISH_integrated), ] %>%
  merFISH_integrated_PFC_enrich_iSpatial_gene_sig %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

DotPlot(object = merFISH_integrated_PFC_enrich_iSpatial, features = unique(merFISH_integrated_PFC_enrich_iSpatial_gene_sig_top$gene)) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  coord_flip() +
  scale_x_discrete(limits=rev) + xlab("") + ylab("") +
  scale_color_gradientn(colours = c(pals::coolwarm(10)[3:10]))
ggsave("../figure_new//GeneExpr_iSpatial_enriched_subregions_dotplot.pdf", width = 5, height = 6)

SpatialPlot(merFISH_integrated_PFC_enrich_iSpatial,
            features = c( "Nnat", "Bcl11b", "Grp", "Npas4", "Rab3b", "Bdnf", "Scn4b"),
            stroke = 0, max.cutoff = "q99",
            pt.size.factor = 1.6,
            ncol = 7,
            images = c("slice.43")) &
  scale_fill_gradientn(colours = c("royalblue2", "gray90", "pink", "red", "darkred")) 


