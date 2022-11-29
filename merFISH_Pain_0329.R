# pain vs control

merFISH_integrated = readRDS(file = "./RDS/merFISH_final.1101.RDS")
merFISH_integrated$subcluster = factor(merFISH_integrated$subcluster, levels = levels(merFISH_integrated_L3$subcluster))

# 
merFISH_integrated_pain = subset(merFISH_integrated, orig.ident %in% c("Pain1", "Pain2", "Pain3", "Pain4", "Pain5", "Pain8", "Pain9", "Pain10", "Pain11"))
# remove Pain3&4
merFISH_integrated_pain = subset(merFISH_integrated, orig.ident %in% c("Pain1", "Pain2", "Pain5", "Pain8", "Pain9", "Pain10", "Pain11"))

merFISH_integrated_pain$subcluster = factor(merFISH_integrated_pain$subcluster, levels = levels(merFISH_integrated_L3$subcluster))

SpatialPlot(merFISH_integrated_pain,group.by = "slice",
            images = "image.22", label = T)

# pain
pain_slice = c(32, 35, 40, 42, 45, 57, 59)
merFISH_integrated_pain$pain = "ctr"
merFISH_integrated_pain$pain[merFISH_integrated_pain$slice %in% pain_slice] = "pain"

# PFC
merFISH_integrated_pain = subset(merFISH_integrated_pain, ABA_PFC == "In")
SpatialPlot(merFISH_integrated_pain, images = "slice.57")

### LR
LR_DEG = lapply(unique(merFISH_integrated_pain$subcluster), function(x){
  merFISH_s = subset(merFISH_integrated_pain, subset = subcluster == x)
  print(x)
  tryCatch({
    degs = FindMarkers(merFISH_s, test.use = "LR", logfc.threshold = 0,
                       latent.vars = "orig.ident", pseudocount.use = 1,
                       ident.1 = "pain", ident.2 = "ctr", group.by = "pain")
    
    gene_cluster = data.frame(gene = rownames(degs),
                              cluster = rep(x, nrow(degs)))
    degs = cbind(degs, gene_cluster)
  }, error = function(e) message(e))
  degs
})
LR_DEG = do.call(rbind, LR_DEG)
LR_DEG_fc = LR_DEG[abs(LR_DEG$avg_log2FC) > log2(1.2) & LR_DEG$p_val_adj < 0.05, ] 
write.csv(LR_DEG, "./DEGs//PFC_pain_LR_DEG_removePain3_4.csv", quote = F)
#write.csv(LR_DEG, "./DEGs//PFC_pain_LR_DEG1101.csv", quote = F)
#LR_DEG = read.csv( "./DEGs//PFC_pain_LR_DEG1101.csv", row.names = 1)

#### MAST
MAST_DEG = lapply(unique(merFISH_integrated_pain$subcluster), function(x){
  merFISH_s = subset(merFISH_integrated_pain, subset = subcluster == x)
  print(x)
  tryCatch({
    degs = FindMarkers(merFISH_s, test.use = "MAST", logfc.threshold = 0,
                       latent.vars = "orig.ident", pseudocount.use = 0.1,
                       ident.1 = "pain", ident.2 = "ctr", group.by = "pain")
    
    gene_cluster = data.frame(gene = rownames(degs),
                              cluster = rep(x, nrow(degs)))
    degs = cbind(degs, gene_cluster)
  }, error = function(e) message(e))
  degs
})
MAST_DEG = do.call(rbind, MAST_DEG)
MAST_DEG_fc = MAST_DEG[abs(MAST_DEG$avg_log2FC) > log2(1.2) & MAST_DEG$p_val_adj < 0.05, ] # 60
write.csv(MAST_DEG, "./DEGs//PFC_pain_MAST_DEG.csv", quote = F)


## plot
LR_DEG.p = LR_DEG
#LR_DEG.p = MAST_DEG
LR_DEG.p$cluster = factor(LR_DEG.p$cluster, levels = levels(Idents(merFISH_integrated_pain)))
LR_DEG.p$DEG = "Non Sig"
LR_DEG.p[LR_DEG.p$avg_log2FC > log2(1) & LR_DEG.p$p_val_adj < 0.05, "DEG"] <- "Up"
LR_DEG.p[LR_DEG.p$avg_log2FC < -log2(1) & LR_DEG.p$p_val_adj < 0.05, "DEG"] <- "Down"
ggplot(LR_DEG.p[LR_DEG.p$cluster %in% levels(LR_DEG.p$cluster)[1:18], ], aes(avg_log2FC, -log10(p_val_adj), color = DEG)) +
  geom_point(size = 1) + ylim(0,10)+
  facet_wrap(~ cluster, nrow = 8) +
  theme_bw()+
  scale_color_manual(values = c("blue", "gray80", "red2"))
ggsave("../figure_new/PFC_Pain_DEG_ExtClusters.pdf", width = 8, height = 12)


library(dplyr)
plotting_df <-
  LR_DEG.p[LR_DEG.p$DEG != "Non Sig", ] %>%
  group_by(cluster, DEG) %>%
  summarise(Freq = n()) %>%
  mutate(Freq = if_else(DEG == "Up", Freq, -Freq))
temp_df <-
  LR_DEG.p[LR_DEG.p$DEG != "Non Sig", ]  %>%
  group_by(cluster) %>%
  summarise(Freq = n())
the_order <- temp_df[order(temp_df$Freq), 1]

ggplot(plotting_df, aes(cluster, y = Freq, group = DEG, fill=DEG)) +
  geom_bar(stat = "identity", width = 0.85) +
  coord_flip() +
  scale_x_discrete(limits = the_order$cluster)+
  cowplot::theme_cowplot() + ylab("# DEGs")+
  scale_fill_manual(values = c( "blue3", "red3"))
ggsave("../figure_new/PFC_Pain_DEG_numbers_allClusters_removePain3_4.pdf", width = 6, height = 8)

#
Pain_DEG_hp = dcast(LR_DEG, cluster~gene, value.var= "avg_log2FC" )
rownames(Pain_DEG_hp) = Pain_DEG_hp$cluster
Pain_DEG_hp = Pain_DEG_hp[,-1]
Pain_DEG_hp = t(Pain_DEG_hp)
Pain_DEG_hp[Pain_DEG_hp > 0.8] <- 0.8
Pain_DEG_hp[Pain_DEG_hp < -0.8] <- -0.8

Pain_DEG_hp1 = Pain_DEG_hp[,1:18]
Pain_DEG_hp1 = Pain_DEG_hp1[rownames(Pain_DEG_hp1) %in% LR_DEG_fc[LR_DEG_fc$cluster %in% colnames(Pain_DEG_hp1), "gene"], ]
pdf("../figure_new/Pain_DEG_fc1.2_heatmap.pdf", 4.5, 8)
pheatmap::pheatmap((Pain_DEG_hp1), scale = "none", cluster_rows = T, cluster_cols = F, border_color = NA,
                   color = colorRampPalette(c("blue3", "white", "white",  "red3"))(50)[1:43], na_col = "gray90")
dev.off()


SpatialPlot(merFISH_integrated_pain, group.by = "subcluster",
            images = paste0("slice.", 34:35))

SpatialPlot(subset(merFISH_integrated_pain, subcluster %in% c("L4/5 IT 2")), 
            label = F,
            features = "Arc",
            slot = "counts",
            pt.size.factor = 3,
            stroke=0.1, images = paste0("slice.", 44:45)) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]), limits=c(0,25) )


# postive control
g = paste0("Hcn", 1:4)
g = paste0("Chrm", 1:5)
g = c("Kcnk5", "Kcnj9")
# mouse 1
SpatialPlot(merFISH_integrated, 
            label = F,
            features = g,
            slot = "counts",
            pt.size.factor = 1.3,
            max.cutoff = "q99",
            stroke=0, images = paste0("slice.", c(44,45,60,59))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  #scale_fill_gradientn(colours =  c("gray90", "gray90", "red", "red4") )
#mouse 2
SpatialPlot(merFISH_integrated, 
            label = F,
            features = g,
            #slot = "counts",
            pt.size.factor = 1.3,
            max.cutoff = "q99",
            stroke=0, images = paste0("slice.", c(43,42,44,45))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  #scale_fill_gradientn(colours =  c("gray90", "gray90", "red", "red4") )



g = "Hcn1"
SpatialFeaturePlot_CCF(merFISH_integrated_pain,
                       feature = g,
                       slice = 44,
                       max.cut = 30,
                       ABA_coronal = 37,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)
SpatialFeaturePlot_CCF(merFISH_integrated_pain,
                            feature = g,
                            slice = 45,
                       max.cut = 30,
                            ABA_coronal = 37,
                            ABA_color = "no",
                            pt.size = 0.3,
                            plot.axis = F)

SpatialFeaturePlot_CCF(merFISH_integrated_pain,
                       feature = "Arc",
                       slice = 34,
                       max.cut = 25,
                       ABA_coronal = 34,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)

SpatialFeaturePlot_CCF(merFISH_integrated_pain,
                       feature = "Arc",
                       slice = 35,
                       max.cut = 25,
                       ABA_coronal = 34,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)

SpatialPlot(merFISH_integrated, 
            label = F,
            features = "Junb",
            slot = "counts",
            pt.size.factor = 1,
            stroke=0, images = paste0("slice.", 46:48)) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() 
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]), limits=c(0,25) )

SpatialFeaturePlot_CCF(merFISH_integrated,
                       feature = "Junb",
                       slice = 46,
                       max.cut = 30,
                       ABA_coronal = 32,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)
SpatialFeaturePlot_CCF(merFISH_integrated,
                       feature = "Junb",
                       slice = 47,
                       max.cut = 30,
                       ABA_coronal = 32,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)
SpatialFeaturePlot_CCF(merFISH_integrated,
                       feature = "Junb",
                       slice = 48,
                       max.cut = 30,
                       ABA_coronal = 32,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)



SpatialPlot(merFISH_integrated_pain, 
            label = F,
            features = "Arc",
            slot = "counts",
            pt.size.factor = 1,
            stroke=0, images = paste0("slice.", c(43,42))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]), limits=c(0,25) )



#
SpatialDimPlot(subset(merFISH_integrated_pain, CellType == "Ext"),
                   pt.size.factor = 1,
                   ncol = 2,
                   stroke=0, images = paste0("slice.", c(33,32,34,35,37,36,39,38,41,40,43,42,44,45))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_color_manual(values = L2_color) & NoAxes()

SpatialDimPlot(subset(merFISH_integrated_pain, CellType == "Ext"),
               pt.size.factor = 1,
               ncol = 2,
               stroke=0, images = paste0("slice.", c(44,45))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_color_manual(values = L2_color) & NoAxes()


SpatialDimPlot(subset(merFISH_integrated_pain, CellType == "Ext"),
               pt.size.factor = 1,
               ncol = 2,
               stroke=0, images = paste0("slice.", c(33,32,34,35,37,36,39,38,41,40,43,42,44,45))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() 


SpatialFeaturePlot(merFISH_integrated_pain, 
            features = "Htr2c",
            #slot = "counts",
            pt.size.factor = 1,
            ncol = 2,
            stroke=0, images = paste0("slice.", c(37,36,39,38,41,40,43,42,44,45))) &
  theme(aspect.ratio = 1) & scale_x_reverse() & coord_flip() &
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]), limits=c(0,4) )

SpatialFeaturePlot(merFISH_integrated_pain, 
                     features = "Arc",
                     #slot = "counts",
                     pt.size.factor = 1,
                     ncol = 2,
                     stroke=0, images = paste0("slice.", c(60,59))) &
    theme(aspect.ratio = 1) & scale_x_reverse() 
  scale_fill_gradientn(colours =  c(rainbow(10)[3:10]), limits=c(0,3) )
  
  SpatialFeaturePlot(merFISH_integrated, 
                     features = c("Arc", "Fosb", "Npas4", "Grin2a"),
                     #slot = "counts",
                     max.cutoff = "q99",
                     pt.size.factor = 1,
                     images = c("image.22"),
                     ncol = 7,
                     stroke=0) &
    theme(aspect.ratio = 1) & scale_x_reverse()  &
    scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  
  SpatialFeaturePlot(merFISH_integrated, 
                     features = c("nCount_RNA"),
                     #slot = "counts",
                     pt.size.factor = 1,
                     images = "image.22",
                     ncol = 1,
                     stroke=0) &
    theme(aspect.ratio = 1) & scale_x_reverse()  &
    scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  
  
  # pain 3 & 4
  SpatialFeaturePlot(merFISH_integrated, 
                     features = NAGs,
                     #slot = "counts",
                     max.cutoff = "q99",
                     pt.size.factor = 1,
                     images = "image.13",
                     ncol = 6,
                     stroke=0) &
    theme(aspect.ratio = 1.5) & scale_x_reverse()  &
    scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  
  SpatialDimPlot(merFISH_integrated, 
                     #slot = "counts",
                     pt.size.factor = 1,
                     images = c("image.12", "image.13"),
                     stroke=0) &
    theme(aspect.ratio = 1.5) & scale_x_reverse()  
    scale_fill_gradientn(colours =  c(rainbow(10)[3:10]))
  

VlnPlot(merFISH_integrated_pain, features = c("Kcnh7"), split.by = c("pain", "orig.ident"), pt.size = -1)

#

# neruonal active genes
#NAGs = read.table("genelist/Neuron_activate_gene.txt")$V1
NAGs = read.table("genelist/rPRG.genes.txt")$V1 # list from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5934296/#SD1
NAGs = NAGs[NAGs%in% rownames(merFISH_integrated_pain)]

merFISH_integrated_pain_ext = subset(merFISH_integrated_pain, CellType == "Ext")
merFISH_integrated_pain_ext$mice = 1
merFISH_integrated_pain_ext$mice[merFISH_integrated_pain_ext$slice %in% 36:45] <- 2
merFISH_integrated_pain_ext$mice[merFISH_integrated_pain_ext$slice %in% 57:60] <- 3

merFISH_integrated_pain_ext = NormalizeData(merFISH_integrated_pain_ext)
merFISH_integrated_pain_ext = ScaleData(merFISH_integrated_pain_ext, vars.to.regress = c("orig.ident"))

mean.exp <- colMeans(x = merFISH_integrated_pain_ext@assays$RNA@scale.data[NAGs, ], na.rm = TRUE)
merFISH_integrated_pain_ext@meta.data$active.score <- mean.exp
FeaturePlot(merFISH_integrated_pain_ext, features = "active.score", raster=FALSE, split.by = "pain") &
  NoAxes()  &
  scale_colour_gradientn(colours = c("lightblue", "yellow", "red", "black"))

smoothScatter(merFISH_integrated_pain_ext$nCount_RNA, merFISH_integrated_pain_ext$active.score)

ggplot(merFISH_integrated_pain_ext@meta.data, aes(x = active.score)) +
  geom_density() + facet_grid(slice ~.)


VlnPlot(merFISH_integrated_pain_ext, features = "active.score", split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
          pt.size = 0, slot = "scale.data", ncol=1) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=median, geom="point", size=8, color="darkred")
ggsave("../figure_new/PFC_Pain_Neuron_activeScore_allClusters.pdf", width = 8, height = 3.5)

p1 = VlnPlot(subset(merFISH_integrated_pain_ext, mice == 1),
        features = "active.score", split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
        pt.size = 0, slot = "scale.data", ncol=1) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=median, geom="point", size=8, color="darkred")
p2 = VlnPlot(subset(merFISH_integrated_pain_ext, mice == 2),
             features = "active.score", split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
             pt.size = 0, slot = "scale.data", ncol=1) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=median, geom="point", size=8, color="darkred")
p3 = VlnPlot(subset(merFISH_integrated_pain_ext, mice == 3),
             features = "active.score", split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
             pt.size = 0, slot = "scale.data", ncol=1) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=median, geom="point", size=8, color="darkred")
p1/p2/p3
ggsave("../figure_new/PFC_Pain_Neuron_activeScore_allClusters_3mice.pdf", width = 8, height = 10)


VlnPlot(merFISH_integrated_pain_ext, features = NAGs, split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
        pt.size = 0, slot = "data", ncol=1) &
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=mean, geom="point", size=8, color="darkred")
ggsave("../figure_new/PFC_Pain_Neuron_IEGs_allClusters.pdf", width = 8, height = 16)

merFISH_integrated_pain_ext_score = merFISH_integrated_pain_ext@meta.data
for(i in levels(Idents(merFISH_integrated_pain_ext))){
  print(i)
  merFISH_integrated_pain_ext_score1 = merFISH_integrated_pain_ext_score[merFISH_integrated_pain_ext_score$subcluster == i, ]
  p = wilcox.test(merFISH_integrated_pain_ext_score1[merFISH_integrated_pain_ext_score1$pain == "pain", "active.score"],
              merFISH_integrated_pain_ext_score1[merFISH_integrated_pain_ext_score1$pain == "ctr", "active.score"])$p.value
  print(p)
}

merFISH_integrated_pain_ext_score$pain_ctr = paste0(merFISH_integrated_pain_ext_score$orig.ident, merFISH_integrated_pain_ext_score$pain)
merFISH_integrated_pain_ext_score$sample_cluster = paste0(merFISH_integrated_pain_ext_score$orig.ident, merFISH_integrated_pain_ext_score$subcluster)

a = merFISH_integrated_pain_ext_score %>% group_by(orig.ident, pain, subcluster, mice) %>%
  summarise_all(mean)

ggplot(as.data.frame(a), aes(pain, active.score)) +
  geom_boxplot(aes(fill = pain), outlier.size = -1) +
  geom_point(aes(shape = pain, group = orig.ident, color = factor(mice)), position = position_dodge(0.3), size = 2) +
  geom_line(aes(color = factor(mice), group = orig.ident), position = position_dodge(0.1)) +
  #stat_summary(fun.data = mean_cl_boot,shape = 18, size = 1, color = "orange") +
  facet_wrap(~subcluster, ncol = 18) +
  scale_fill_manual(values = c("gray70", "red2")) +
  scale_color_manual(values = c("blue3", "orange2", "green4")) +
  cowplot::theme_cowplot()
ggsave("../figure_new/PFC_Pain_Neuron_activeScore_allClusters_pseudoBulk_removePain3_4.pdf", width = 15, height = 4)

a = as.data.frame(a)
pain_effect = a[a$pain == "pain", ]
pain_effect$active.score = pain_effect$active.score - a[a$pain != "pain", "active.score"]  
pain_effect1 <-  pain_effect %>%
  group_by(subcluster) %>%
  summarise( 
    n=n(),
    mean=mean(active.score),
    sd=sd(active.score)
  ) %>%
  mutate( se=sd/sqrt(n))

ggplot(pain_effect1, aes(subcluster, mean)) +
  geom_bar(stat="identity", aes(fill = subcluster), outlier.size = -1) +
  scale_fill_manual(values = L2_color) +
  geom_point(data = pain_effect, aes(subcluster, active.score, color = factor(mice)), position = position_dodge(0.3), size = 2) +
  geom_errorbar( aes(x=subcluster, ymin=mean-se, ymax=mean+se), width=0.3, colour="red2", alpha=1, size=1) +
  scale_color_manual(values = c("blue3", "orange2", "green4")) +
  cowplot::theme_cowplot() +
  theme(legend.position = "None")

ggplot(pain_effect, aes(subcluster, active.score)) +
  geom_boxplot(aes(fill = subcluster),outlier.size = -1) +
  geom_point(aes(color = factor(mice)) , position = position_dodge(0.3), size = 2) +
  scale_fill_manual(values = L2_color) +
  scale_color_manual(values = c("blue3", "orange2", "green4")) +
  ylab("Normalized IEG score:\n Pain - Ctr") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  theme(legend.position = "None")
ggsave("../figure_new/PFC_Pain_Neuron_activeScore_Different_allClusters_pseudoBulk_removePain3_4.pdf", width = 10, height = 4)


ggplot(pain_effect, aes(subcluster, active.score)) +
  geom_bar(data = pain_effect1, aes(subcluster, mean, fill = subcluster), stat="identity", outlier.size = -1) +
  geom_point(position = position_dodge(0.3), size = 2) +
  scale_fill_manual(values = L2_color) +
  stat_summary(fun.data = mean_cl_boot,shape = 18, size = 1, color = "orange") +
  cowplot::theme_cowplot()



g = ggplot(as.data.frame(a), aes(pain, active.score)) +
  #geom_boxplot(aes(fill = pain), outlier.size = -1) +
  geom_point(aes(color = orig.ident), position = position_dodge(0.3), size = 2) +
  geom_line(aes(color = orig.ident, group = orig.ident), position = position_dodge(0.3)) +
  stat_summary(fun.data = mean_cl_boot,shape = 18, size = 1, color = "orange") +
  facet_wrap(~subcluster, ncol = 18) +
  theme_cowplot()
plotly::ggplotly(g)


for(i in levels(Idents(merFISH_integrated_pain_ext))){
  print(i)
  a1 = as.data.frame(a)
  a1 = a1[a1$subcluster == i, ]
  p = t.test(a1[a1$pain == "pain", "active.score"],
                  a1[a1$pain == "ctr", "active.score"],
                  paired = T)$p.value
  print(p)
}



# select cluster to plot
VlnPlot(subset(merFISH_integrated_pain_ext, subcluster %in% c("L2/3 IT 4", "L5 ET 1", "L5 IT 2","L5 IT 3")),
        features = "active.score", split.by = "pain", split.plot =F,cols = c( "gray70",  "red2"),
        y.max = 3,
        pt.size = 0) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=mean, geom="point", size=8, color="darkred")


FeaturePlot(merFISH_integrated_pain, features = c( "Arc"), raster=FALSE, split.by = "pain",
            slot = "scale.data", max.cutoff = "q99") &
  NoAxes()  &
  scale_colour_gradientn(colours =  c(rainbow(10)[3:10]))



# subregion
merFISH_integrated_pain_ALL = subset(merFISH_integrated, orig.ident %in% c("Pain1", "Pain2", "Pain3", "Pain4", "Pain5", "Pain8", "Pain9", "Pain10", "Pain11"))
# remove 3 4
merFISH_integrated_pain_ALL = subset(merFISH_integrated, orig.ident %in% c("Pain1", "Pain2", "Pain5", "Pain8", "Pain9", "Pain10", "Pain11"))

# pain
merFISH_integrated_pain_ALL$pain = "ctr"
merFISH_integrated_pain_ALL$pain[merFISH_integrated_pain_ALL$slice %in% pain_slice] = "pain"

merFISH_integrated_pain_ALL_ext = subset(merFISH_integrated_pain_ALL, CellType == "Ext")

merFISH_integrated_pain_ALL_ext = ScaleData(merFISH_integrated_pain_ALL_ext, vars.to.regress = c( "orig.ident"))

mean.exp <- colMeans(x = merFISH_integrated_pain_ALL_ext@assays$RNA@scale.data[NAGs, ], na.rm = TRUE)
merFISH_integrated_pain_ALL_ext@meta.data$active.score <- mean.exp

NAG_subregion = merFISH_integrated_pain_ALL_ext@meta.data[,c("ABA_metaRegion","ABA_regions", "L3_cluster", "orig.ident", "pain", "active.score")]
NAG_subregion = NAG_subregion[NAG_subregion$ABA_metaRegion %in% c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm", "MOs", "TTd"), ]

#NAG_subregion$pain = paste( NAG_subregion$pain, NAG_subregion$orig.ident)

NAG_subregion1 = reshape2::dcast(NAG_subregion, ABA_metaRegion~pain, mean)
NAG_subregion1 = na.omit(NAG_subregion1)
rownames(NAG_subregion1) = NAG_subregion1$ABA_metaRegion
NAG_subregion1 = NAG_subregion1[,-1]
pdf("../figure_new/NAG_PFC_subregion_remove3_4.pdf", 3, 6)
pheatmap::pheatmap(NAG_subregion1, cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   color = pals::coolwarm(50))
dev.off()

NAG_subregion2 = reshape2::dcast(NAG_subregion, ABA_regions~pain, mean)
NAG_subregion2 = na.omit(NAG_subregion2)
rownames(NAG_subregion2) = NAG_subregion2$ABA_regions
NAG_subregion2 = NAG_subregion2[,-1]
pdf("../figure_new/NAG_PFC_subregion_withLayer_remove3_4.pdf", 5, 12)
pheatmap::pheatmap(NAG_subregion2, cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   color = pals::coolwarm(30))
dev.off()

#
NAG_subregion = merFISH_integrated_pain_ALL_ext@meta.data[,c("ABA_metaRegion","ABA_regions", "subcluster", "orig.ident", "pain", "active.score")]
NAG_subregion = NAG_subregion[NAG_subregion$ABA_metaRegion %in% c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm"), ]

NAG_subregion$subcluster = paste(  NAG_subregion$ABA_metaRegion, NAG_subregion$subcluster)
rm_small_cells = table(NAG_subregion$subcluster)
NAG_subregion = NAG_subregion[NAG_subregion$subcluster %in% names(rm_small_cells[rm_small_cells >= 100]), ]

#NAG_subregion$pain = paste( NAG_subregion$pain, NAG_subregion$orig.ident)

NAG_subregion3 = reshape2::dcast(NAG_subregion[NAG_subregion$ABA_metaRegion %in% c("ACAd", "PL"), ], subcluster~pain, mean)
rownames(NAG_subregion3) = NAG_subregion3$subcluster
NAG_subregion3 = NAG_subregion3[,-1]
pdf("../figure_new/NAG_PFC_subregion_subcluster_PL_ACAd_remove3_4.pdf", 4, 10)
pheatmap::pheatmap((NAG_subregion3), cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   color = pals::coolwarm(30))
dev.off()

pdf("../figure_new/NAG_PFC_subregion_subcluster_selected.pdf", 3, 8)
pheatmap::pheatmap((NAG_subregion3[abs(NAG_subregion3$pain - NAG_subregion3$ctr) >= 0.1, ]), cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   color = pals::viridis(50))
dev.off()

#
merFISH_integrated_pain_ext$active.score.scale = merFISH_integrated_pain_ext$active.score
merFISH_integrated_pain_ext$active.score.scale[merFISH_integrated_pain_ext$active.score < 0.3 ] <- 0 
merFISH_integrated_pain_ext$active.score.scale[merFISH_integrated_pain_ext$active.score >= 0.3 & merFISH_integrated_pain_ext$active.score < 1 ] <- 1
merFISH_integrated_pain_ext$active.score.scale[merFISH_integrated_pain_ext$active.score >= 1 ] <- 2

SpatialPlot(merFISH_integrated_pain_ext,
            features = "active.score",
            max.cutoff = 2.5, stroke = 0,
            images = c("slice.37", "slice.36")) &
  scale_fill_gradientn(colours = c(rainbow(10)[4:10]))


SpatialPlot(merFISH_integrated_pain_ext,
            features = "active.score.scale",
            max.cutoff = 2.5, stroke = 0,
            images = c("slice.37", "slice.36")) &
  scale_fill_gradientn(colours = c( "gray90", "darkorange", "red2"))

pdf("../figure_new/NAG_score_spatial_CCF.pdf", 12, 12)
SpatialFeaturePlot_CCF(merFISH_integrated_pain_ext,
                       feature = "active.score",
                       slice = 44,
                       max.cut = 2.5,
                       ABA_coronal = 37,
                       ABA_color = "No",
                       pt.size = 0.3,
                       plot.axis = F)

SpatialFeaturePlot_CCF(merFISH_integrated_pain_ext,
                       feature = "active.score",
                       slice = 45,
                       max.cut = 2.5,
                       ABA_coronal = 37,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)
dev.off()

merFISH_integrated_pain_ext$active.score1 = merFISH_integrated_pain_ext$active.score
merFISH_integrated_pain_ext$active.score1[merFISH_integrated_pain_ext$active.score1 < 0.5] <- 0

pdf("../figure_new/NAG_score_spatial_CCF1.pdf", 12, 12)
SpatialFeaturePlot_CCF(merFISH_integrated_pain_ALL_ext,
                       feature = "active.score1",
                       slice = 60,
                       max.cut = 2,
                       ABA_coronal = 39,
                       ABA_color = "no",
                       pt.size = 0.3,
                       plot.axis = F)
SpatialFeaturePlot_CCF(merFISH_integrated_pain_ext,
                       feature = "active.score1",
                       slice = 59,
                       max.cut = 2,
                       ABA_coronal = 37,
                       ABA_color = "No",
                       pt.size = 0.3,
                       plot.axis = F)
dev.off()

# 
merFISH_Fentanyl = subset(merFISH_integrated, orig.ident == "Fentanyl")
#merFISH_Fentanyl = subset(merFISH_integrated, orig.ident == "Fentanyl" & ABA_PFC == "In")

merFISH_Fentanyl = ScaleData(merFISH_Fentanyl, vars.to.regress = c("abs_volume"))

merFISH_Fentanyl$treatment = "Ctr"
merFISH_Fentanyl@meta.data[merFISH_Fentanyl$slice == 46, "treatment"] <- "Morphine"
merFISH_Fentanyl@meta.data[merFISH_Fentanyl$slice == 47, "treatment"] <- "Fentanyl"

mean.exp <- colMeans(x = merFISH_Fentanyl@assays$RNA@scale.data[NAGs, ], na.rm = TRUE)
merFISH_Fentanyl@meta.data$active.score <- mean.exp
FeaturePlot(merFISH_Fentanyl, features = "active.score", raster=FALSE, split.by = "slice") &
  NoAxes()  &
  scale_colour_gradientn(colours = c("lightblue", "yellow", "red", "black"))


VlnPlot(subset(merFISH_Fentanyl, CellType == "Ext"), features = "active.score", split.by = "treatment", split.plot =F,
        cols = c( "gray80", "orange", "red2"),
        pt.size = 0, slot = "scale.data", ncol=1) +
  stat_summary(aes(color = "pain"), position = position_dodge(1),
               shape = 95, fun.y=mean, geom="point", size=8, color="darkred")
ggsave("../figure_new/PFC_Fentanyl_Neuron_activeScore_allClusters.pdf", width = 8, height = 3.5)



###### remove "Pain1" "Pain2" samples. They are not exactly match
# 
merFISH_integrated_pain = subset(merFISH_integrated_L3, orig.ident %in% c("Pain3", "Pain4", "Pain5", "Pain8", "Pain9"))

# pain
pain_slice = c(36, 38, 40, 42, 45)
merFISH_integrated_pain$pain = "ctr"
merFISH_integrated_pain$pain[merFISH_integrated_pain$slice %in% pain_slice] = "pain"


LR_DEG = lapply(unique(merFISH_integrated_pain$subcluster), function(x){
  merFISH_s = subset(merFISH_integrated_pain, subset = subcluster == x)
  print(x)
  tryCatch({
    degs = FindMarkers(merFISH_s, test.use = "LR", logfc.threshold = 0,
                       latent.vars = "orig.ident", pseudocount.use = 0.1,
                       ident.1 = "pain", ident.2 = "ctr", group.by = "pain")
    
    gene_cluster = data.frame(gene = rownames(degs),
                              cluster = rep(x, nrow(degs)))
    degs = cbind(degs, gene_cluster)
  }, error = function(e) message(e))
  degs
})
LR_DEG = do.call(rbind, LR_DEG)
LR_DEG_fc = LR_DEG[abs(LR_DEG$avg_log2FC) > log2(1.2) & LR_DEG$p_val_adj < 0.05, ] # 143
write.csv(LR_DEG, "./DEGs/pain_5samples_LR_DEG.csv", quote = F)


# plot
LR_DEG.p = LR_DEG
LR_DEG.p$DEG = "Non Sig"
LR_DEG.p[LR_DEG.p$avg_log2FC > log2(1) & LR_DEG.p$p_val_adj < 0.05, "DEG"] <- "Up"
LR_DEG.p[LR_DEG.p$avg_log2FC < -log2(1) & LR_DEG.p$p_val_adj < 0.05, "DEG"] <- "Down"
ggplot(LR_DEG.p, aes(avg_log2FC, -log10(p_val_adj), color = DEG)) +
  geom_point(size = 1) + ylim(0,10)+
  facet_wrap(~ cluster, nrow = 8) +
  theme_bw()+
  scale_color_manual(values = c("blue", "gray80", "red2"))

ggsave("../figure_new/pain_5samples_DEG_allClusters.pdf", width = 10, height = 8)


library(dplyr)
plotting_df <-
  LR_DEG.p[LR_DEG.p$DEG != "Non Sig", ] %>%
  group_by(cluster, DEG) %>%
  summarise(Freq = n()) %>%
  mutate(Freq = if_else(DEG == "Up", Freq, -Freq))
temp_df <-
  LR_DEG.p[LR_DEG.p$DEG != "Non Sig", ]  %>%
  group_by(cluster) %>%
  summarise(Freq = n())
the_order <- temp_df[order(temp_df$Freq), 1]

ggplot(plotting_df, aes(cluster, y = Freq, group = DEG, fill=DEG)) +
  geom_bar(stat = "identity", width = 0.85) +
  coord_flip() +
  scale_x_discrete(limits = the_order$cluster)+
  theme_cowplot() + ylab("# DEGs")+
  scale_fill_manual(values = c( "blue3", "red3"))
ggsave("../figure_new/pain_5samples_DEG_numbers_allClusters.pdf", width = 6, height = 6)

#
Pain_DEG_hp = dcast(LR_DEG, cluster~gene, value.var= "avg_log2FC" )
rownames(Pain_DEG_hp) = Pain_DEG_hp$cluster
Pain_DEG_hp = Pain_DEG_hp[,-1]
Pain_DEG_hp = t(Pain_DEG_hp)
Pain_DEG_hp[Pain_DEG_hp > 0.8] <- 0.8
Pain_DEG_hp[Pain_DEG_hp < -0.8] <- -0.8

Pain_DEG_hp1 = Pain_DEG_hp[,1:18]
Pain_DEG_hp1 = Pain_DEG_hp1[rowSums(is.na(Pain_DEG_hp1)) < 12, ]
Pain_DEG_hp1 = Pain_DEG_hp1[rownames(Pain_DEG_hp1) %in% LR_DEG_fc[LR_DEG_fc$cluster %in% colnames(Pain_DEG_hp1), "gene"], ]
pdf("../figure_new/Ppain_5samples_DEG_fc1.2_heatmap.pdf", 4.5, 8)
pheatmap::pheatmap(Pain_DEG_hp1, scale = "none", cluster_rows = T, cluster_cols = T, border_color = NA,
                   color = colorRampPalette(c("blue3", "white", "red3"))(50), na_col = "gray90")
dev.off()



#

# cFos
cFos_percent = read.table("./DEGs/Pain_ctr_cFos_percent.txt", header = T)
cFos_percent = melt(cFos_percent)


ggplot(cFos_percent, aes(variable, value, color = variable)) +
  geom_point( position = position_dodge(1), size = 2) +
  stat_summary(fun.data = mean_cl_boot,shape = 18, geom="errorbar", width = .2, color = "red") +
  stat_summary(fun.y=mean, geom="point", shape = 18, size = 2, color="red") +
  #scale_fill_manual(values = L2_color) +
  scale_color_manual(values = c("blue3", "orange2", "green4")) +
  ylab("Normalized IEG score:\n Pain - Ctr") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  theme(legend.position = "None")

df.summary <- cFos_percent %>%
  group_by(variable) %>%
  summarise(
    sd = sd(value, na.rm = TRUE),
    value = mean(value)
  )
df.summary

ggplot(cFos_percent, aes(variable, value, fill = variable)) +
  geom_col(data = df.summary, width = 0.7) +
  geom_jitter( position = position_jitter(0.2), size=2, color = "black") + 
  geom_errorbar( aes(ymin = value-sd, ymax = value+sd), 
                 data = df.summary, width = 0.2, size=1) +
  scale_fill_manual(values = c("gray70", "red2")) +
  theme_cowplot() + ylab("Percentage of cFos+ in Pou3f1+ cells") + xlab("")
ggsave("../figure_new/cFos_percentage_Pou3f1.pdf", width = 4, height = 5)

wilcox.test(cFos_percent[cFos_percent$variable == "Control", 2],
            cFos_percent[cFos_percent$variable == "Pain", 2],)
