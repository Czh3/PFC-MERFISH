# predict excitatory neuron projecting

load("../scRNA/GSE161936_Lui_Nguyen_Seuratobjects.RData")

integrated_Exc = subset(merFISH_integrated, CellType == "Ext")

DefaultAssay(Fig3Rbp4retro) <- "RNA"
Fig3Rbp4retro@assays$integrated = NULL

#
integrated_Exc_sample = subset(integrated_Exc, slice %in% c("44"))
SpatialPlot(integrated_Exc_sample, images = "slice.44")
SpatialPlot(merFISH_integrated, features = "nCount_RNA",
            pt.size.factor = 1,
            stroke = 0, images = "image.10") &
  #scale_fill_gradientn(colours =  as.vector(pals::gnuplot(20))) &
  scale_fill_gradientn(colours =  c("gray90", "gray90", "red2")) &
  theme(panel.background = element_rect(fill = 'black', color = 'black'),
        panel.grid = element_line(color = 'black'))  &
  theme(aspect.ratio = 2) 
  

#
genes_select = intersect(rownames(integrated_Exc_sample), rownames(Fig3Rbp4retro))

Idents(Fig3Rbp4retro) = "Cre.line"
Fig3Rbp4retro = DietSeurat(Fig3Rbp4retro)
Fig3Rbp4retro = subset(Fig3Rbp4retro, Tissue.type != "OFC")
#Fig3Rbp4retro <- subset(Fig3Rbp4retro, downsample = 300)
integrated_Exc_sample = DietSeurat(integrated_Exc_sample)

# normalize
Fig3Rbp4retro = NormalizeData(Fig3Rbp4retro[genes_select, ])
integrated_Exc_sample = NormalizeData(integrated_Exc_sample[genes_select,])

#
integrated_Exc = merge(Fig3Rbp4retro, integrated_Exc_sample)
integrated_Exc@meta.data$orig.ident[is.na(integrated_Exc@meta.data$orig.ident)] <- "scRNA"
integrated_Exc@meta.data$tech = ifelse(is.na(integrated_Exc@meta.data$abs_volume), "scRNA", "merFISH")

integrated_Exc = NormalizeData(integrated_Exc)
VariableFeatures(integrated_Exc) <- genes_select

integrated_Exc = ScaleData(integrated_Exc, vars.to.regress = c("nFeature_RNA", "tech"), verbose = FALSE)
integrated_Exc = RunPCA(integrated_Exc, npcs = 50,  verbose = FALSE)

DimPlot(object = integrated_Exc, reduction = "pca", pt.size = .1, group.by = "tech") 

# harmony
require(harmony)
integrated_Exc <- RunHarmony(
  object = integrated_Exc,
  group.by.vars = c('tech'),
  dims.use = 1:30,
  max.iter.harmony = 15,
  plot_convergence = TRUE,
) # regress batch effect & tech


DimPlot(object = integrated_Exc, reduction = "harmony", pt.size = .1, group.by = "tech") 

integrated_Exc <- RunUMAP(integrated_Exc, reduction = 'harmony', dims = 1:30)

p1 = DimPlot(subset(integrated_Exc, subset = tech == "scRNA"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "Cre.line", label = F) &
  NoAxes()
p2 = DimPlot(subset(integrated_Exc, subset = tech == "merFISH"), reduction = "umap",
             repel = TRUE, shuffle = TRUE, group.by = "subcluster", label = F) &
  NoAxes()
p1|p2
ggsave("../figure_new/Projecting_cobedding_UMAP_split_new.pdf", width = 16, height = 8)


DimPlot(object = integrated_Exc, pt.size = 0.5, group.by = "tech",
        shuffle = T) + NoAxes()
ggsave("../figure_new/Projecting_cobedding_UMAP_new.pdf", width = 8, height = 8)


saveRDS(integrated_Exc, "./RDS/integrated_Exc_projecting.RDS")
integrated_Exc = readRDS("./RDS/integrated_Exc_projecting.RDS")


DimPlot(integrated_Exc, reduction = "umap", group.by = "tech", pt.size = 0.1, 
        repel = TRUE, shuffle = TRUE) & NoAxes()


ElbowPlot(integrated_Exc, ndims = 50, reduction = "pca")

# train model
library(e1071)

cells_embed = integrated_Exc@reductions$harmony@cell.embeddings[, 1:30]

cells_embed_trn = as.data.frame(cells_embed[rownames(cells_embed) %in% colnames(Fig3Rbp4retro), ])
cells_embed_trn$projecting = as.factor(integrated_Exc$Cre.line[rownames(cells_embed_trn)])

set.seed(1)
train_cells = sample(rownames(cells_embed_trn), nrow(cells_embed_trn) * 0.5)
# training data
cells_embed_train = cells_embed_trn[train_cells, ]
# test data
cells_embed_test = cells_embed_trn[!rownames(cells_embed_trn) %in% train_cells, ]
# predict MERFISH data
cells_embed_predict= as.data.frame(cells_embed[rownames(cells_embed) %in% colnames(integrated_Exc_sample), ])


set.seed(1)
tune.out=tune(svm , projecting~.,data=cells_embed_train, kernel ="radial", scale = F,
              ranges =list(cost=c(1, 5, 10, 20),
                           gamma=c(0.001, 0.01, 0.1)))

summary(tune.out)
# cost gamma
#   10  0.01
bestmod = tune.out$best.model

# train
svmfit = svm(projecting ~ ., data = cells_embed_train, kernel = "radial", scale = F,
             cost = 10, gamma = 0.01, probability=T)

# test
pred <- predict(svmfit, cells_embed_test[, -31], probability=TRUE)
table(pred, cells_embed_test$projecting)
#mean(pred == cells_embed_test$projecting)


library(pROC)
roc.multi = multiclass.roc(cells_embed_test[, 31], attr(pred, "probabilities"))
auc(roc.multi) # 0.9134

pred_prob = attr(pred, "probabilities")

project_col = c("#A3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032","gray60", "navy")


pdf("../figure_new/projecting_ROC_new.pdf", 6, 6)
plot.roc(ifelse(cells_embed_test$projecting == "Cav_cPFC", 1, 0), pred_prob[,3], col = project_col[1])

lapply(2:7, function(x){
  target = levels(cells_embed_test$projecting)[x]
  gt = ifelse(cells_embed_test$projecting == target, 1, 0)
  plot.roc(gt, pred_prob[,target], col = project_col[x], add = T)
})
legend(.4, .5, legend = levels(cells_embed_test$projecting),
       col=project_col, lty=1, lwd = 2, box.lty=0)
dev.off()




# predict
pred <- predict(svmfit, cells_embed_predict, probability=TRUE)

integrated_Exc$predict_projecting = integrated_Exc$Cre.line
integrated_Exc$predict_projecting[names(pred)] <- as.character(pred)

DimPlot(integrated_Exc, reduction = "umap", group.by = "predict_projecting", split.by = "tech", pt.size = 0.1, 
        repel = TRUE, shuffle = TRUE) & NoAxes()
DimPlot(integrated_Exc[, integrated_Exc$predict_projecting != "Rbp4"], reduction = "umap", group.by = "predict_projecting", pt.size = 0.1, 
        repel = TRUE, shuffle = TRUE) & NoAxes()


project_col = c("gray90","#A3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "navy")
integrated_Exc$predict_projecting = factor(integrated_Exc$predict_projecting,
                                           levels = c("Rbp4", "Cav_cPFC", "Cav_Dstriatum",
                                                      "Cav_Hypo", "Cav_NAc", "Cav_PAG", "Retro_Amyg"))

SpatialPlot(integrated_Exc,
            group.by = "predict_projecting", pt.size = 1, stroke =  0,
            images = "slice.44") &
  theme(aspect.ratio = 1) &
  scale_fill_manual(values = project_col)
ggsave("../figure_new/Projecting_spatial_slice44.pdf", width = 8, height = 8)

SpatialPlot(integrated_Exc,
            group.by = "predict_projecting", pt.size = 1, stroke =  0,
            images = "slice.44") &
  theme(aspect.ratio = 1) &
  scale_fill_manual(values = project_col)


integrated_Exc_plot = integrated_Exc@meta.data
integrated_Exc_plot = integrated_Exc_plot[!is.na(integrated_Exc_plot$slice) & integrated_Exc_plot$slice == 44,]

ggplot(integrated_Exc_plot, aes(reloc_x, reloc_y, color = predict_projecting)) +
  geom_point(size = ifelse(integrated_Exc_plot$predict_projecting=="Rbp4", 0.3, 0.5)) +
  scale_color_manual(values = project_col) +
  theme_void() + coord_fixed() +
  theme(aspect.ratio = 1)
ggsave("../figure_new/Projecting_spatial_slice44.pdf", width = 8, height = 8)



Idents(integrated_Exc) = "predict_projecting"

SpatialDimPlot(integrated_Exc,
               cells.highlight = CellsByIdentities(object = integrated_Exc)[2:7], 
               images = "slice.44",
               crop = T,
               stroke=0,
               pt.size.factor = 1.2,
               ncol = 3,
               cols.highlight = c("#DE2D26", "grey90"),
               facet.highlight = TRUE)
mySpatialDimPlot(integrated_Exc,
                 cells.highlight.size.factor = c(1.5, 1.8),
                 cells.highlight = CellsByIdentities(object = integrated_Exc)[2:7], 
                 images = "slice.44",
                 crop = T,
                 stroke=0,
                 ncol = 3,
                 cols.highlight = c("#DE2D26", "grey90"),
                 facet.highlight = TRUE) 
ggsave("../figure_new/Projecting_spatial_slice44_each.pdf", width = 18, height = 12)




projecting = integrated_Exc@meta.data[rownames(cells_embed_predict), c("subcluster", "predict_projecting")]
projecting$predict_projecting = as.character(projecting$predict_projecting)
projecting = projecting[projecting$predict_projecting != "Rbp4", ]

projecting = round(prop.table(table(projecting), margin = 2),digits = 2)
projecting = as.data.frame(projecting)

#projecting = as.data.frame(table(projecting))

projecting$subcluster = factor(projecting$subcluster, levels = levels(Idents(merFISH_integrated))[1:18])
projecting$predict_projecting = stringr::str_replace(projecting$predict_projecting, "Cav_", "")

library(ggalluvial)
ggplot(data = projecting,
       aes(axis1 = subcluster, axis2 = predict_projecting, y = Freq)) +
  geom_alluvium(aes(fill = predict_projecting), width = 1/4) +
  scale_fill_manual(values = c(project_col[-1]))+
  geom_stratum(fill = c(L2_color[1:18], rev(project_col[-1])), width = 1/4) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() 
ggsave("../figure_new//Ext_projecting_regions_alluviumPlot_new.pdf", width = 8, height = 8)

# without normalization
projecting = integrated_Exc@meta.data[rownames(cells_embed_predict), c("subcluster", "predict_projecting")]
projecting$predict_projecting = as.character(projecting$predict_projecting)
projecting = projecting[projecting$predict_projecting != "Rbp4", ]

projecting = table(projecting)
#projecting_total = colSums(projecting)
#projecting = apply(projecting, 1, function(x) x/projecting_total * as.numeric(table(Fig3Rbp4retro$Cre.line)[-6]))
#projecting = as.data.frame(as.table(t(projecting)))

projecting = as.data.frame((projecting))
projecting$subcluster = factor(projecting$subcluster, levels = levels(Idents(merFISH_integrated))[1:18])
projecting$predict_projecting = stringr::str_replace(projecting$predict_projecting, "Cav_", "")

ggplot(data = projecting,
       aes(axis1 = subcluster, axis2 = predict_projecting, y = Freq)) +
  geom_alluvium(aes(fill = predict_projecting), width = 1/4) +
  scale_fill_manual(values = c(project_col[-1]))+
  geom_stratum(fill = c(L2_color[1:18], rev(project_col[-1])), width = 1/4) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() 


projecting = integrated_Exc@meta.data[rownames(cells_embed_predict), c("L3_cluster", "predict_projecting")]
projecting$predict_projecting = as.character(projecting$predict_projecting)
projecting = projecting[projecting$predict_projecting != "Rbp4", ]

projecting = round(prop.table(table(projecting), margin = 2),digits = 2)
projecting = as.data.frame(projecting)

#projecting = as.data.frame(table(projecting))

projecting$predict_projecting = stringr::str_replace(projecting$predict_projecting, "Cav_", "")

ggplot(data = projecting,
       aes(axis1 = L3_cluster, axis2 = predict_projecting, y = Freq)) +
  geom_alluvium(aes(fill = predict_projecting), width = 1/4) +
  scale_fill_manual(values = c(project_col[-1]))+
  geom_stratum(fill = c(L2_color[seq(1,17,2)[1:7]], rev(project_col[-1])), width = 1/4) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() 
ggsave("../figure_new//Ext_mainType_projecting_regions_alluviumPlot_new.pdf", width = 8, height = 8)


# ABA

integrated_Exc_plot1 = integrated_Exc_plot[integrated_Exc_plot$tech != "scRNA" & integrated_Exc_plot$ABA_PFC == "In", ]

integrated_Exc_plot1$ABA_metaRegion = factor(integrated_Exc_plot1$ABA_metaRegion, levels = c("ACAd", "PL", "ILA", "DP"))
p1 = ggplot(integrated_Exc_plot1, aes(predict_projecting, fill = ABA_metaRegion)) +
  geom_bar(position="fill") +
  ggsci::scale_fill_npg() +
  theme_cowplot() +
  theme(axis.text.x=element_text(angle=30, hjust=1))
p2 = ggplot(integrated_Exc_plot1, aes(predict_projecting, fill = L3_cluster)) +
  geom_bar(position="fill") +
  ggsci::scale_fill_ucscgb() +
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=30, hjust=1))
integrated_Exc_plot1$region_cluster = paste(integrated_Exc_plot1$ABA_metaRegion, integrated_Exc_plot1$L3_cluster)
p3 = ggplot(integrated_Exc_plot1, aes(predict_projecting, fill = region_cluster)) +
  geom_bar(position="fill", size=2) +
  ggsci::scale_fill_ucscgb() +
   ggsci::scale_color_npg() +
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=30, hjust=1))

p1/p2
ggsave("../figure_new//Ext_projecting_regions_barPlot_new.pdf", width = 8, height = 10)

ggplot(integrated_Exc_plot1, aes(x=1, fill = ABA_metaRegion)) +
  geom_bar(position="fill", width = 2) + 
  ggsci::scale_fill_ucscgb() +
  theme_void()+
  theme(aspect.ratio = 1) +
coord_polar(theta="y") + facet_wrap(.~predict_projecting)


integrated_Exc_plot1$layer = stringr::str_replace(integrated_Exc_plot1$L3_cluster, " .*", "") 
ggplot(integrated_Exc_plot1, aes(predict_projecting, fill = layer)) +
  geom_bar(position="fill", size=2) +
  coord_polar("y", start=0)+
  ggsci::scale_fill_ucscgb() +
  ggsci::scale_color_npg() +
  theme_cowplot()+
  theme(axis.text.x=element_text(angle=30, hjust=1)) 
a = table(integrated_Exc_plot1$layer, integrated_Exc_plot1$predict_projecting) 
apply(a, 1, function(x) x/colSums(a))


#
plot_target_to = function(ct){
  projecting = integrated_Exc@meta.data[rownames(cells_embed_predict), c("subcluster", "predict_projecting")]
  projecting = projecting[projecting$predict_projecting != "Rbp4", ]
  projecting = round(prop.table(table(projecting), margin = 2),digits = 2)
  projecting = as.data.frame((projecting))
  projecting = projecting[!is.nan(projecting$Freq), ]
  projecting = projecting[projecting$Freq != 0, ]
  projecting = projecting[projecting$subcluster == ct, ]
  projecting$fraction = projecting$Freq / sum(projecting$Freq)
  
  projecting$ymax <- cumsum(projecting$fraction)
  projecting$ymin <- c(0, head(projecting$ymax, n=-1))
  
  # Compute label position
  projecting$labelPosition <- (projecting$ymax + projecting$ymin) / 2
  
  # Compute a good label
  projecting$label <-  sprintf("%.1f%%", projecting$fraction * 100)
  
  
  ggplot(projecting, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=predict_projecting)) +
    geom_rect() +
    geom_text( x=3.5, aes(y=labelPosition, label=label), size=5) + 
    scale_fill_manual(values = c(as.character(pals::kelly()[3:8]))) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void()
}

p1 = plot_target_to("L6 CT 1")
p2 = plot_target_to("L6 CT 2")
p3 = plot_target_to("L6 CT 3")
p1|p2|p3
ggsave("../figure/Cell_projectint_targets_L6_CT.pdf", width = 14, height = 4)

plot_target_to("L5 ET 1")

SpatialDimPlot(integrated_Exc_sample,
               cells.highlight = CellsByIdentities(object = integrated_Exc_sample, idents = c("L6 CT 1", "L6 CT 2", "L6 CT 3")), 
               images = "slice.8",
               crop = T,
               stroke=.1,
               pt.size.factor = 2.3,
               ncol = 3,
               cols.highlight = c("#DE2D26", "grey80"),
               facet.highlight = TRUE) &
  theme(aspect.ratio = 1) &
  coord_flip() & scale_x_reverse() 


#

