# cell interaction
merFISH_integrated_L3 = readRDS(file = "./RDS/merFISH_integrated_L3_final_0719.seurat.rds")

merFISH_integrated = readRDS(file = "./RDS/merFISH_final.0827.RDS")


cal_CCI = function(obj, s){
  
  merFISH_integrated_test = suppressWarnings(subset(obj, slice == s))
  
  embeddings = data.matrix(merFISH_integrated_test@images[[paste0('slice.', s)]]@coordinates[,1:2])
  colnames(embeddings) = c("Loc.1", "loc.2")
  
  merFISH_location =  CreateDimReducObject(
    embeddings = embeddings,
    key = "LOC_",
    assay = "RNA",
    global = TRUE
  )
  merFISH_integrated_test[['LOC']] = merFISH_location

  merFISH_integrated_test = FindNeighbors(merFISH_integrated_test, reduction = "LOC", dims = 1:2,k.param = 31,
                                          annoy.metric = "euclidean", return.neighbor = T)
  
  # cell cell interaction
  CCI = lapply(levels(Idents(merFISH_integrated)), function(x){
    # obs.
    cell_query = colnames(merFISH_integrated_test)[merFISH_integrated_test$subcluster == x] 
    NN = Seurat::TopNeighbors(merFISH_integrated_test@neighbors$RNA.nn, cell_query, n = 31)[-1]
    # random
    cell_random = sample(colnames(merFISH_integrated_test), 30 * length(cell_query), replace = T)
    if(length(NN) == 0){
      res = rep(x, length(levels(Idents(merFISH_integrated))))
      res = cbind(res, 0, 0)
      rownames(res) = levels(Idents(merFISH_integrated))
      res
    }else{
      cbind(x, table(merFISH_integrated_test$subcluster[NN])+1, table(merFISH_integrated_test$subcluster[cell_random])+1)
    }
  })
  CCI = do.call(rbind, CCI)
  CCI = data.frame(CCI)
  CCI$V2 = as.numeric(CCI$V2)
  CCI$V3 = as.numeric(CCI$V3)
  CCI
}

CCI_allSlice = lapply(unique(merFISH_integrated_L3$slice), cal_CCI, obj = merFISH_integrated_L3)

CCI_allSlice = do.call(cbind, CCI_allSlice)


CCI_allSlice$Obs = apply(CCI_allSlice, 1, function(x) mean(as.numeric(x[seq(2, length(x), 3)])))
CCI_allSlice$Exp = apply(CCI_allSlice, 1, function(x) mean(as.numeric(x[seq(3, length(x), 3)])))
CCI_allSlice$log2ObsExp = log2(CCI_allSlice$Obs/CCI_allSlice$Exp)
CCI_allSlice$log2ObsExp[CCI_allSlice$log2ObsExp > 2] <- 2
CCI_allSlice$log2ObsExp[CCI_allSlice$log2ObsExp < -2] <- -2
CCI_allSlice$pval = apply(CCI_allSlice, 1, function(x){
  wilcox.test(as.numeric(x[seq(2, length(x), 3)]), as.numeric(x[seq(3, length(x), 3)]))$p.value
})
CCI_allSlice$FDR = p.adjust(CCI_allSlice$pval, method = "BH")
CCI_allSlice$y = levels(merFISH_integrated_L3$subcluster)

CCI_allSlice = CCI_allSlice[, c(1, (ncol(CCI_allSlice)-5):ncol(CCI_allSlice))]
CCI_allSlice$x = factor(CCI_allSlice$x, levels = (levels(merFISH_integrated_L3$subcluster)[-19]))
CCI_allSlice$y = factor(CCI_allSlice$y, levels = rev(levels(merFISH_integrated_L3$subcluster)[-19]))

CCI_allSlice$FDR[CCI_allSlice$FDR < 1e-6 ] <- 1e-6
ggplot(na.omit(CCI_allSlice), aes(x, y, color = log2ObsExp, size = -log10(FDR)))+
  geom_point() +
  #scale_color_gradientn(colours = pals::coolwarm(50)) +
  scale_color_gradientn(colours = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(50))+
  theme_bw(base_size = 18) + xlab("") + ylab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) 
ggsave("../figure_new//CellCellInteraction.heatmap.pdf", width = 13, height = 11)


CCI_allSlice$Obs = apply(CCI_allSlice, 1, function(x) mean(as.numeric(x[seq(2, length(x), 3)])))
CCI_allSlice$Exp = apply(CCI_allSlice, 1, function(x) mean(as.numeric(x[seq(3, length(x), 3)])))
CCI_allSlice$log2ObsExp = log2(CCI_allSlice$Obs/CCI_allSlice$Exp)

CCI_allSlice.df = reshape2::dcast(na.omit(CCI_allSlice), x ~ y, value.var = "log2ObsExp")
rownames(CCI_allSlice.df) = CCI_allSlice.df$x
CCI_allSlice.df = CCI_allSlice.df[,-1]
CCI_allSlice.df = as.matrix(CCI_allSlice.df)
CCI_allSlice.df = CCI_allSlice.df[, rev(colnames(CCI_allSlice.df))]
diag(CCI_allSlice.df) <- 0
pheatmap::pheatmap(CCI_allSlice.df, 
                   #cluster_rows = F, cluster_cols = F,
                   border_color = NA,
                   color = colorRampPalette(c("blue3", "blue", "white", "red2"))(100),
                   clustering_method = "ward.D"
                   )

# In/Out PFC
merFISH_integrated_PFC = subset(merFISH_integrated, ABA_PFC == "In")
merFISH_integrated_PFC$subcluster = factor(merFISH_integrated_PFC$subcluster, levels = levels(Idents(merFISH_integrated_PFC)))

CCI_allSlice_PFC = lapply(unique(merFISH_integrated_PFC$slice), cal_CCI, obj = merFISH_integrated_PFC)

CCI_allSlice_PFC = do.call(cbind, CCI_allSlice_PFC)


CCI_allSlice_PFC$Obs = apply(CCI_allSlice_PFC, 1, function(x) mean(as.numeric(x[seq(2, length(x), 3)])))
CCI_allSlice_PFC$Exp = apply(CCI_allSlice_PFC, 1, function(x) mean(as.numeric(x[seq(3, length(x), 3)])))
CCI_allSlice_PFC$log2ObsExp = log2(CCI_allSlice_PFC$Obs/CCI_allSlice_PFC$Exp)
CCI_allSlice_PFC$log2ObsExp[CCI_allSlice_PFC$log2ObsExp > 2] <- 2
CCI_allSlice_PFC$log2ObsExp[CCI_allSlice_PFC$log2ObsExp < -2] <- -2
CCI_allSlice_PFC$pval = apply(CCI_allSlice_PFC, 1, function(x){
  wilcox.test(as.numeric(x[seq(2, length(x), 3)]), as.numeric(x[seq(3, length(x), 3)]))$p.value
})
CCI_allSlice_PFC$FDR = p.adjust(CCI_allSlice_PFC$pval, method = "BH")
CCI_allSlice_PFC$y = levels(Idents(merFISH_integrated))

CCI_allSlice_PFC = CCI_allSlice_PFC[, c(1, (ncol(CCI_allSlice_PFC)-5):ncol(CCI_allSlice_PFC))]
CCI_allSlice_PFC$x = factor(CCI_allSlice_PFC$x, levels = levels(Idents(merFISH_integrated)))
CCI_allSlice_PFC$y = factor(CCI_allSlice_PFC$y, levels = rev(levels(Idents(merFISH_integrated))))

CCI_allSlice_PFC$FDR[CCI_allSlice_PFC$FDR < 1e-4 ] <- 1e-4
ggplot(na.omit(CCI_allSlice_PFC), aes(x, y, color = log2ObsExp, size = -log10(FDR)))+
  geom_point() +
  #scale_color_gradientn(colours = pals::coolwarm(50)) +
  scale_color_gradientn(colours = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(50))+
  theme_bw(base_size = 18) + xlab("") + ylab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) 
ggsave("../figure_new//CellCellInteraction_PFC.heatmap.pdf", width = 13, height = 11)


CCI_allSlice_PFC.cluster = dcast(CCI_allSlice_PFC, x ~ y, value.var = "log2ObsExp")
rownames(CCI_allSlice_PFC.cluster) = CCI_allSlice_PFC.cluster$x
CCI_allSlice_PFC.cluster = CCI_allSlice_PFC.cluster[,-1]
CCI_allSlice_PFC.cluster = CCI_allSlice_PFC.cluster[levels(Idents(merFISH_integrated_L3))[1:37],
                                                    levels(Idents(merFISH_integrated_L3))[1:37]]
diag(CCI_allSlice_PFC.cluster) <- NA
CCI_allSlice_PFC.cluster[CCI_allSlice_PFC.cluster > 1.5] <- 1.5
CCI_allSlice_PFC.cluster[CCI_allSlice_PFC.cluster < -1.5] <- -1.5
pheatmap::pheatmap(CCI_allSlice_PFC.cluster, #clustering_method = "ward.D2",
                   )


# cell complexity
cal_CC = function(s){
  
  merFISH_integrated_test = suppressWarnings(subset(merFISH_integrated_L3, slice == s))
  
  embeddings = data.matrix(merFISH_integrated_test@images[[paste0('slice.', s)]]@coordinates[,1:2])
  colnames(embeddings) = c("Loc.1", "loc.2")
  
  merFISH_location =  CreateDimReducObject(
    embeddings = embeddings,
    key = "LOC_",
    assay = "RNA",
    global = TRUE
  )
  merFISH_integrated_test[['LOC']] = merFISH_location

  merFISH_integrated_test = FindNeighbors(merFISH_integrated_test, reduction = "LOC", dims = 1:2, k.param = 31,
                                          annoy.metric = "euclidean", return.neighbor = T)
  
  # Neighbor complexity 
  CC = parallel::mclapply(colnames(merFISH_integrated_test), function(x){
    NN = Seurat::TopNeighbors(merFISH_integrated_test@neighbors$RNA.nn, x, n = 31)[-1]
    c(x, sum(table(merFISH_integrated_test$subcluster[NN]) != 0 ))
  }, mc.cores = 8)
  
  CC = do.call(rbind, CC)
  CC = as.data.frame(CC)
  CC$V2 = as.numeric(CC$V2)
  CC$Cell = merFISH_integrated_test$subcluster[CC$V1]
  CC
}


CC_allSlice = lapply(unique(merFISH_integrated_L3$slice), cal_CC)

CC_allSlice = do.call(rbind, CC_allSlice)
CC_allSlice$CT = merFISH_integrated_L3$CellType[CC_allSlice$V1]

ggplot(CC_allSlice, aes(x = Cell, y = V2, fill = Cell))+
  geom_violin(bw=0.5)+
  stat_summary(fun=mean, geom="point", size=2, color="red")+
  cowplot::theme_cowplot() + ylab("Neighbor complexity")+
  scale_fill_manual(values = L2_color) +
  theme(legend.position = "None") +
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1)) 
ggsave("../figure_new//CellNeighborComplexity.violinplot.pdf", width = 12, height = 3)


ggplot(CC_allSlice, aes(x = CT, y = V2, fill = CT))+
  geom_violin(bw=0.5)+
  stat_summary(fun=mean, geom="point", size=2, color="red")+
  cowplot::theme_cowplot() + ylab("Neighbor complexity")+
  scale_fill_manual(values = L2_color) +
  theme(legend.position = "None") +
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 30, hjust=1)) 



library(ggridges)
ggplot(CC_allSlice[CC_allSlice$CT == "Ext", ], aes(y = Cell, x=V2, fill = Cell))+
  #stat_bin(aes(y=..density..),geom="step", bins = 20)+
  geom_density_ridges(bandwidth = .5, quantile_lines = TRUE, quantile_fun=function(x,...)mean(x))+
  cowplot::theme_cowplot() + xlab("Neighbor complexity") + ylab("")+
  scale_fill_manual(values = L2_color)
ggsave("../figure_new//CellNeighborComplexity_Ext.violinplot.pdf", width = 5, height = 5)

CC_allSlice$project = stringr::str_replace(CC_allSlice$Cell, " [0-9]", "")
CC_allSlice$project = stringr::str_replace(CC_allSlice$project, ".* ", "")

ggplot(CC_allSlice[CC_allSlice$project == "IT", ], aes(y = Cell, x=V2, fill = Cell))+
  #stat_bin(aes(y=..density..),geom="step", bins = 20)+
  geom_density_ridges(bandwidth = .5, quantile_lines = TRUE, quantile_fun=function(x,...)mean(x))+
  cowplot::theme_cowplot() + xlab("Neighbor complexity") + ylab("")+
  scale_fill_manual(values = L2_color)
ggsave("../figure_new//CellNeighborComplexity_IT.violinplot.pdf", width = 5, height = 4)






# moleclar distance

harmony_embed = merFISH_integrated_L3@reductions$harmony@cell.embeddings
harmony_embed = as.data.frame(harmony_embed)
harmony_embed = aggregate(x= harmony_embed,     
                          by = list(merFISH_integrated_L3$subcluster),      
                          FUN = mean)
rownames(harmony_embed) = harmony_embed$Group.1
harmony_embed = harmony_embed[, -1]
harmony_embed.cor = cor(t(harmony_embed))
harmony_embed.cor = reshape2::melt(harmony_embed.cor)
harmony_embed.cor$value[harmony_embed.cor$value > 0.8] <- 0.8
ggplot(harmony_embed.cor, aes(Var1, Var2, color = value, size = value))+
  geom_point() +
  #scale_color_gradientn(colours = pals::coolwarm(50)) +
  scale_color_gradientn(colours = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(50))+
  theme_bw(base_size = 18) + xlab("") + ylab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) 




merFISH_avg_expr = AverageExpression(merFISH_integrated_L3, group.by = "subcluster", slot = "scale.data")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 

merFISH_avg_expr.cor = cor(merFISH_avg_expr, method = "spearman")
pheatmap::pheatmap(merFISH_avg_expr.cor)

merFISH_avg_expr.cor = reshape2::melt(merFISH_avg_expr.cor)
merFISH_avg_expr.cor$Var2 =  factor(merFISH_avg_expr.cor$Var2, levels = rev(levels(merFISH_integrated_L3$subcluster)[-19]))
ggplot(merFISH_avg_expr.cor, aes(Var1, Var2, color = value, size = abs(value)))+
  geom_point() +
  #scale_color_gradientn(colours = pals::coolwarm(50)) +
  scale_color_gradientn(colours = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(50)[10:50])+
  theme_bw(base_size = 18) + xlab("") + ylab("")+
  theme(axis.text=element_text(colour="black"),
        axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) 
ggsave("../figure_new//CellType_exprCorrelation.heatmap.pdf", width = 12, height = 10)




