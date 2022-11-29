# compare with MOp


MOp = readRDS("../../merFISH_MOp/MOp_mouse1sample3.RDS")

SpatialPlot(MOp, group.by = "subclass" )

genes = intersect(rownames(MOp), rownames(merFISH_integrated))

# All
merFISH_integrated_L2_1 = merFISH_integrated
merFISH_integrated_L2_1 = NormalizeData(merFISH_integrated_L2_1[genes,])
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "L3_cluster", slot = "data")

MOp_1 = NormalizeData(MOp[genes,])
merFISH_MOp_expr = AverageExpression(MOp_1, group.by = "subclass", slot = "data")

PFC_avg_expr.scale = t(scale(t(merFISH_avg_expr$RNA[genes,])))
MOp_avg_expr.scale = t(scale(t(merFISH_MOp_expr$RNA[genes, ])))

Expr_cor = cor(PFC_avg_expr.scale, MOp_avg_expr.scale)
pdf("../figure_new/compareMOp.cluster.pdf", 7, 6)
pheatmap::pheatmap(Expr_cor[,c(3,4,6,9,10,5,7,8,11,12,18,21,22,23,1,2,13,14,15,24)], cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("gray90", "white", "red2"))(50), border_color = NA)
dev.off()



# Ext
merFISH_integrated_L2_1 = subset(merFISH_integrated, CellType == "Ext")
merFISH_integrated_L2_1 = NormalizeData(merFISH_integrated_L2_1[genes,])
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "subcluster", slot = "data")

MOp_1 = subset(MOp, class_label == "Glutamatergic")
MOp_1 = NormalizeData(MOp_1[genes,])
merFISH_MOp_expr = AverageExpression(MOp_1, group.by = "label", slot = "data")

PFC_avg_expr.scale = t(scale(t(merFISH_avg_expr$RNA[genes,])))
MOp_avg_expr.scale = t(scale(t(merFISH_MOp_expr$RNA[genes, ])))

Expr_cor = cor(PFC_avg_expr.scale, MOp_avg_expr.scale)
pdf("../figure_new/compareMOp.Exc_cluster.pdf", 10, 8)
pheatmap::pheatmap(Expr_cor[,c(1:7,11,12,8:10, 16,17,15,13,14,19,21,18,20,22,23,
                               24:32,37:39,33:36)], cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("gray90", "white", "red2"))(50), border_color = NA)
dev.off()


# Inh
merFISH_integrated_L2_1 = subset(merFISH_integrated, CellType == "Inh")
merFISH_integrated_L2_1 = NormalizeData(merFISH_integrated_L2_1[genes,])
merFISH_avg_expr = AverageExpression(merFISH_integrated_L2_1, group.by = "subcluster", slot = "data")

MOp_1 = subset(MOp, class_label == "GABAergic")
MOp_1 = NormalizeData(MOp_1[genes,])
merFISH_MOp_expr = AverageExpression(MOp_1, group.by = "label", slot = "data")

PFC_avg_expr.scale = t(scale(t(merFISH_avg_expr$RNA[genes,])))
MOp_avg_expr.scale = t(scale(t(merFISH_MOp_expr$RNA[genes, ])))

Expr_cor = cor(PFC_avg_expr.scale, MOp_avg_expr.scale)
pdf("../figure_new/compareMOp.Inh_cluster.pdf", 7, 8)
pheatmap::pheatmap(Expr_cor, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("gray90", "white", "red2"))(50), border_color = NA)
dev.off()

summary(MOp$nCount_RNA)
summary(merFISH_integrated$nCount_RNA)
boxplot(log10(MOp$nCount_RNA), log10(merFISH_integrated$nCount_RNA))

boxplot(log10(MOp_1$nCount_RNA), log10(merFISH_integrated_L2_1$nCount_RNA))


