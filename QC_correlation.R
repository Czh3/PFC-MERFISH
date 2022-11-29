#
# correlation with bulk RNAseq
bulk_expr = read.table("/nfs4/chaozhang/proj/Neuron/Aritra/bulk/stringtie/merge.TPM.txt", row.names = 1, header = T)

merFISH_integrated = readRDS(file = "./RDS/merFISH_final.0827.RDS")

merFISH_avg_expr = AverageExpression(merFISH_integrated, group.by = "orig.ident", slot = "counts")
merFISH_avg_expr = merFISH_avg_expr$RNA
merFISH_avg_expr = as.data.frame(merFISH_avg_expr) 
ggplot(as.data.frame(merFISH_avg_expr), aes(log(mouse3.4 + 0.01), log(Pain1+ 0.01))) +
  geom_point(size=1, color="blue4") + 
  geom_abline(slope = 1, color = "gray60")+ 
  xlab("MERFISH rep1 (counts/cell)") + ylab("MERFISH rep2 (counts/cell)") +
  cowplot::theme_cowplot() + ggtitle("Spearman cor: 0.960")
ggsave("../figure_new/merFISH_rep_cor.scatter.pdf", width = 4, height = 4)
cor(log(merFISH_avg_expr$mouse3.4 + 0.01), log(merFISH_avg_expr$Pain1 + 0.01), method = "spearman")

merFISH_avg_expr.cor = cor(merFISH_avg_expr[,c(2:20,1)])
pheatmap::pheatmap(merFISH_avg_expr.cor, cluster_rows = F, cluster_cols = F,
                   display_numbers = T, number_format = "%.2f", number_color = "black",
                   border_color = NA,
                   color = colorRampPalette(c("blue", "white", "red2"))(50))

merFISH_avg_expr = as.data.frame(rowMeans(merFISH_avg_expr))
genes = intersect(rownames(merFISH_avg_expr), rownames(bulk_expr))

merFISH_bulk = cbind(merFISH_avg_expr[genes, ], bulk_expr[genes, ])
colnames(merFISH_bulk) = c("merFISH", "bulk_rep1", "bulk_rep2")
merFISH_bulk$bulk = rowMeans(merFISH_bulk[, c(2:3)])

ggplot(merFISH_bulk, aes(log(bulk + 1), log(merFISH+ 0.01))) +
  geom_point(size=1, color="blue4") + 
  cowplot::theme_cowplot() + ggtitle("Spearman cor: 0.81")
ggsave("../figure_new//merFISH_bulk_cor.scatter.pdf", width = 4, height = 4)
cor(log(merFISH_bulk$merFISH + 0.01), log(merFISH_bulk$bulk + 1), method = "spearman")


