library(Seurat)
library(SeuratData)
library(dplyr)


# stable cluster
one_bootstrap = function(obj, res, origin_res_label, algorithm, k.param = 20, reduction=reduction, dims=dims){

  obj_random = obj[, sample(colnames(obj), length(colnames(obj)) * 0.5)]
  obj_random = FindNeighbors(obj_random, k.param = k.param,  reduction=reduction, dims=dims)
  obj_random = FindClusters(obj_random, algorithm = algorithm, resolution = res,
                            method = "igraph",
                            group.singletons = F, verbose = FALSE)
  
  returns = as.data.frame(Idents(obj_random))
  returns$origin = obj[[origin_res_label]][colnames(obj_random),1]
  returns = data.frame(table(returns))
  cell_number = table(Idents(obj_random))
  #cell_number = table(obj[[origin_res_label]][colnames(obj_random),1])
  returns$cell_n = cell_number[as.character(returns$Idents.obj_random.)]
  #returns$cell_n = cell_number[as.character(returns$origin)]
  returns$percent = returns$Freq / returns$cell_n
  
  returns1 <- returns %>%
    group_by(origin) %>%
    dplyr::summarise(percent = max(percent))
  return(data.frame(returns1))
}

stable_cluster = function(obj, res, origin_res_label, algorithm, k.param = 20, reduction="pca", dims = 1:30){
  returns2 = lapply(1:10, function(x) one_bootstrap(obj, res, origin_res_label, algorithm, k.param = k.param, reduction=reduction, dims=dims) )

  returns3 = do.call(cbind.data.frame, returns2)
  returns3 = returns3[,c(1,seq(2,20,2))]
  returns3$mean = rowMeans(returns3[,2:11])
  returns3 = returns3[,c(1,12)]
  return(returns3)
}

# determine K
determineK = function(obj, res, algorithm, reduction=reduction, dims=dims){
  k = c(10, 15, 20, 30, 40, 50)
  cat(k)
  returns2 = lapply(k, function(x) {
    obj <- FindNeighbors(obj, 
                        k.param = x,
                        reduction=reduction, dims=dims)
    
    obj <- FindClusters(obj, 
                        resolution = res,
                        method = "igraph",
                        algorithm = algorithm,
                        group.singletons = F,
                        verbose = FALSE)
    origin_res_label =  paste0("RNA_snn_res.", res)
    
    return1 = one_bootstrap(obj, res, origin_res_label, algorithm, k.param = x, reduction=reduction, dims=dims)
    return1 = return1[return1$percent >= 0.4, 1]
    return(c(sum(table(obj[[origin_res_label]][,1]) > 10),
             sum(obj[[origin_res_label]][,1] %in% return1) ))
  })
  returns3 = do.call(cbind.data.frame, returns2)
  rownames(returns3) = c("cluters", "cells")
  colnames(returns3) = k
  returns3 = as.data.frame(t(returns3))
  returns3$k = factor(rownames(returns3), levels = rownames(returns3))
  g = ggplot(returns3, aes(cells, cluters, color = k)) +
    geom_point(size=3)+
    theme_classic() +
    scale_color_manual(values = as.vector(pals::alphabet()))
  return(g)
}




