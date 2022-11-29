
setwd("/nfs4/chaozhang/proj/Neuron/Aritra/merFISH/script/")

#devtools::install_github("tractatus/wholebrain") 

library(wholebrain)

match_CCF = function(slice, 
                     slice.col = 1,
                     ABA_coronal = 50,
                     ABA_color = "left",
                     x.scale = 8.5, y.scale = 8.5,
                     x.move = 0, y.move = 0,
                     rotation = 0,
                     flip.x = FALSE,
                     flip.y = FALSE,
                     pt.size = 0.3,
                     plot.axis = TRUE){
  
  ## plot CCF
  
  k = ABA_coronal
  numPaths <- EPSatlas$plates[[k]]@summary@numPaths


  style = as.character(EPSatlas$plate.info[[k]]$style)
  
  xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x)
  xmax <- max(EPSatlas$plates[[k]][[1]]@paths$path@x)
  
  plot(EPSatlas$plates[[k]][[1]]@paths$path@x, 
       EPSatlas$plates[[k]][[1]]@paths$path@y, col = 0, 
       #xlim = c(-xmax, xmax), ylim = c(0, 68234.56),
       xlim = c(-50000, 50000), ylim = c(0, 60000),
       axes = plot.axis, ylab = "", xlab = "", asp = 1, main = "")
  
  
  
  plot.data = merFISH_integrated@meta.data[merFISH_integrated$slice == slice,  ]
  
  if(slice.col == 1){
    L3_color = c( pals::parula(7), rep("gray", 6),
                  "#AA0DFE","#993F00","#F99379","#604E97", "maroon2", "yellow3")
  }else if(slice.col == 2){
    L3_color = as.character(pals::polychrome()[3:30])
  }else{
    L3_color = as.character(pals::glasbey()[1:30])
  }

  L3_color = L3_color[plot.data$L3_cluster]
  names(L3_color) = rownames(plot.data)
  L3_color = scales::alpha(L3_color, 0.8)
  
  
  x0 = (plot.data$centroid_1 - median(plot.data$centroid_1)) 
  y0 = (plot.data$centroid_2 - median(plot.data$centroid_2))
  
  
  a = rotation / 180 * pi
  
  x1 = x0 * cos(a) - y0 * sin(a)
  y1 = x0 * sin(a) + y0 * cos(a)
  
  x1 = x1 * x.scale + x.move
  y1 = y1 * y.scale + y.move
  
  if(flip.x == TRUE){
    x1 = mean(x1) - x1 + mean(x1) 
  }
  if(flip.y == TRUE){
    y1 = mean(y1) - y1 + mean(y1) 
  }
  
  points(x1,
         y1,
         pch = 16, col = L3_color, cex = pt.size, add = T)
  
  
  lapply(1:numPaths, function(N) {
    if(ABA_color == "right"){
      polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - xmin,
              EPSatlas$plates[[k]][[N]]@paths$path@y, 
              col = style[N], border = "black", lwd = 1, 
              lty = 3)
    }else{
      polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - xmin,
              EPSatlas$plates[[k]][[N]]@paths$path@y, 
              border = "black", lwd = 1, 
              lty = 3)
    }
    
    if(ABA_color == "left"){
      polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                  xmin ), EPSatlas$plates[[k]][[N]]@paths$path@y, 
              col = style[N], border = "black", lwd = 1, 
              lty = 3)
    }else{
      polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                  xmin ), EPSatlas$plates[[k]][[N]]@paths$path@y, 
              border = "black", lwd = 1, 
              lty = 3)
    }
    
    
  })
  
  
  # align to ABA regions
  
  namelist<-as.numeric(as.character(EPSatlas$plate.info[[k]]$structure_id))
  colorlist<-as.character(EPSatlas$plate.info[[k]]$style)
  neuronregion<-rep(0,length(x1))
  neuroncolor<-rep('#000000',length(x1))
  right.hemisphere<-rep(NA, length(x1))
  
  
  for(i in 1:numPaths) {
    temp = point.in.polygon(x1, y1, EPSatlas$plates[[k]][[i]]@paths$path@x - xmin, EPSatlas$plates[[k]][[i]]@paths$path@y)
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-FALSE
    
    temp = point.in.polygon(x1, y1, -(EPSatlas$plates[[k]][[i]]@paths$path@x - xmin), EPSatlas$plates[[k]][[i]]@paths$path@y)
    neuronregion[which(temp==1)]<-namelist[i]
    neuroncolor[which(temp==1)]<-as.character(colorlist[i])
    right.hemisphere[which(temp==1)]<-TRUE
    
  }
  
  plot.data$reloc_x = x1
  plot.data$reloc_y = y1
  plot.data$ABA_id<- neuronregion
  plot.data$ABA_color<- neuroncolor
  plot.data$ABA_hemisphere<- ifelse(right.hemisphere, "Right", "Left")
  plot.data$ABA_regions = as.character(acronym.from.id(neuronregion))
  plot.data$ABA_metaRegion = stringr::str_replace(plot.data$ABA_regions, "[0-9].*", "")
  plot.data$ABA_PFC = ifelse(plot.data$ABA_metaRegion %in% c("ACAd", "ACAv", "PL", "DP", "ILA", "ORBm"), "In", "Out")
  return(plot.data)
}

data(EPSatlas, envir=environment())
data(atlasIndex, envir=environment())
merFISH_integrated_L3 = readRDS("./RDS/merFISH_integrated_L3_final_0719.seurat.rds")
merFISH_integrated = readRDS("./RDS/merFISH_final.0826.RDS")

ABA_region_plot = function(data){
  ggplot(data, aes(reloc_x, reloc_y, color = ABA_metaRegion))+
    geom_point(pch = 16, size = 1) +
    theme_void() + coord_fixed()+
    scale_color_manual(values = as.character(pals::alphabet()))
}

slice33 = match_CCF(slice = 33,
                    ABA_coronal = 31,
                    #ABA_style = "color",
                    rotation = 182,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 3000,
                    y.move = 42000,
                    pt.size = 0.3)
ABA_region_plot(slice33)

#
mySpatialPlot_example(merFISH_integrated_L3, c(as.character(pals::polychrome()[3:34]), as.character(pals::kelly())), margin_x=2000, margin_y=2000)

slice27 = match_CCF(slice = 27,
                    ABA_coronal = 40,
                    ABA_col = "left",
                    rotation = 93,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 7000,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)
ABA_region_plot(slice27)


slice31 = match_CCF(slice = 31,
                    ABA_coronal = 36,
                    ABA_col = "left",
                    rotation = 273,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 11000,
                    y.move = 41000,
                    pt.size = 0.3,
                    plot.axis = F)

slice30 = match_CCF(slice = 30,
                    ABA_coronal = 29,
                    ABA_col = "left",
                    rotation = 275,
                    x.scale = 8.3,
                    y.scale = 8.3,
                    x.move = 11000,
                    y.move = 42500,
                    pt.size = 0.3,
                    plot.axis = F)

slice28 = match_CCF(slice = 28,
                    ABA_coronal = 27,
                    rotation = 93,
                    x.scale = 9,
                    y.scale = 9.4,
                    x.move = 10000,
                    y.move = 41000,
                    pt.size = 0.3,
                    plot.axis = F)


slice29 = match_CCF(slice = 29,
                    ABA_coronal = 25,
                    rotation = -10,
                    x.scale = 10.5,
                    y.scale = 10.5,
                    x.move = 10000,
                    y.move = 44000,
                    pt.size = 0.3,
                    plot.axis = F)



#
slice23 = match_CCF(slice = 23,
                    ABA_coronal = 28,
                    rotation = 205,
                    x.scale = 7.8,
                    y.scale = 8.1,
                    x.move = 12000,
                    y.move = 42300,
                    flip.y = T,
                    pt.size = 0.3,
                    plot.axis = F)

slice24 = match_CCF(slice = 24,
                    ABA_coronal = 37,
                    rotation = 90,
                    x.scale = 8.0,
                    y.scale = 8.5,
                    x.move = 9000,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)

slice25 = match_CCF(slice = 25,
                    #slice.col = 2,
                    ABA_coronal = 39,
                    rotation = 85,
                    x.scale = 8.2,
                    y.scale = 8.2,
                    x.move = 6000,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)


slice26 = match_CCF(slice = 26,
                    slice.col = 1,
                    ABA_coronal = 25,
                    rotation = -90,
                    x.scale = 8.2,
                    y.scale = 8.2,
                    x.move = 10000,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)

#
slice8 = match_CCF(slice = 8,
                    slice.col = 1,
                    ABA_coronal = 31,
                    rotation = 185,
                    x.scale = 8.2,
                    y.scale = 8.2,
                    x.move = 11000,
                    y.move = 43000,
                    pt.size = 0.3,
                    plot.axis = F)

slice9 = match_CCF(slice = 9,
                   slice.col = 1,
                   ABA_coronal = 38,
                   rotation = 175,
                   flip.y = T,
                   x.scale = 8.8,
                   y.scale = 8.8,
                   x.move = 8000,
                   y.move = 44000,
                   pt.size = 0.3,
                   plot.axis = F)
#
slice10 = match_CCF(slice = 10,
                    slice.col = 2,
                    ABA_coronal = 31,
                    ABA_color = "no",
                    rotation = 210,
                    x.scale = 9,
                    y.scale = 9,
                    x.move = 4000,
                    y.move = 40000,
                    pt.size = 0.3,
                    plot.axis = F)

slice11 = match_CCF(slice = 11,
                   slice.col = 1,
                   ABA_coronal = 40,
                   ABA_color = "no",
                   rotation = 0,
                   x.scale = 10,
                   y.scale = 10,
                   x.move = -2000,
                   y.move = 43000,
                   pt.size = 0.3,
                   plot.axis = F)

slice12 = match_CCF(slice = 12,
                    slice.col = 1,
                    ABA_coronal = 39,
                    ABA_color = "no",
                    rotation = -3,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = -2000,
                    y.move = 43000,
                    pt.size = 0.3,
                    plot.axis = F)

#
slice13 = match_CCF(slice = 13,
                    slice.col = 1,
                    ABA_coronal = 40,
                    ABA_color = "no",
                    rotation = 155,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = 2000,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)
slice14 = match_CCF(slice = 14,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "no",
                    rotation = 155,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = 000,
                    y.move = 40000,
                    pt.size = 0.3,
                    plot.axis = F)
# pain
slice32 = match_CCF(slice = 32,
                    slice.col = 1,
                    ABA_coronal = 36,
                    ABA_color = "no",
                    rotation = 185,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = 3000,
                    y.move = 39000,
                    pt.size = 0.3,
                    plot.axis = F)
slice33 = match_CCF(slice = 33,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "no",
                    rotation = 180,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = 3000,
                    y.move = 40000,
                    pt.size = 0.3,
                    plot.axis = F)

slice34 = match_CCF(slice = 34,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "no",
                    rotation = 9,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = -3800,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)

slice35 = match_CCF(slice = 35,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "no",
                    rotation = 9,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = -3800,
                    y.move = 38000,
                    pt.size = 0.3,
                    plot.axis = F)

slice36 = match_CCF(slice = 36,
                    slice.col = 1,
                    ABA_coronal = 41,
                    ABA_color = "no",
                    rotation = 173,
                    x.scale = 9.8,
                    y.scale = 9.8,
                    x.move = 2700,
                    y.move = 46500,
                    pt.size = 0.3,
                    plot.axis = F)
slice37 = match_CCF(slice = 37,
                    slice.col = 1,
                    ABA_coronal = 44,
                    ABA_color = "no",
                    rotation = 177,
                    x.scale = 9.8,
                    y.scale = 9.6,
                    x.move = 2500,
                    y.move = 48500,
                    pt.size = 0.3,
                    plot.axis = F)

slice38 = match_CCF(slice = 38,
                    slice.col = 1,
                    ABA_coronal = 42,
                    ABA_color = "no",
                    rotation = 193,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 1800,
                    y.move = 48000,
                    pt.size = 0.3,
                    plot.axis = F)
slice39 = match_CCF(slice = 39,
                    slice.col = 1,
                    ABA_coronal = 47,
                    ABA_color = "no",
                    rotation = 13,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = -2100,
                    y.move = 50000,
                    pt.size = 0.3,
                    plot.axis = F)

slice40 = match_CCF(slice = 40,
                    slice.col = 1,
                    ABA_coronal = 29,
                    ABA_color = "no",
                    rotation = 2,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = -3100,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)
slice41 = match_CCF(slice = 41,
                    slice.col = 1,
                    ABA_coronal = 29,
                    ABA_color = "no",
                    rotation = 108,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 500,
                    y.move = 41000,
                    pt.size = 0.3,
                    plot.axis = F)

slice42 = match_CCF(slice = 42,
                    slice.col = 1,
                    ABA_coronal = 37,
                    ABA_color = "no",
                    rotation = -3,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = -2200,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)
slice43 = match_CCF(slice = 43,
                    slice.col = 1,
                    ABA_coronal = 39,
                    ABA_color = "no",
                    rotation = 94,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = -1200,
                    y.move = 43500,
                    pt.size = 0.3,
                    plot.axis = F)

slice44 = match_CCF(slice = 44,
                    slice.col = 1,
                    ABA_coronal = 37,
                    ABA_color = "no",
                    rotation = -77,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = -1400,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)
slice45 = match_CCF(slice = 45,
                    slice.col = 1,
                    ABA_coronal = 37,
                    ABA_color = "no",
                    rotation = 185,
                    x.scale = 8.5,
                    y.scale = 8.5,
                    x.move = 2000,
                    y.move = 43000,
                    pt.size = 0.3,
                    plot.axis = F)

#
slice46 = match_CCF(slice = 46,
                    slice.col = 1,
                    ABA_coronal = 32,
                    ABA_color = "no",
                    rotation = 175,
                    flip.x = T,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 10500,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)

slice47 = match_CCF(slice = 47,
                    slice.col = 1,
                    ABA_coronal = 32,
                    ABA_color = "no",
                    rotation = 190,
                    flip.x = T,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 11500,
                    y.move = 43000,
                    pt.size = 0.3,
                    plot.axis = F)
slice48 = match_CCF(slice = 48,
                    slice.col = 1,
                    ABA_coronal = 32,
                    ABA_color = "no",
                    rotation = 185,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 13000,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)
#

slice49 = match_CCF(slice = 49,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "left",
                    rotation = -12,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 8500,
                    y.move = 45000,
                    pt.size = 0.3,
                    plot.axis = F)
slice50 = match_CCF(slice = 50,
                    slice.col = 1,
                    ABA_coronal = 34,
                    ABA_color = "left",
                    rotation = 0,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 8600,
                    y.move = 47000,
                    pt.size = 0.3,
                    plot.axis = F)

slice51 = match_CCF(slice = 51,
                    slice.col = 1,
                    ABA_coronal = 25,
                    ABA_color = "left",
                    rotation = 175,
                    x.scale = 8.9,
                    y.scale = 8.9,
                    x.move = 12000,
                    y.move = 43000,
                    pt.size = 0.3,
                    plot.axis = F)
slice52 = match_CCF(slice = 52,
                    slice.col = 1,
                    ABA_coronal = 24,
                    ABA_color = "left",
                    rotation = 0,
                    x.scale = 9,
                    y.scale = 9,
                    x.move = 10500,
                    y.move = 43500,
                    pt.size = 0.3,
                    plot.axis = F)
slice53 = match_CCF(slice = 53,
                    slice.col = 1,
                    ABA_coronal = 28,
                    ABA_color = "left",
                    rotation = 82,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = 11500,
                    y.move = 40500,
                    pt.size = 0.3,
                    plot.axis = F)



slice_region = do.call(rbind, lapply(paste0("slice", c(8:14, 23:53)), get))

slice_region = slice_region[colnames(merFISH_integrated),]

# add pain: 54:60
slice54 = match_CCF(slice = 54,
                    slice.col = 1,
                    ABA_coronal = 33,
                    ABA_color = "left",
                    rotation = 1,
                    x.scale = 9,
                    y.scale = 9,
                    x.move = 12500,
                    y.move = 41000,
                    flip.x = T,
                    pt.size = 0.3,
                    plot.axis = F)

slice55 = match_CCF(slice = 55,
                    slice.col = 1,
                    ABA_coronal = 32,
                    ABA_color = "left",
                    rotation = -123,
                    x.scale = 8.5,
                    y.scale = 9,
                    x.move = 11000,
                    y.move = 43000,
                    flip.x = T,
                    pt.size = 0.3,
                    plot.axis = F)

slice56 = match_CCF(slice = 56,
                    slice.col = 3,
                    ABA_coronal = 33,
                    ABA_color = "left",
                    rotation = 183,
                    x.scale = 9,
                    y.scale = 9,
                    x.move = 12000,
                    y.move = 40500,
                    pt.size = 0.3,
                    plot.axis = F)

slice57 = match_CCF(slice = 57,
                    slice.col = 1,
                    ABA_coronal = 32,
                    ABA_color = "all",
                    rotation = 102,
                    x.scale = 7.5,
                    y.scale = 7.5,
                    x.move = -1500,
                    y.move = 41500,
                    pt.size = 0.3,
                    plot.axis = F)

slice58 = match_CCF(slice = 58,
                    slice.col = 1,
                    ABA_coronal = 37,
                    ABA_color = "all",
                    rotation = 8,
                    x.scale = 8,
                    y.scale = 8,
                    x.move = -4500,
                    y.move = 44000,
                    pt.size = 0.3,
                    plot.axis = F)

slice59 = match_CCF(slice = 59,
                    slice.col = 1,
                    ABA_coronal = 37,
                    ABA_color = "all",
                    rotation = 88,
                    x.scale = 7.5,
                    y.scale = 7.5,
                    x.move = -500,
                    y.move = 42000,
                    pt.size = 0.3,
                    plot.axis = F)

slice60 = match_CCF(slice = 60,
                    slice.col = 1,
                    ABA_coronal = 39,
                    ABA_color = "all",
                    rotation = 4,
                    x.scale = 8.9,
                    y.scale = 8.9,
                    x.move = -2600,
                    y.move = 46500,
                    pt.size = 0.3,
                    plot.axis = F)

new_pain_slice = do.call(rbind, lapply(paste0("slice", c(54:60)), get))

slice_region = merFISH_integrated@meta.data
slice_region = rbind(slice_region[!rownames(slice_region) %in% rownames(new_pain_slice), ], new_pain_slice)

slice_region = slice_region[colnames(merFISH_integrated),]


merFISH_integrated@meta.data = slice_region

merFISH_integrated$CellType = as.character(merFISH_integrated$L3_cluster)
merFISH_integrated$CellType[merFISH_integrated$CellType %in% levels(merFISH_integrated$L3_cluster)[1:7]] <- "Ext"
merFISH_integrated$CellType[merFISH_integrated$CellType %in% levels(merFISH_integrated$L3_cluster)[8:12]] <- "Inh"

saveRDS(merFISH_integrated, "./RDS/merFISH_final.1101.RDS")

