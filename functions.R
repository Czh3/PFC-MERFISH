### some functions for spatial transcriptome


# gate points
library(gatepoints)
gate_points = function(obj, image = image){
  plot(obj@images[[image]]@coordinates[, 1:2], pch=20, cex=0.1)
  selectedPoints <- fhs(obj@images[[image]]@coordinates[, 1:2], mark = TRUE)
  c(selectedPoints)
}


# gene ontology enrichments
library(clusterProfiler)
require('org.Mm.eg.db')
library(enrichplot)
GO_enrichment = function(genes, all.genes){

  ego = enrichGO(gene     = genes,
                 keyType = "SYMBOL",
                 universe = all.genes,
                 OrgDb    = org.Mm.eg.db,
                 qvalueCutoff = 1,
                 pvalueCutoff = 1,
                 ont      = "ALL",
                 minGSSize = 2)
  ego = as.data.frame(ego[1:10, ])
  ego$Description = factor(ego$Description, levels = ego[order(ego$pvalue, decreasing = T), "Description"] ) 
  ggplot(ego, aes(-log10(pvalue), Description))+
    geom_bar(stat="identity") + ylab("") +
    theme_cowplot()
}

GO_enrichment_markers = function(markers, all.genes){
  enrich_GO <- compareCluster(gene~cluster, data=markers, ont="ALL",
                              universe = all.genes, keyType = "SYMBOL",
                              fun='enrichGO', OrgDb='org.Mm.eg.db')
  dotplot(enrich_GO)
}



Spatial_correlation = function(obj, number, res=0.05){
  set.seed(1)
  obj_metadata = obj@meta.data
  
  sample_region = function(image){
    
    images = sample(names(obj@images), 1)
    image.coordinates = obj@images[[images]]@coordinates
    x = sample(image.coordinates$x, 1)
    y = sample(image.coordinates$y, 1)

    n_bin_x = (max(image.coordinates$x) - min(image.coordinates$x)) * res
    n_bin_y = (max(image.coordinates$y) - min(image.coordinates$y)) * res
    cell_in = rownames(image.coordinates[image.coordinates$x > x & image.coordinates$x < x+n_bin_x  & image.coordinates$y > y & image.coordinates$y < y+n_bin_y,])
    
    as.vector(table(Idents(obj)[cell_in]))
  }
  
  cor_mat= sapply(1:number, sample_region)
  rownames(cor_mat) = levels(Idents(obj))
  cor_mat = cor_mat[,colSums(cor_mat) != 0]
  cor_mat = cor(t(cor_mat))
  pheatmap::pheatmap(cor_mat,
                     cluster_rows = T, cluster_cols = T, border_color = NA,clustering_method = "ward.D",
                     color = colorRampPalette(c("blue", "white", "red"))(30))
  return(cor_mat)
}

signature_plot1 = function(obj, genelist){
  obj = AddModuleScore(
    object = obj,
    features = list(genelist),
    ctrl = length(genelist),
    #replace = FALSE,
    slot = "data",
    name = 'sig'
  )
  p1 = FeaturePlot(obj, features = "sig1", raster=FALSE) + ggtitle("signature") +
    scale_colour_gradientn(colours = c( "grey95", "pink", "red", "darkred"))
  p1
}

signature_plot = function(object, gene.set, set.name){
  
  gene.set = gene.set[gene.set %in% rownames(object)]
  mean.exp <- colMeans(x = object@assays$RNA@data[gene.set, ], na.rm = TRUE)
  gene.set.name = set.name

  object@meta.data[gene.set.name] <- mean.exp
  max.cut = quantile(mean.exp, prob = 0.999)
  
  p1 = FeaturePlot(object = object, features = gene.set.name, pt.size = 0.1, max.cutoff = max.cut) &
    NoAxes()  &
    scale_color_gradientn(colours = c(rainbow(10)[5:10], "red")) 
    #scale_colour_gradientn(colours = c( "grey95", "pink", "red", "darkred"))
  p2 = SpatialFeaturePlot(object, features = gene.set.name, pt.size.factor=1, stroke=NA, slot="data", images = "image.2", max.cutoff = max.cut) +
    theme(legend.position = "right") +
    scale_fill_gradientn(colours = c(rainbow(10)[5:10], "red")) 
    #scale_fill_gradientn(colours = c("gray90", "pink", "red", "darkred"))
  
  gl.txt = data.frame(gene.set)
  gl.txt$x = rep(1:100, each=20)[1:length(gene.set)]
  gl.txt$y = rep(1:20, 100)[1:length(gene.set)]
  p3 = ggplot(gl.txt, aes(x, y, label = gene.set))+
    geom_text(size=2.5) +
    theme_void() + ggtitle(paste(set.name, " genes:"))
  
  print(p1|p2|p3)
}

plot_cell_location = function(object, color.by, color = NA, feature_metadata_path, fov, legend){
  if(is.na(color)){
    my_col = c("#e0829f",  "#7abd3d",  "#7269dd",  "#ccb635",  "#c458c7",  "#51c472",  "#cd4a97",  "#488f34",  "#844da3",  "#8e8c22",  "#5e74bf",  "#db8a2d",  "#5fa1d8",
               "#d74a2c",  "#3abec8",  "#d83b69",  "#58c29f",  "#c54b4b",  "#35815b",  "#c58dd6",  "#526a26",  "#9c4b6e",  "#85b26f",  "#b25621",  "#bbb564",  "#e6856c",
               "#878645",  "#9c5e3b",  "#d7a065",  "#866620", "#f12357", "#30a61a", "#4918a5", "#e5421d", "#3a47d9", "#e233c0", "#8f65f4", "#9527b7", "#743fd1", "#c052e4")
  } else {
    my_col = color
  }
  
  feature_metadata = read.csv(paste0(feature_metadata_path, "/feature_metadata.csv"), header = T, row.names = 1,  stringsAsFactors = F)
  feature_metadata = feature_metadata[feature_metadata$primary_fovID %in% fov, ]
  print(dim(feature_metadata))
  obj = object[, stringr::str_replace(colnames(object), "_.+", "") %in% rownames(feature_metadata)]
  obj@meta.data$boundaryX = feature_metadata[stringr::str_replace(colnames(obj), "_.+", ""), "boundaryX"]
  obj@meta.data$boundaryY = feature_metadata[stringr::str_replace(colnames(obj), "_.+", ""), "boundaryY"]
  obj_meta = obj@meta.data
  opar <- par(no.readonly = TRUE)
  par(mar = c(1, 1, 1, 15), xpd=TRUE)
  
  plot.new( )
  plot.window( xlim=c(min(obj_meta$centroid_1)-30, max(obj_meta$centroid_1)+30), 
               ylim=c(min(obj_meta$centroid_2)-30, max(obj_meta$centroid_2)+30) )
  
  colors = my_col[1:length(unique(object@meta.data[, color.by]))]
  names(colors) = sort(unique(object@meta.data[, color.by]))
  apply(obj_meta, 1, function(l) {
    x = l["boundaryX"][1]
    y = l["boundaryY"][1]
    x = strsplit(x, ";")[[1]]
    y = strsplit(y, ";")[[1]]
    polygon(x, y, border = colors[l[color.by]], lwd=.2, fillOddEven = FALSE, col = colors[l[color.by]])
  })
  
  if(legend){
    # Add a legend
    #n_col = ceiling(length(colors)/20)
    #legend("right", title=color.by, inset = c(-(n_col/8+0.1), 0), ncol = n_col,  border=NA, bty = "n",  x.intersp = 0.1,
    #       names(colors), fill=colors, cex=0.8, xpd = TRUE)
    
    n_col = ifelse(length(colors) > 30, 3, 1)
    legend("right", title=color.by, inset = c(-.2, 0), ncol = n_col,  border=NA, bty = "n",  x.intersp = 0.1,
           names(colors), fill=colors, cex=0.8, xpd = TRUE)
  }
  par(opar)
  on.exit(par(opar))
}






plot_expr_spots = function(object, genes, color = NA, data_path, fov, plot_cell = T){
  genes1 = genes
  
  geneName = read.csv(paste0(data_path, "/geneNames.csv"), header = F)$V1
  genes = match(genes, geneName)
  
  names(color) = genes
  
  barcode_metadata = data.table::fread(paste0(data_path, "/barcode_metadata.csv"), header = T, data.table = F)
  barcode_metadata = barcode_metadata[barcode_metadata$barcode_id %in% genes & barcode_metadata$fov_id %in% fov, ]
  
  feature_metadata = read.csv(paste0(data_path, "/feature_metadata.csv"), header = T, row.names = 1,  stringsAsFactors = F)
  feature_metadata = feature_metadata[feature_metadata$primary_fovID %in% fov, ]

  obj = object[, stringr::str_replace(colnames(object), "_.+", "") %in% rownames(feature_metadata)]
  obj@meta.data$boundaryX = feature_metadata[stringr::str_replace(colnames(obj), "_.+", ""), "boundaryX"]
  obj@meta.data$boundaryY = feature_metadata[stringr::str_replace(colnames(obj), "_.+", ""), "boundaryY"]
  obj_meta = obj@meta.data
  opar <- par(no.readonly = TRUE)
  par(mar = c(1, 1, 1, 1), bg="black")
  
  plot.new( )
  plot.window( xlim=c(min(obj_meta$centroid_1)-30, max(obj_meta$centroid_1)+30), 
               ylim=c(min(obj_meta$centroid_2)-30, max(obj_meta$centroid_2)+30) )
  
  if(plot_cell){
    apply(obj_meta, 1, function(l) {
      x = l["boundaryX"][1]
      y = l["boundaryY"][1]
      x = strsplit(x, ";")[[1]]
      y = strsplit(y, ";")[[1]]
      #polygon(x, y, border = "gray90", lwd=1, fillOddEven = FALSE, col = "black")
      polygon(x, y, border = "white", lwd=2, fillOddEven = FALSE, col = "black")
    })
  }
  
  barcode_metadata$color = color[as.character(barcode_metadata$barcode_id)]
  points(barcode_metadata$abs_position_1, barcode_metadata$abs_position_2,
         pch = 20,  cex = 0.1,
         col = barcode_metadata$color)
  
  legend("topleft",
         genes1, fill=color, cex=0.8, xpd = TRUE, text.col="white")
  
  on.exit(par(opar))
}




# spatial plot
MySpatialPlot = function(
  obj, 
  feature = NULL,
  group.by = NULL,
  images = images,
  pt.size.factor = 1.5,
  stroke=0.5
){
  y1 = min(obj@meta.data$centroid_1) - 1
  y2 = max(obj@meta.data$centroid_1) + 1
  x1 = min(obj@meta.data$centroid_2) - 1
  x2 = max(obj@meta.data$centroid_2) + 1
  
  lim1 = min(x1, y1)
  lim2 = max(x2, y2)
  
  SpatialPlot(obj, group.by = group.by,
              feature = feature,
              images = images,
              pt.size.factor = pt.size.factor, crop = T,
              stroke=stroke) &
    ylim(lim1, lim2) & xlim(lim1, lim2) &
    theme(legend.position = "right")
}

mySpatialDimPlot = function (object, group.by = NULL, features = NULL, images = NULL, 
                             cols = NULL, image.alpha = 1, crop = TRUE, slot = "data", 
                             min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL, cells.highlight.size.factor = c(1, 2),
                             cols.highlight = c("#DE2D26", "grey50"), facet.highlight = FALSE, 
                             label = FALSE, label.size = 5, label.color = "white", label.box = TRUE, 
                             repel = FALSE, ncol = NULL, combine = TRUE, pt.size.factor = 1.6, 
                             alpha = c(1, 1), stroke = 0.25, interactive = FALSE, do.identify = FALSE, 
                             identify.ident = NULL, do.hover = FALSE, information = NULL) 
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity", 
            call. = FALSE, immediate. = TRUE)
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(object = object, image = images[1], 
                             group.by = group.by, alpha = alpha))
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(ISpatialFeaturePlot(object = object, feature = features[1], 
                                 image = images[1], slot = slot, alpha = alpha))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, 
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index], 
                             data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index], 
                             data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && 
                                           is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature", 
                     no = "grouping")
      warning("'do.hover' requires only one ", type, ", using ", 
              features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells", 
              call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) * 
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, 
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && 
        is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, 
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      
      mySingleSpatialPlot = function (data, image, cols = NULL, image.alpha = 1, pt.alpha = NULL, 
                crop = TRUE, pt.size.factor = NULL, cells.highlight.size.factor = cells.highlight.size.factor, stroke = 0.25, col.by = NULL, 
                alpha.by = NULL, cells.highlight = NULL, cols.highlight = c("#DE2D26", 
                                                                            "grey50"), geom = c("spatial", "interactive", "poly"), 
                na.value = "grey50") 
      {
        geom <- match.arg(arg = geom)
        if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
          warning("Cannot find '", col.by, "' in data, not coloring", 
                  call. = FALSE, immediate. = TRUE)
          col.by <- NULL
        }
        col.by <- col.by %iff% paste0("`", col.by, "`")
        alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
        if (!is.null(x = cells.highlight)) {
          highlight.info <- Seurat:::SetHighlight(cells.highlight = cells.highlight, 
                                         cells.all = rownames(x = data), sizes.highlight = pt.size.factor, 
                                         cols.highlight = cols.highlight[1], col.base = cols.highlight[2])
          order <- highlight.info$plot.order
          data$highlight <- highlight.info$highlight
          col.by <- "highlight"
          levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), 
                                                     y = order))
          data <- data[order(data$ident), ]
        }

        plot <- ggplot(data = data, aes_string(x = colnames(x = data)[2], 
                                               y = colnames(x = data)[1], fill = col.by, alpha = alpha.by))
        plot <- switch(EXPR = geom, spatial = {
          if (is.null(x = pt.alpha)) {
            plot <- plot + Seurat:::geom_spatial(point.size.factor = ifelse(data$highlight == "Unselected", cells.highlight.size.factor[1], cells.highlight.size.factor[2]), 
                                        data = data, image = image, image.alpha = image.alpha, 
                                        crop = crop, stroke = stroke, )
          } else {
            plot <- plot + Seurat:::geom_spatial(point.size.factor = ifelse(data$highlight == "Unselected", cells.highlight.size.factor[1], cells.highlight.size.factor[2]),  
                                        data = data, image = image, image.alpha = image.alpha, 
                                        crop = crop, stroke = stroke, alpha = pt.alpha)
          }
          plot + coord_fixed() + theme(aspect.ratio = 1)
        }, interactive = {
          plot + Seurat:::geom_spatial_interactive(data = tibble(grob = list(GetImage(object = image, 
                                                                             mode = "grob"))), mapping = aes_string(grob = "grob"), 
                                          x = 0.5, y = 0.5) + geom_point(mapping = aes_string(color = col.by)) + 
            xlim(0, ncol(x = image)) + ylim(nrow(x = image), 
                                            0) + coord_cartesian(expand = FALSE)
        }, poly = {
          data$cell <- rownames(x = data)
          data[, c("x", "y")] <- NULL
          data <- merge(x = data, y = GetTissueCoordinates(object = image, 
                                                           qhulls = TRUE), by = "cell")
          plot + geom_polygon(data = data, mapping = aes_string(fill = col.by, 
                                                                group = "cell")) + coord_fixed() + theme_cowplot()
        }, stop("Unknown geom, choose from 'spatial' or 'interactive'", 
                call. = FALSE))
        if (!is.null(x = cells.highlight)) {
          plot <- plot + scale_fill_manual(values = cols.highlight)
        }
        if (!is.null(x = cols) && is.null(x = cells.highlight)) {
          if (length(x = cols) == 1 && (is.numeric(x = cols) || 
                                        cols %in% rownames(x = brewer.pal.info))) {
            scale <- scale_fill_brewer(palette = cols, na.value = na.value)
          }
          else if (length(x = cols) == 1 && (cols %in% c("alphabet", 
                                                         "alphabet2", "glasbey", "polychrome", "stepped"))) {
            colors <- DiscretePalette(length(unique(data[[col.by]])), 
                                      palette = cols)
            scale <- scale_fill_manual(values = colors, na.value = na.value)
          }
          else {
            scale <- scale_fill_manual(values = cols, na.value = na.value)
          }
          plot <- plot + scale
        }
        plot <- plot + NoAxes() + theme(panel.background = element_blank())
        return(plot)
      }
      
      plot <- mySingleSpatialPlot(data = cbind(coordinates, 
                                             data[rownames(x = coordinates), features[j], 
                                                  drop = FALSE]), image = image.use, image.alpha = image.alpha, 
                                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                                  features[j]
                                }
                                else {
                                  NULL
                                }, pt.alpha = if (!is.null(x = group.by)) {
                                  alpha[j]
                                }
                                else {
                                  NULL
                                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                                  "poly"
                                }
                                else {
                                  "spatial"
                                }, cells.highlight = highlight.use, cols.highlight = cols.highlight, 
                                pt.size.factor = pt.size.factor, stroke = stroke, 
                                cells.highlight.size.factor = cells.highlight.size.factor,
                                crop = crop)
      if (is.null(x = group.by)) {
        plot <- plot + scale_fill_gradientn(name = features[j], 
                                            colours = SpatialColors(n = 100)) + theme(legend.position = "top") + 
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight), 
                                                       yes = features[j], no = "highlight"), geom = if (inherits(x = image.use, 
                                                                                                                 what = "STARmap")) {
                                                         "GeomPolygon"
                                                       }
                              else {
                                "GeomSpatial"
                              }, repel = repel, size = label.size, color = label.color, 
                              box = label.box, position = "nearest")
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <- patchwork::wrap_plots(plots = plots, ncol = length(x = images))
  }
  else if (length(x = images == 1) && combine) {
    plots <- patchwork::wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

##
  
determinePCA = function(obj){
  std.PCA = c()
  for(i in 10){
    obj_random = obj
    obj_random@assays$RNA@scale.data = matrix(sample(obj_random@assays$RNA@scale.data), nrow = nrow(obj_random))
    colnames(obj_random@assays$RNA@scale.data) = colnames(obj@assays$RNA@scale.data)
    rownames(obj_random@assays$RNA@scale.data) = rownames(obj@assays$RNA@scale.data)
    
    obj_random = RunPCA(obj_random, npcs = 50, verbose = FALSE)
    std.PCA = c(obj_random@reductions$pca@stdev[1])
  }
  std.PCA = median(std.PCA)
  
  ElbowPlot(obj, ndims=50) +
    geom_hline(yintercept = std.PCA, color="red")
  
  num.PCA = sum(obj@reductions$pca@stdev >= std.PCA)
  message(paste0("The principle components to choose:", num.PCA))
  return(std.PCA)
}

# detect doublets
cell_purity = function(obj, 
                       reduction = reduction,
                       dims = dims){
  obj = Seurat::FindNeighbors(obj, 
                              k.param = 10,
                              reduction=reduction,
                              dims = dims,
                              return.neighbor = T)
  
  purity = parallel::mclapply(colnames(obj), function(i){
    a = table(Idents(obj)[Seurat::TopNeighbors(obj@neighbors$RNA.nn, i, n = 10)])
    max(a)/sum(a)
  }, mc.cores = 10)
  
  purity = do.call(c, purity)
  return(purity)
}

find_cluster_marker = function(obj){
  makers = lapply(levels(obj), function(i) {
    j = seq(as.numeric(i)-2, as.numeric(i)+2)
    j = j[j %in% levels(obj)]
    res = FindMarkers(obj, ident.1 = i, ident.2 = j, only.pos = T)
    res$cluster = i
    res$gene = rownames(res)
    res
  })
  makers = do.call(rbind, makers)
  makers
}




mySpatialPlot_example =  function(obj, L2_color, margin_x=2000, margin_y=1500, pt.size=1){
  #
  obj.p = obj@meta.data[obj@meta.data$slice == "29",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  
  p1 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = subcluster), size = pt.size) +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim-margin_y, ydim+margin_y) +
    theme_void()+
    theme(legend.position = "none") +
    scale_color_manual(values = L2_color) + coord_flip() 
    
  obj.p = obj@meta.data[obj@meta.data$slice == "28",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  
  p2 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = subcluster), size = pt.size) +
    theme_void()+
    scale_x_reverse() +
    theme(legend.position = "none") +
    scale_color_manual(values = L2_color) +
    xlim(xdim+margin_x, xdim-margin_x) + ylim(ydim-margin_y, ydim+margin_y)
  
  
  obj.p = obj@meta.data[obj@meta.data$slice == "30",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  
  p3 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = subcluster), size = pt.size) +
    scale_color_manual(values = L2_color) +
    theme_void()+
    theme(legend.position = "none") + scale_y_reverse() +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim+margin_y, ydim-margin_y)
  

  obj.p = obj@meta.data[obj@meta.data$slice == "31",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  
  p4 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = subcluster), size = pt.size) +
    scale_color_manual(values = L2_color) +
    theme_void()+
    theme(legend.position = "none") + scale_y_reverse() +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim+margin_y, ydim-margin_y)
  
  obj.p = obj@meta.data[obj@meta.data$slice == "27",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  
  p5 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = subcluster), size = pt.size) +
    scale_color_manual(values = L2_color) +
    theme_void()+
    theme(legend.position = "none") + scale_x_reverse() +
    xlim(xdim+margin_x, xdim-margin_x) + ylim(ydim-margin_y-400, ydim+margin_y-400)

  p1|p2|p3|p4|p5
}



mySpatialPlot_example2 =  function( obj, cells_highligth, margin_x=2000, margin_y=1500, pt.size=1){
  #
  obj.p = obj@meta.data[obj@meta.data$slice == "29",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  obj.p$hl = "No"
  obj.p$hl[obj.p$subcluster == cells_highligth] = "Yes"
  
  p1 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = hl), size = pt.size) +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim-margin_y, ydim+margin_y) +
    theme_void()+
    theme(legend.position = "none") +
    scale_color_manual(values = c("gray", "red2")) + coord_flip() 
  
  obj.p = obj@meta.data[obj@meta.data$slice == "28",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  obj.p$hl = "No"
  obj.p$hl[obj.p$subcluster == cells_highligth] = "Yes"
  
  p2 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = hl), size = pt.size) +
    theme_void()+
    scale_x_reverse() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("gray", "red2")) +
    xlim(xdim+margin_x, xdim-margin_x) + ylim(ydim-margin_y, ydim+margin_y)
  
  
  obj.p = obj@meta.data[obj@meta.data$slice == "30",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  obj.p$hl = "No"
  obj.p$hl[obj.p$subcluster == cells_highligth] = "Yes"
  
  p3 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = hl), size = pt.size) +
    scale_color_manual(values = c("gray", "red2")) +
    theme_void()+
    theme(legend.position = "none") + scale_y_reverse() +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim+margin_y, ydim-margin_y)
  
  
  obj.p = obj@meta.data[obj@meta.data$slice == "31",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  obj.p$hl = "No"
  obj.p$hl[obj.p$subcluster == cells_highligth] = "Yes"
  
  p4 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = hl), size = pt.size) +
    scale_color_manual(values = c("gray", "red2")) +
    theme_void()+
    theme(legend.position = "none") + scale_y_reverse() +
    xlim(xdim-margin_x, xdim+margin_x) + ylim(ydim+margin_y, ydim-margin_y)
  
  obj.p = obj@meta.data[obj@meta.data$slice == "27",]
  xdim = mean(obj.p$centroid_2)
  ydim = mean(obj.p$centroid_1)
  obj.p$hl = "No"
  obj.p$hl[obj.p$subcluster == cells_highligth] = "Yes"
  
  p5 = ggplot(obj.p, aes(centroid_2, centroid_1)) +
    geom_point(aes(color = hl), size = pt.size) +
    scale_color_manual(values = c("gray", "red2")) +
    theme_void()+
    theme(legend.position = "none") + scale_x_reverse() +
    xlim(xdim+margin_x, xdim-margin_x) + ylim(ydim-margin_y-400, ydim+margin_y-400)
  
  p1|p2|p3|p4|p5
}





## ABA
SpatialFeaturePlot_CCF = function(obj,
                                  slice, 
                                  feature,
                                 slot = "counts",
                                 max.cut = NULL,
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
  
  
  plot.data <- FetchData(
    object = obj,
    vars = c("reloc_x", "reloc_y", 'slice', "ABA_metaRegion", feature),
    slot = "counts"
  )
  
  plot.data = plot.data[plot.data$slice == slice, ]
  
  if(!is.null(max.cut)){
    plot.data[plot.data[,5] > max.cut, 5] <- max.cut
  }
  
  
  #rbPal <- c(rainbow(100)[31:100])
  #rbPal <- as.character(pals::viridis(20))
  rbPal <- colorRampPalette(c("gray90", "red2"))(50)
  plot.data$Col <- rbPal[as.numeric(cut(plot.data[,5],breaks = 70))]
  
  
  
  points(plot.data$reloc_x,
         plot.data$reloc_y,
         pch = 16, col = plot.data$Col, cex = pt.size, add = T)
  
  
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
  return(0)
}

SpatialDimPlot_CCF = function(obj,
                                  slice, 
                                  highlight,
                                  slot = "counts",
                                  max.cut = NULL,
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
  
  
  plot.data <- FetchData(
    object = obj,
    vars = c("reloc_x", "reloc_y", 'slice', "ABA_metaRegion", "subcluster"),
    slot = "counts"
  )
  
  plot.data = plot.data[plot.data$slice == slice, ]
  
  

  
  points(plot.data$reloc_x,
         plot.data$reloc_y,
         pch = 16, col = ifelse(plot.data$subcluster == highlight, "red", "gray90"), cex = pt.size, add = T)
  
  
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
  return(0)
}

