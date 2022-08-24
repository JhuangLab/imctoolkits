#' Visualization of 2D data of IMC
#'
#' @name plot2d
#'
#' @param object A singlecellexperiment object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot_meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#' @importFrom stats aggregate
#'
#' @export
#' @return ggplot2 figure
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "jhuanglabHyperion")
#' sce <- readRDS(file = sce.file)
#'
#' # Default plot
#' plot2d(sce)
#'
#' # PCA plot
#' plot2d(sce, item.use = c("PC_1", "PC_2"))
#' plot2d(sce, item.use = c("PC_1", "PC_2"), color.by = "cluster.id")
#' plot2d(sce, item.use = c("PC_1", "PC_2"), color.by = "stage")
#' plot2d(sce, item.use = c("PC_2", "PC_3"), color.by = "stage")
#' plot2d(sce, item.use = c("PC_2", "PC_3"), color.by = "CD43",
#'        category = "numeric")
#' plot2d(sce, item.use = c("PC_2", "PC_3"), color.by = "CD43",
#'        category = "numeric")
#'
#' # tSNE plot
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"))
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"), color.by = "stage")
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"), color.by = "cluster.id",
#'        alpha = 0.5, main = "tSNE Plot")
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"), color.by = "cluster.id",
#'        alpha = 1, main = "tSNE Plot", show.cluser.id = TRUE)
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"), color.by = "CD43",
#'        category = "numeric", size = 3)
#' plot2d(sce, item.use = c("tSNE_1", "tSNE_2"), color.by = "stage")
#'
#' # Diffusion Map plot
#' plot2d(sce, item.use = c("DC_1", "DC_2"))
#' plot2d(sce, item.use = c("DC_1", "DC_2"), color.by = "stage")
#' plot2d(sce, item.use = c("DC_2", "DC_3"), color.by = "cluster.id",
#'        alpha = 0.5, main = "Diffusion Map Plot")
#' plot2d(sce, item.use = c("DC_2", "DC_3"), color.by = "cluster.id",
#'        alpha = 1, main = "Diffusion Map Plot", show.cluser.id = TRUE)
#' plot2d(sce, item.use = c("DC_1", "DC_2"), color.by = "CD43",
#'        category = "numeric", size = 3)
#'
#' # UMAP plot
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"))
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"), color.by = "stage")
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id",
#'        alpha = 0.5, main = "UMAP Plot")
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id",
#'        alpha = 1, main = "UMAP Plot", show.cluser.id = TRUE)
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD43",
#'        category = "numeric", size = 3)
#' plot2d(sce, item.use = c("UMAP_1", "UMAP_2"), color.by = "stage")
#'
#' # Marker Plot
#' plot2d(sce, item.use = c("CD43", "CD90"), color.by = "cluster.id")
#' plot2d(sce, item.use = c("CD34", "CD90"), color.by = "CD43",
#'        category = "numeric", size = 3)
#'
#' # Pseudotime
#' plot2d(sce, item.use = c("pseudotime", "CD43"), color.by = "stage")
#'
#'
#'
plot2d <- function(object,
                   assay = "log2counts_censored",
                   item.use = c("tSNE_", "tSNE_"),
                   color.by = "stype",
                   order.by = NULL,
                   size = 1,
                   alpha = 1,
                   category = "categorical",
                   show.cluser.id = FALSE,
                   show.cluser.id.size = 4,
                   main = "2D plot of hypersion",
                   plot.theme = NULL) {
  if(is.null(plot.theme)){
    plot.theme <-  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  }
  gray_blue <- c("#d3d3d3", "#1506ff")
  # update and fetch plot meta information
  plot_meta <- fetch_reducedims(object, assay = assay, verbose = FALSE)
  # check color.by parameter in plot_meta data.frame
  if(!color.by  %in% colnames(plot_meta)){
    stop(Sys.time(), ": ", color.by, " is not in features.")
  }
  # check item.use parameter in plot_meta data.frame
  if ( !all(item.use %in% colnames(plot_meta)) ) {
    if (all(item.use %in% paste0("tSNE_", 1:10))) {
      stop(Sys.time(), " item.use is not in plot_meta of sce, please run runTSNE first.")
    } else if (all(item.use %in% paste0("PCA_", 1:10))) {
      stop(Sys.time(), " item.use is not in plot_meta of sce, please run runFastPCA first.")
    } else if (all(item.use %in% paste0("DC_", 1:10))) {
      stop(Sys.time(), " item.use is not in plot_meta of sce, please run runDiffusionMap first.")
    } else if (all(item.use %in% paste0("UMAP_", 1:10))) {
      stop(Sys.time(), " item.use is not in plot_meta of sce, please run runUMAP first.")
    } else {
      stop(Sys.time(), " item.use is not in plot_meta of sce.")
    }
  }

  if (length(item.use) < 2) stop(Sys.time(), " item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[seq_len(2)]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(plot_meta))
  color.by.idx <- match(color.by, colnames(plot_meta))

  plot.x = plot.y = NULL

  plot_data <- data.frame(plot.x = plot_meta[, item.use.idx[1]],
                          plot.y = plot_meta[, item.use.idx[2]],
                          color.by = plot_meta[, color.by.idx])

  if ((length( unique(plot_data$color.by) ) > 256) & (category != "numeric")) {
    warning(Sys.time(), " color.by is categorical and has more than 256 elements. It will be used as numeric instead.")
    category = "numeric"
  }

  if (is.null(category)) {
    if (is.numeric(plot_data$color.by)) category="numeric" else category="categorical"
  }
  if (category == "categorical") {
    if (is.null(order.by)) {
      plot_data$color.by <- factor(plot_data$color.by)
    } else {
      plot_data$color.by <- factor(as.character(plot_data$color.by), levels = order.by)
    }
  } else if (category == "numeric") {
    if (!is.numeric(plot_data$color.by)) plot_data$color.by <- as.numeric(factor(plot_data$color.by))
  } else {
    warning(Sys.time(), " Unidentified parameters of category")
  }
  # plot
  gg <- ggplot(plot_data) + geom_point(aes(x=plot.x, y=plot.y, color = color.by), size = size, alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))
  gg <- gg + labs(color = color.by)
  if (show.cluser.id & (category == "categorical")) {
    pos <- aggregate(  plot_data[, seq_len(2)], list( pos = plot_data$color.by ), mean)

    for ( i in seq_along(pos$pos)) {
      gg <- gg + annotate(geom="text", x = pos$plot.x[i], y = pos$plot.y[i],
                          label = pos$pos[i],
                          size = show.cluser.id.size)
    }
  }
  return(gg)
}

#'
#' Visualization violin plot of sce
#'
#' @name plotViolin
#'
#' @param object A sce object
#' @param marker character. Markers used to plot
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot_meta, or it can be just "density" (the default value).
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param text.angle numberic. Text angle of the violin plot
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#' @importFrom stats aggregate
#'
#' @export
#' @return ggplot2 figure
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "jhuanglabHyperion")
#' sce <- readRDS(file = sce.file)
#'
#' plotViolin(sce, marker = "CD34")
#' plotViolin(sce, marker = "CD34", order.by = "pseudotime")
#'
#'
plotViolin <- function(object,
                       marker,
                       color.by = "cluster.id",
                       order.by = NULL,
                       size = 1,
                       text.angle = 0,
                       main = "Violin plot sce",
                       plot.theme = theme_bw()) {

  # update plot meta information
  plot_meta <- fetchPlotMeta(object, verbose = FALSE)

  if (missing(marker)) stop(Sys.time(), " marker is missing.")
  # check item.use parameter in plot_meta data.frame
  if (length(marker) > 1) {
    warning(Sys.time(), " marker has more than two elements. Only the first two will be used")
    marker <- marker[1]
  }
  if ( marker %in% colnames(plot_meta) ) {
    plot_meta <- data.frame(plot_meta, marker = plot_meta[which(object@meta.data$dowsample == 1), marker])
  } else {
    stop(Sys.time(), " marker name is not correct")
  }


  # check color.by parameter in plot_meta data.frame
  if ( !all(color.by %in% colnames(plot_meta)) ) stop(Sys.time(), " color.by is not in plot_meta of sce.")

  if (length(color.by) > 1) {
    warning(Sys.time(), " color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  color.by.idx <- match(color.by, colnames(plot_meta))

  marker.by = NULL

  plot_data <- data.frame(marker.by = plot_meta$marker,
                          color.by = plot_meta[, color.by.idx])

  if (length( unique(plot_data$color.by) ) > 128) {
    stop(Sys.time(), " color.by is categorical and has more than 128 elements.")
  }

  if (is.null(order.by)) {
    plot_data$color.by <- factor(plot_data$color.by)
  } else if (order.by == "pseudotime") {
    sub <- plot_meta[, c("pseudotime", color.by)]
    colnames(sub) <- c("pseudotime", "color.by.tag")
    sub <- aggregate(sub, list(color.by = sub$color.by.tag), mean)
    plot_data$color.by <- factor(as.character(plot_data$color.by), levels = sub$color.by.tag[order(sub$pseudotime)])
  }
  else {
    plot_data$color.by <- factor(as.character(plot_data$color.by), levels = order.by)
  }
  # plot
  gg <- ggplot(plot_data, aes(x = color.by, y= marker.by, fill = color.by)) + geom_violin(scale = "width")
  gg <- gg + plot.theme
  gg <- gg + stat_summary(fun.y=mean, geom="point", size = size, color="black")
  gg <- gg + labs(y = marker, x = color.by, title = paste0(main))
  gg <- gg + labs(fill = color.by)
  gg <- gg + theme(axis.text.x = element_text(angle = text.angle, hjust = 1, vjust = 1))
  return(gg)
}




#'
#' Visualization pie plot of cluster data of sce
#'
#' @name plotPieCluster
#'
#' @param object A sce object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param cex.size numeric. Size of the dot
#' @param size.by.cell.number logical. Whether to show size of cell number.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#' @return ggplot2 figure
#'
#' @export
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "jhuanglabHyperion")
#' sce <- readRDS(file = sce.file)
#'
#' # Runs only have more than two stages
#' plotPieCluster(sce, cex.size = 0.5)
#'
#' plotPieCluster(sce, item.use = c("PC_1", "PC_2"), cex.size = 0.5)
#' plotPieCluster(sce, item.use = c("PC_2", "PC_3"), cex.size = 0.5)
#'
#' plotPieCluster(sce, item.use = c("tSNE_1", "tSNE_2"), cex.size = 20)
#'
#' plotPieCluster(sce, item.use = c("DC_1", "DC_2"), cex.size = 0.5)
#'
#' plotPieCluster(sce, item.use = c("UMAP_1", "UMAP_2"), cex.size = 1)
#' plotPieCluster(sce, item.use = c("UMAP_1", "UMAP_2"), cex.size = 1)
#'
#'
plotPieCluster <- function(object,
                           item.use = c("PC_1", "PC_2"),
                           cex.size = 1,
                           size.by.cell.number = TRUE,
                           main = "2D pie plot of sce",
                           plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " object is missing")
  if (is.null(object@network)) stop(Sys.time(), " network is missing, please run runCluster first!")
  if (length(unique(object@meta.data$stage)) <= 1) stop(Sys.time(), " plotPieCluster only fits elements in stage over 2!")

  # update plot meta information
  plot_data <- fetchClustMeta(object, verbose = FALSE)


  if (length(item.use) < 2) stop(Sys.time(), " item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[seq_len(2)]
  }
  item.use.idx <- match(item.use, colnames(object@cluster))

  plot.cols <- paste0(unique(object@meta.data$stage), ".percent")

  pos.x = pos.y = cluster = cell.number.percent = NULL
  plot_data <- data.frame(plot_data,
                          pos.x = object@cluster[, item.use.idx[1]],
                          pos.y = object@cluster[, item.use.idx[2]])

  gg <- ggplot()
  if (size.by.cell.number) {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = cell.number.percent*cex.size),
                               data = plot_data, cols = plot.cols, color=NA) + coord_equal()
  } else {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = 0.1*cex.size),
                               data = plot_data, cols = plot.cols, color=NA) + coord_equal()
  }

  gg <- gg + plot.theme
  gg <- gg + labs(x = "", y = "", title = main)

  return(gg)

}

#'
#' Visualization of cluster data of sce
#'
#' @name plotCluster
#'
#' @param object An sce object
#' @param item.use character. Items use to 2D plot, axes x and y must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot_meta, or it can be just "density" (the default value).
#' @param size.by character. Size of the dot
#' @param order.by vector. Order of color theme.
#' @param size numeric. Size of the dot
#' @param alpha numberic. Transparency (0-1) of the dot, default is 1.
#' @param category character. numeric or categorical
#' @param show.cluser.id logical. Whether to show cluster id in the plot.
#' @param show.cluser.id.size numeric. Size of the cluster id.
#' @param main character. Title of the plot.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @import ggplot2
#' @return ggplot2 figure
#'
#' @export
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "jhuanglabHyperion")
#' sce <- readRDS(file = sce.file)
#'
#' plotCluster(sce)
#'
#' plotCluster(sce, item.use = c("PC_1", "PC_2"))
#' plotCluster(sce, item.use = c("PC_2", "PC_3"))
#' plotCluster(sce, item.use = c("PC_2", "PC_3"), color.by = "CD43", category = "numeric")
#' plotCluster(sce, item.use = c("PC_2", "PC_3"), color.by = "CD43", category = "numeric")
#'
#' plotCluster(sce, item.use = c("tSNE_1", "tSNE_2"))
#' plotCluster(sce, item.use = c("tSNE_1", "tSNE_2"), show.cluser.id = TRUE)
#'
#' plotCluster(sce, item.use = c("DC_1", "DC_2"))
#'
#' plotCluster(sce, item.use = c("UMAP_1", "UMAP_2"))
#'
#'
plotCluster <- function(object,
                        item.use = c("PC_1", "PC_2"),
                        color.by = "cluster",
                        size.by = "cell.number.percent",
                        order.by = NULL,
                        size = 1,
                        alpha = 1,
                        category = "categorical",
                        show.cluser.id = FALSE,
                        show.cluser.id.size = 4,
                        main = "2D plot of cluster in sce",
                        plot.theme = theme_bw()) {

  # update plot meta information
  plot_meta.data <- fetchClustMeta(object, verbose = FALSE)
  plot_meta.data <- cbind(plot_meta.data, object@cluster)

  # check item.use parameter in plot_meta data.frame
  if ( !all(item.use %in% colnames(plot_meta.data)) ) {
    stop(Sys.time(), " item.use is not in cluster data of sce.")
  }

  # check color.by parameter in plot_meta data.frame
  if ( !all(color.by %in% colnames(plot_meta.data)) ) stop(Sys.time(), " color.by is not in cluster data of sce.")

  # check size.by parameter in plot_meta data.frame
  if ( !all(size.by %in% colnames(plot_meta.data)) ) stop(Sys.time(), " size.by is not in cluster data of sce.")

  if (length(item.use) < 2) stop(Sys.time(), " item.use is less than two elements.")
  if (length(item.use) > 2) {
    warning(Sys.time(), " item.use has more than two elements. Only the first two will be used")
    item.use <- item.use[seq_len(2)]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  if (length(size.by) > 1) {
    warning(Sys.time(), " size.by has more than one elements. Only the first one will be used")
    size.by <- size.by[1]
  }

  item.use.idx <- match(item.use, colnames(plot_meta.data))
  color.by.idx <- match(color.by, colnames(plot_meta.data))
  size.by.idx <- match(size.by, colnames(plot_meta.data))

  plot.x = plot.y = cluster = cell.number.percent = NULL
  plot_data <- data.frame(plot.x = plot_meta.data[, item.use.idx[1]],
                          plot.y = plot_meta.data[, item.use.idx[2]],
                          color.by = plot_meta.data[, color.by.idx],
                          size.by = plot_meta.data[, size.by.idx])

  if ((length( unique(plot_data$color.by) ) > 256) & (category != "numeric")) {
    warning(Sys.time(), " color.by is categorical and has more than 256 elements. It will be used as numeric instead.")
    category = "numeric"
  }

  if (is.null(category)) {
    if (is.numeric(plot_data$color.by)) category="numeric" else category="categorical"
  }
  if (category == "categorical") {
    if (is.null(order.by)) {
      plot_data$color.by <- factor(plot_data$color.by)
    } else {
      plot_data$color.by <- factor(as.character(plot_data$color.by), levels = order.by)
    }
  } else if (category == "numeric") {
    if (!is.numeric(plot_data$color.by)) plot_data$color.by <- as.numeric(factor(plot_data$color.by))
  } else {
    warning(Sys.time(), " Unidentified parameters of category")
  }

  # plot
  gg <- ggplot(plot_data) + geom_point(aes(x=plot.x, y=plot.y, color = color.by, size = size*size.by), alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = item.use[1], y = item.use[2], title = paste0(main))
  gg <- gg + labs(color = color.by)

  if (show.cluser.id) {
    for ( i in seq_along(rownames(plot_data))) {
      plot.x.anno = plot_data$plot.x
      plot.y.anno = plot_data$plot.y
      gg <- gg + annotate(geom="text", x = plot.x.anno[i], y = plot.y.anno[i],
                          label = rownames(plot_data)[i],
                          size = show.cluser.id.size)
    }
  }

  return(gg)

}

#'
#' Visualization of cluster data of sce
#'
#' @name boxplot_clinical
#'
#' @param df An sce object
#' @param sel_group character. Items use to 2D plot, axes x and y must be numeric.
#' @param celltype character. Dot or mesh color by which character. It can be one of the column
#'     of plot_meta, or it can be just "density" (the default value).
#' @return ggplot2 figure
#'
#' @export
boxplot_clinical <- function(df = df, sel_group = "stage", celltype = "DC",
                             axistheme = axistheme1, compare = compare_list){
  colors <- list(color2 = c("#F1D1CB", "#ECA2AE"),
                 color3 = c("#A4D8D9", "#CCDCEB", "#B5BAD7"),
                 color6 = c("#A4D8D9", "#CCDCEB", "#B5BAD7", "#F1D1CB", "#ECA2AE", "#FCD092"))
  dot_colors <- list(dot_color2 = c("#E3A498", "#D9465D"),
                     dot_color3 = c("#56B5B7", "#99B7D8", "#7079B1"),
                     dot_color6 = c("#56B5B7", "#99B7D8", "#7079B1", "#E3A498", "#D9465D", "#F5B04E"))
  plot.theme <- theme_bw() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank()) +
                             theme(legend.position = "none")
  g_num <- length(unique(df[[sel_group]]))
  color_index <- case_when(g_num == 2 ~ 1,
                           g_num == 3 ~ 2,
                           g_num == 6 ~ 3,
                           TRUE ~ 1)
  dot_color_index <- case_when(g_num == 2 ~ 1,
                               g_num == 3 ~ 2,
                               g_num == 6 ~ 3)
  comparisons_n <- length(compare) + 2
  ydata <- df %>% dplyr::select(!!sym(celltype))
  yrange <- max(ydata) - min(ydata)
  p <- ggplot(df, aes_string(sel_group, celltype))+
    geom_boxplot(outlier.shape= NA, lwd = 0.2, color = "#000000", fill = colors[[color_index]]) +
    geom_jitter(shape = 16, size = 0.5, position = position_jitter(0.2), aes(color = !!sym(sel_group)))+
    scale_colour_manual(values = dot_colors[[dot_color_index]]) +
    plot.theme +
    axistheme +
    stat_compare_means(method = "t.test",
                       comparisons = compare,
                       hide.ns = TRUE,
                       bracket.size = 0.2,
                       vjust = 0.8,
                       aes(label = ..p.signif..))
  return(p)
}
