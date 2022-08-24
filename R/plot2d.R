#' Visualization of 2D data of IMC, modified from Cytotree R package.
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
#' sce.file <- system.file("extdata/sce.rds", package = "imctoolkits")
#' sce <- readRDS(file = sce.file)
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


