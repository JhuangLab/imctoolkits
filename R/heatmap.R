#'
#' Visualization heatmap of data of sce, part of codes are modified from Cytotree R package.
#' Thanks Yuting Dai.
#'
#' @name heatmap_anno
#'
#' @param sce A sce object
#' @param assay_name vector. Colors used in heatmap.
#' @param feature vector. markers to plot on the heatmap
#' @param scale character. Whether the values should be centered and scaled in either
#'    the row direction or the column direction, or none. Corresponding values are
#'    "row", "column" and "none"
#' @param col_fun logical. Whether columns should be clustered
#' @param ... options to pass on to the \code{\link[ComplexHeatmap]{ComplexHeatmap}} function.
#'
#' @export
#' @return ComplexHeatmap Heatmap figure
heatmap_anno <- function(sce, assay_name = "log2counts_censored", feature = "som.id",
                         scale = T, col_fun = NULL, ...){
  colors <- jhtools::show_me_the_colors("single_cell1")
  heatmap_data_scale <- fetch_heatmap_dat(sce, assay = assay_name, cluster.id = feature,
                                          scale = scale)
  cluster_prop <- heatmap_data_scale$cluster_prop
  row_ha = ComplexHeatmap::rowAnnotation(cluster_prop = anno_barplot(cluster_prop[ ,"cluster_proportion"],
                                                                     gp = gpar(fill = colors)))
  heatmap_plot <- heatmap_data_scale$heatmap_dat %>% data.matrix()
  if(is.null(col_fun)){
    col_fun = circlize::colorRamp2(c(-4, 0, 1, 2, 4), c("#94C6E5", "white", "#F4AB89", "#BD2333", "#6D001E"))
  }
  ht <- ComplexHeatmap::Heatmap(heatmap_plot, col = col_fun,
                                right_annotation = row_ha,
                                row_labels = rownames(heatmap_plot), ...)
  return(ht)
}

#'
#' Visualization heatmap of data of sce
#'
#' @name heatmap_anno_combined
#'
#' @param sce A sce object
#' @param assay_name vector. Colors used in heatmap.
#' @param feature vector. markers to plot on the heatmap
#' @param scale character. Whether the values should be centered and scaled in either
#'    the row direction or the column direction, or none. Corresponding values are
#'    "row", "column" and "none"
#' @param col_fun logical. Whether columns should be clustered
#' @param ... options to pass on to the \code{\link[ComplexHeatmap]{ComplexHeatmap}} function.
#'
#' @export
#' @return ComplexHeatmap Heatmap figure
#' @examples
#' \dontrun{
#' all_anno <- readr::read_rds("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/rds/combined/all_anno.rds")
#' pp <- c(0, 0.117, 0.25, 0.863,
#'          0.255, 0.001, 0.44, 0.98,
#'         0.705, 0.097, 0.3, 0.888,
#'         0.705, 0.097, 0.3, 0.888,
#'         0.85, 0.63, 0.02, 0.1)
#' pp <- matrix(pp, byrow = T, nrow = 5)
#' source("/cluster/home/yjliu_jh/projects/hyperion/code/plot/heatmap/main_heatmap_function.R")
#' heatmap_func(all_anno, plots_positions = pp)
#' }
heatmap_anno_combined <- function(sce, keycol = "cell_type10",
                                  nbhood_result,
                                  comm_clust,
                                  plots_positions
){
  # generate commonly used data
  meta_all <- metadata(sce)$sampleinfo
  coldata <- as.data.frame(colData(sce))
  # generate heatmap matrix
  c_mat_lst <- imctoolkits::fetch_heatmap_dat(sce, cluster.id = keycol)
  c_mat <- c_mat_lst[[1]]
  cluster_prop <- c_mat_lst[[2]]
  c_mat <- c_mat[!rownames(c_mat) %in% c("Unknown", "negative"),  ]
  # the ordering functions may change according to need, currently it's hclust
  cl_cell <- hclust(dist(c_mat, method = "euclidean"), method = "complete")
  label_cell <- cutree(cl_cell, k = 9, h = NULL) ## final order will include meta-clusters of cells
  order_cell <- names(label_cell)[order.dendrogram(as.dendrogram(cl_cell))]
  heatmap_dat <- c_mat[rev(order_cell), ]
  # colors setting
  colors <- jhtools::show_me_the_colors("single_cell1")
  colors2 <- jhtools::show_me_the_colors("single_cell2")
  names(colors2) <- 1:40
  # cell amount annotation
  col_fun = circlize::colorRamp2(c(-1, 0, 1, 2, 3), c("#94C6E5", "white", "#F4AB89", "#BD2333", "#6D001E"))
  prop_data <- cluster_prop[ ,3]
  names(prop_data) <- cluster_prop[ ,1]
  prop_data <- prop_data[match(rev(order_cell), names(prop_data))]
  cell_amount_anno = ComplexHeatmap::rowAnnotation(frac = anno_barplot(prop_data,
                                                                       gp = gpar(fill = colors)))
  # the main heatmap
  p0 <- Heatmap(heatmap_dat,
                cluster_rows = F, heatmap_legend_param = list(title = "exp"),
                show_row_names = F, col = col_fun, show_heatmap_legend = F,
                show_row_dend = F, show_column_dend = F
  )
  p0 <- attach_annotation(p0, cell_amount_anno, side = "right")
  grob <- grid.grabExpr(draw(p0))
  # ========= other annotations =========
  # the annotation heatmap 1
  coldat <- left_join(coldata, meta_all[, c("sample_id", "base_excision_eval")], by = "sample_id")
  coldat <- coldat[!coldat[, keycol] %in% c("Unknown", "negative"), ]
  group_dat <- coldat %>% group_by(across(keycol)) %>% mutate(allcount = n()) %>%
                group_by(across(c(keycol, "base_excision_eval", "allcount"))) %>% summarise(count = n()) %>%
                na.omit() %>% mutate(frac = count / allcount)
  group_dat$cell_type_order <- factor(group_dat[, keycol], levels = order_cell)

  # plot
  suppressWarnings(pcirc <- ggplot(group_dat, aes(y = cell_type_order, x = base_excision_eval)) +
                     geom_point(aes(color = cell_type_order, alpha = frac, size = count)) +
                     scale_color_manual(values = colors) +
                     scale_size(range = c(1, 15), name = "Fraction of clinical in celltype group") +
                     theme(
                       axis.text.y = element_text(color = colors),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       axis.line = element_line(colour = "black"),
                       legend.position = "none",
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank()
                     ))

  # load the data
  # nbhood_result <- read_csv(anno_1, col_types = cols())
  # comm_clust <- readr::read_rds(anno_2, col_types = cols())

  # the barplot of neighborhood
  coldat2 <- left_join(coldata, nbhood_result[, c("cell_id", "neighborhood10")], by = "cell_id")
  coldat2 <- coldat2[!coldat2[, keycol] %in% c("Unknown", "negative"), ]
  group_nbhood <- coldat2 %>% group_by(across(keycol), neighborhood10) %>% summarise(count = n()) %>% na.omit()
  group_nbhood$cell_type_order <- factor(group_nbhood[, keycol], levels = order_cell)

  # plot
  pbar <- ggplot(group_nbhood, aes(x = cell_type_order, y = count, fill = as.factor(neighborhood10))) +
    geom_bar(stat = "identity", show.legend = TRUE) +
    scale_fill_manual(values = brewer.pal(n = 10, name = "Set3")) +
    coord_flip() +
    theme(
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    guides(fill = guide_legend(title = "TME clusters"))

  # the barplot of community
  pcm <- ggplot(comm_clust, aes(cell_type_order, V1, fill = cluster)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = colors2) +
    coord_flip() +
    theme_classic() +
    labs(y = "Average cells per community", fill = "Clusters") +
    theme(axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_blank())
  # all legends can be modified by Legend() and arranged by packLegend()
  # here only one legends is manipulated in that way
  heat_legend <- Legend(title = "exp", at = -1:3, col_fun = col_fun)
  pp <- plots_positions
  # combination of sub-plots

  p4 <- ggdraw() +
    draw_plot(pcirc, pp[1, 1], pp[1, 2], pp[1, 3], pp[1, 4]) +   ## x, y, width, height
    draw_plot(grob, pp[2, 1], pp[2, 2], pp[2, 3], pp[2, 4])+
    draw_plot(pcm, pp[3, 1], pp[3, 2], pp[3, 3], pp[3, 4]) +
    draw_plot(grid.grabExpr(draw(heat_legend)), pp[5, 1], pp[5, 2], pp[5, 3], pp[5, 4])

  ggsave(glue("{out_dir}/annotated_heatmap_community.pdf"), p4, height = 7, width = 13)

  p5 <- ggdraw() +
    draw_plot(pcirc, pp[1, 1], pp[1, 2], pp[1, 3], pp[1, 4]) +   ## x, y, width, height
    draw_plot(grob, pp[2, 1], pp[2, 2], pp[2, 3], pp[2, 4])+
    draw_plot(pbar, pp[4, 1], pp[4, 2], pp[4, 3], pp[4, 4]) +
    draw_plot(grid.grabExpr(draw(heat_legend)), pp[5, 1], pp[5, 2], pp[5, 3], pp[5, 4])
  ggsave(glue("{out_dir}/annotated_heatmap_nbhood_clus.pdf"), p5, height = 7, width = 13)

}



#'
#' Visualization heatmap of data of sce, part of codes are modified from Cytotree R package.
#' Thanks Yuting Dai.
#'
#' @name plot_heatmap
#'
#' @param sce A sce object
#' @param assay_name vector. Colors used in heatmap.
#' @param feature vector. markers to plot on the heatmap
#' @param scale character. Whether the values should be centered and scaled in either
#'    the row direction or the column direction, or none. Corresponding values are
#'    "row", "column" and "none"
#' @param col_fun logical. Whether columns should be clustered
#' @param ... options to pass on to the \code{\link[ComplexHeatmap]{ComplexHeatmap}} function.
#'
#' @export
#' @return ComplexHeatmap Heatmap figure
#'
plot_heatmap <- function(tiffs, res, sce, filenam){

  res <- res %>% dplyr::filter(group_by %in% tiffs)
  n_tif <- length(unique(res$group_by))

  res <- res %>% group_by(from_label, to_label) %>% summarise(n_interact = sum(sigval == 1),n_avoid = sum(sigval == -1)) %>% ungroup()

  res <- res %>% dplyr::mutate(pct_interact = n_interact/n_tif,
                               pct_avoid = n_avoid/n_tif) %>%
    dplyr::mutate(n = ifelse(n_interact >= n_avoid, n_interact, n_avoid),
                  pct = ifelse(pct_interact >= pct_avoid, pct_interact, -pct_avoid))

  # pearson correlation of cell type proportions
  dt_cell1 <- colData(sce)[, c("cell_id", "sample_tiff_id", "cell_type10")] %>%
    as_tibble %>% dplyr::rename(cell_type = cell_type10) %>% group_by(sample_tiff_id,cell_type) %>%
    summarise(n_cell = n()) %>% ungroup()

  dt_cell2 <- colData(sce)[, c("cell_id", "sample_tiff_id", "cell_type10")] %>%
    as_tibble %>% dplyr::rename(cell_type = cell_type10) %>%
    group_by(sample_tiff_id) %>% summarise(n_total = n()) %>% ungroup()

  dt_cell <- left_join(dt_cell1, dt_cell2, by = "sample_tiff_id") %>%
    dplyr::mutate(prop = n_cell/n_total) %>% dplyr::select(-c(n_cell, n_total)) %>%
    pivot_wider(names_from = cell_type, values_from = prop, values_fill = 0)

  mat <- dt_cell %>% as.data.frame() %>% `rownames<-`(.$sample_tiff_id) %>% .[,-1] %>% as.matrix

  m_cor <- cor(mat)
  dt_cor <- as_tibble(m_cor, rownames = "from_label") %>% pivot_longer(!("from_label"),names_to = "to_label",values_to = "cor")

  # plot
  pt <- left_join(res, dt_cor, by = c("from_label","to_label")) %>% replace(is.na(.), 0)

  p <- ggplot(pt, aes(from_label, to_label)) +
    geom_tile(aes(fill = cor)) +
    geom_point(aes(color = pct, size = n)) +
    scale_size(range = c(1, 5), name = "Sig Interaction/\nAvoidance Counts") +
    scale_fill_gradient2(low = alpha("blue",0.3), high = alpha("red",0.3), limit = c(-1,1), name="Coorelation of Cell\nPercentage Per ROIs") +
    scale_color_gradient2(low = "blue", high = "red", limit = c(-1,1), name="Average Interaction\nPercentage Among ROIs") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  p$data$from_label <- factor(p$data$from_label, levels = c("Bcell","CD4T","CD8T","DC","M_MDSC","Monocyte","Macrophage_CD163p","Macrophage_HLADRp","Macrophage_HLADRp_CD163p","Macrophage_HLADRn_CD163n","MMT_HLADRp","MMT_HLADRp_CD163p","MMT_HLADRn_CD163n","mCAF","CAF_col1p","Vimentin_cell","Epithelial_tumor","Ki67_cell","Epithelial_normal","PSC","Endothelial","B7_H4_cell","CCR6_cell"))
  p$data$to_label <- factor(p$data$to_label, levels = c("Bcell","CD4T","CD8T","DC","M_MDSC","Monocyte","Macrophage_CD163p","Macrophage_HLADRp","Macrophage_HLADRp_CD163p","Macrophage_HLADRn_CD163n","MMT_HLADRp","MMT_HLADRp_CD163p","MMT_HLADRn_CD163n","mCAF","CAF_col1p","Vimentin_cell","Epithelial_tumor","Ki67_cell","Epithelial_normal","PSC","Endothelial","B7_H4_cell","CCR6_cell"))

  ggsave(glue::glue("{filenam}_interaction_heatmap.pdf"), p, height = 8, width = 10)
  readr::write_csv(pt, glue::glue("{filenam}_interaction_heatmap.csv"))
  return(p)
}
