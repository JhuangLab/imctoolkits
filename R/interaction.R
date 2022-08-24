#' @title interaction_heatmap
#' @name interaction_heatmap using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
interaction_heatmap <- function(tiffs, res, sce, filename,
                                cell_types = c("Bcell","CD4T","CD8T","DC","M_MDSC","Monocyte","Macrophage_CD163p","Macrophage_HLADRp","Macrophage_HLADRn_CD163n","Macrophage_HLADRp_CD163p","MMT_CD163p","MMT_HLADRp","MMT_HLADRn_CD163n","MMT_HLADRp_CD163p","mCAF","CAF_col1p","Vimentin_cell","Epithelial_tumor","Ki67_cell","Epithelial_normal","PSC","Endothelial","B7_H4_cell","CCR6_cell")
                                ){
  res <- res %>% dplyr::filter(group_by %in% tiffs)
  n_tif <- length(unique(res$group_by))
  res <- res %>% group_by(from_label, to_label) %>%
    summarise(n_interact = sum(sigval == 1),n_avoid = sum(sigval == -1)) %>% ungroup()
  res <- res %>% dplyr::mutate(pct_interact = n_interact/n_tif,
                               pct_avoid = n_avoid/n_tif) %>%
                 dplyr::mutate(n = ifelse(n_interact >= n_avoid, n_interact, n_avoid),
                               pct = ifelse(pct_interact >= pct_avoid, pct_interact, -pct_avoid))
  # pearson correlation of cell type proportions
  dt_cell1 <- colData(sce)[, c("cell_id", "sample_tiff_id", "cell_type10")] %>%
    as_tibble %>% dplyr::rename(cell_type = cell_type10) %>% group_by(sample_tiff_id, cell_type) %>%
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
        scale_fill_gradient2(low = alpha("blue",0.3), high = alpha("red",0.3),
                             limit = c(-1,1), name="Coorelation of Cell\nPercentage Per ROIs") +
        scale_color_gradient2(low = "blue", high = "red", limit = c(-1,1),
                              name="Average Interaction\nPercentage Among ROIs") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p$data$from_label <- factor(p$data$from_label,
                              levels = cell_types)
  p$data$to_label <- factor(p$data$to_label,
                            levels = cell_types)
  ggsave(glue::glue("{filename}_interaction_heatmap.pdf"), p, height = 8, width = 10)
  readr::write_csv(pt, glue::glue("{filename}_interaction_heatmap.csv"))
  return(p)
}
