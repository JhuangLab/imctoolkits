#' @title fetch_reducedims
#' @name fetch_reducedims using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- fetch_reducedims(sce)
#' sce <- fetch_reducedims(sce)
#'}

fetch_reducedims <- function(sce, assay = "log2counts_censored", verbose = FALSE){
  # if(!"tsne" %in% names(SingleCellExperiment::reducedDims(sce))){
  #   stop(Sys.time(), " tsne is not in reducedDims(sce). Nothing to fetch.")
  # }
  if("dfm_trans" %in% names(reducedDims(sce))){
    i <- which("dfm_trans" == names(reducedDims(sce)))
    reducedDims(sce) <- reducedDims(sce)[-i]
    #dat <- dat %>% dplyr::select(-colnames(reducedDims(sce)$dfm_trans))
  }
  dat <- purrr::reduce(reducedDims(sce), cbind) %>% as.matrix() %>% as.data.frame()
  dat <- dat %>% tibble::rownames_to_column(var="cell_id")
  dat_exp <- t(assay(sce, assay)) %>% as.data.frame() %>% tibble::rownames_to_column(var="cell_id")
  dat_metadata <- colData(sce) %>% as.matrix()  %>% as.data.frame() %>% dplyr::left_join(dat, by = "cell_id") %>%
                  dplyr::left_join(dat_exp, by = "cell_id")
  return(dat_metadata)
}

#' @title fetch_plot_reducedims
#' @name fetch_plot_reducedims using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param marker string, "CD45".
#' @param iterm.use string vector, c("UMAP_1", "UMAP_2").
#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- fetch_plot_reducedims(sce)
#' sce <- fetch_plot_reducedims(sce)
#'}

fetch_plot_reducedims <- function(sce, assay = "logcounts", marker = "CD45", iterm.use = c("UMAP_1", "UMAP_2"), verbose = FALSE){
  # if(!"tsne" %in% names(SingleCellExperiment::reducedDims(sce))){
  #   stop(Sys.time(), " tsne is not in reducedDims(sce). Nothing to fetch.")
  # }
  dat <- purrr::reduce(reducedDims(sce), cbind) %>% as.matrix() %>% as.data.frame() %>% dplyr::select(iterm.use)
  if("dfm_trans" %in% names(reducedDims(sce))){
    dat <- dat %>% dplyr::select(-colnames(reducedDims(sce)$dfm_trans))
  }
  dat <- dat %>% tibble::rownames_to_column(var="cell_id")
  dat_exp <- t(assay(sce, assay)) %>% as.data.frame() %>% dplyr::select(marker) %>%
             tibble::rownames_to_column(var="cell_id")
  dat_metadata <- colData(sce) %>% as.matrix()  %>% as.data.frame() %>% dplyr::left_join(dat, by = "cell_id") %>%
    dplyr::left_join(dat_exp, by = "cell_id")
  return(dat_metadata)
}

#' @title fetch_heatmap_dat
#' @name fetch_heatmap_dat using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param assay string, "logcounts", "counts".
#' @param scale bool
#' @export
#' @return A list.
fetch_heatmap_dat <- function(sce, assay = "log2counts_censored", cluster.id = "phenograph.id",
                              scale = T, fil_markers = NULL, prefix = ""){
  sce_trans <- as.data.frame(t(assay(sce, assay))) %>% dplyr::mutate(ID = rownames(.))
  if(!is.null(fil_markers)){
    sce_trans <- dplyr::select(sce_trans, - all_of(fil_markers))
  }
  tib_tmp <- tibble(ID = rownames(colData(sce)), colData(sce)[, cluster.id]) %>% as.data.frame()
  colnames(tib_tmp) <- c("ID", "sID")
  sce_trans <- sce_trans %>%
    dplyr::left_join(tib_tmp, by = "ID") %>%
                   dplyr::select(-ID)
  counts_dat <- sce_trans %>% group_by(sID) %>% summarise(count_cells = n()) %>%
                dplyr::mutate(cluster_proportion = count_cells / sum(count_cells)) %>%
                as.data.frame()
  row.names(counts_dat) <- glue("{prefix}{counts_dat$sID}")
  if(scale){
    sce_heatmap <- sce_trans %>% group_by(sID) %>%
      summarise_if(is.numeric, mean, na.rm = TRUE) %>% tibble::column_to_rownames(var = "sID") %>% scale()
  }else{
    sce_heatmap <- sce_trans %>% group_by(sID) %>%
      summarise_if(is.numeric, mean, na.rm = TRUE) %>% tibble::column_to_rownames(var = "sID")
  }
  rownames(sce_heatmap) <- glue("{prefix}{rownames(sce_heatmap)}")
  return(list(heatmap_dat = sce_heatmap, cluster_prop = counts_dat))
}








