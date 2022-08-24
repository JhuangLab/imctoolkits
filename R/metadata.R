#' @title init_metadata
#' @name init_metadata using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- init_metadata(sce, sample_id = "p155_pdac")
#' sce <- init_metadata(sce, sample_id = "p155_pdac", stype = "tumor")
#'
init_metadata <- function(sce, sample_id, stype = "tumor"){
  dat_tmp <- colData(sce) %>% as.data.frame()
  cell_id <- glue("{sample_id}@{dat_tmp$sample_id}@{dat_tmp$ObjectNumber}")
  dat_tmp <- dat_tmp %>% dplyr::rename(sample_tiff_id = sample_id) %>%
             mutate(cell_id = cell_id) %>% mutate(cell_type = "") %>%
             mutate(sample_id = sample_id) %>%
             mutate(downsample = 1) %>%
             mutate(stype = stype) %>% relocate(sample_id)
  rownames(dat_tmp) <- dat_tmp$cell_id
  colData(sce) <- DataFrame(dat_tmp)
  return(sce)
}

#' @title update_metadata
#' @name update_metadata using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- update_metadata(sce)
#' sce <- update_metadata(sce, sample_id = "p155_pdac", stype = "tumor")
#'
update_metadata <- function(sce, downsample = NULL, cell_type = NULL, sample_id = NULL,
                            stype = NULL, df = NULL){
  if(!is.null(sample_id)){
    colData(sce)$sample_id <- sample_id
  }
  if(!is.null(downsample)){
    colData(sce)$downsample <- downsample
  }
  if(!is.null(stype)){
    colData(sce)$stype <- stype
  }
  if(!is.null(df)){
    colData(sce) <- DataFrame(df)
  }
  if(!is.null(cell_type)){
    colData(sce)$cell_type <- cell_type
  }
  return(sce)
}

#' @title set_metadata_info
#' @name set_metadata_info using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- set_metadata_info(sce, sample_id = "p155_pdac")
#' sce <- set_metadata_info(sce, sample_id = "p155_pdac", stype = "tumor")
#'

set_metadata_info <- function(sce, sample_id, stype = "tumor"){
  if(length(metadata(sce)) == 0){
    metadata(sce)$cell_info <- colData(sce) %>% as.data.frame()
    metadata(sce)$cell_info$stype <- stype
    metadata(sce)$cell_info <- dplyr::rename(metadata(sce)$cell_info, sample_tiff_id = sample_id)
    metadata(sce)$cell_info$sample_id <- sample_id
    metadata(sce)$cell_info <- dplyr::relocate(metadata(sce)$cell_info, sample_id)
    metadata(sce)$cell_info$dowsample <- 1
    dat_tmp <- metadata(sce)$cell_info
    rownames(metadata(sce)$cell_info) <- glue("{dat_tmp$sample_tiff_id}@{dat_tmp$ObjectNumber}")
  }
  if(length(metadata(sce)) > 0 && length(sample_id) > 1 && length(stype)>1){
    metadata(sce)$cell_info$stype <- stype
    metadata(sce)$cell_info$sample_id <- sample_id
  }
  else{
    message("metadata is already setted. We cannot set again to avoid disorder. \n
            You maybe should use update_metadata function.")
  }
  return(sce)
}

#' @title update_metadata_info
#' @name update_metadata_info using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- update_metadata_info(sce)
#' sce <- update_metadata_info(sce, sample_id = "p155_pdac", stype = "tumor")
#'
update_metadata_info <- function(sce, sample_id = NULL, stype = NULL){
  if(length(metadata(sce)) != 0){
    col_dat <- colData(sce) %>% as.data.frame()
    col_dat_indicate <- paste0(col_dat$sample_id, col_dat$ObjectNumber)
    cell_info_indicate <- paste0(metadata(sce)$cell_info$sample_tiff_id, metadata(sce)$cell_info$ObjectNumber)
    fil <- cell_info_indicate %in% col_dat_indicate
    metadata(sce)$cell_info[!fil, "sample_id"] <- 0
    if(!is.null(sample_id)){
      metadata(sce)$cell_info$sample_id <- sample_id
    }
    if(!is.null(stype)){
      metadata(sce)$cell_info$stype <- stype
    }
  }else{
    message("metadata is not setted. Please run set_metadata function first.")
  }
  return(sce)
}

#' @title add_cell_ratio
#' @name add_cell_ratio using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
add_cell_ratio <- function(sce, celltypes = "cell_type10"){
  col_dat <- colData(sce) %>% as.data.frame()
  meta_data <- metadata(sce)$sampleinfo %>% dplyr::select(sample_id:stype2)
  col_dat <- col_dat %>% group_by(sample_id) %>% mutate(n_cell = n()) %>%
    group_by(!!sym(celltypes), .add = T) %>% mutate(cell_ratio = n()/n_cell )
  cell_ratio <- col_dat %>% ungroup() %>% dplyr::select(sample_id, !!sym(celltypes), cell_ratio) %>%
    dplyr::distinct() %>% pivot_wider(names_from = !!sym(celltypes), values_from = cell_ratio) %>%
    mutate(across(where(is.numeric), replace_na, 0))
  # cell_ratio_tmp <- dplyr::select(cell_ratio, -sample_id) %>% names()
  # meta_data <- meta_data %>% dplyr::select(- any_of(cell_ratio_tmp))
  # meta_data <- meta_data %>% left_join(cell_ratio, by = "sample_id")
  metadata(sce)$cell_ratio <- cell_ratio
  return(sce)
}

#' @title add_cell_per_mm
#' @name add_cell_per_mm using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
add_cell_per_mm <- function(sce, celltypes = "cell_type10"){
  dt <- colData(sce) %>%
    as.data.table %>%
    .[, stype2 := factor(stype2)]
  # calculate cell numbers per mm2
  dt_n <- dt[, .N, by = .(sample_id, cell_type10)]
  sp <- dt$sample_id %>% unique
  cell <- dt$cell_type10 %>% unique
  dt_full <- data.table(sample_id = rep(sp, each = length(cell)),
                        cell_type10 = rep(cell, length(sp)))
  dt1 <- merge(dt_full, dt_n, all.x = TRUE)
  dt1[is.na(N), N := 0]
  dt2 <- dt[, area_tiff := width_px * height_px] %>%
    .[, .(sample_id, sample_tiff_id, area_tiff, stype2)] %>%
    unique %>%
    .[, sum(area_tiff), by = .(sample_id, stype2)]
  pt <- merge(dt1, dt2, by = "sample_id") %>%
    .[, cells_per_mm2 := N/V1 * 1000000]
  cell_per_mm <- pt %>% dplyr::select(-any_of(c("N", "stype2", "V1"))) %>%
    tidyr::pivot_wider(names_from = cell_type10,
                       values_from = cells_per_mm2,
                       values_fill = 0)
  metadata(sce)$cell_per_mm <- cell_per_mm
  return(sce)
}

