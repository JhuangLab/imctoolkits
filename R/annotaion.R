#' @title write_anno_by_cid
#' @name write_anno_by_cid using this to avoid function name confiction.
#' @param fastpg data frame, a singcellexperiment object.
#' @param cluster_anno a data frame.
#' @param fn_out string
#' @export
#'
write_anno_by_cid <- function(fastpg, cluster_anno, fn_out){
  anno <- inner_join(fastpg, cluster_anno, by = "cluster_id")
  cn <- c("cell_id", "cluster_id", "cell_type")
  anno <- dplyr::select(anno, all_of(cn)) %>% dplyr::filter(!is.na(cell_type))
  readr::write_csv(anno, fn_out, progress = F)
}

#' @title write_anno_sce_by_cid
#' @name write_anno_sce_by_cid using this to avoid function name confiction.
#' @param fastpg data frame, a singcellexperiment object.
#' @param cluster_anno a data frame.
#' @param cluster_id "fastpg.id"
#' @param fn_out string
#' @export
#' @examples
#' \dontrun{
#'  write_anno_sce_by_cid(sce, cluster_anno, fn_out = fn_out)
#' }
#' write_anno_sce_by_cid(sce, cluster_anno, fn_out = fn_out)
write_anno_sce_by_cid <- function(sce, cluster_anno, cluster_id = "fastpg.id", fn_out, unknow_include = F){
  fastpg <- colData(sce)[, c("cell_id", cluster_id)] %>% as.data.frame()
  colnames(fastpg) <- c("cell_id", "cluster_id")
  anno <- inner_join(fastpg, cluster_anno, by = "cluster_id")
  cn <- c("cell_id", "cluster_id", "cell_type")
  if(unknow_include){
    anno <- dplyr::select(anno, all_of(cn)) %>% dplyr::filter(!is.na(cell_type))
  }else{
    anno <- dplyr::select(anno, all_of(cn)) %>% dplyr::filter(!is.na(cell_type)) %>%
      dplyr::filter(!cell_type == "Unknown")
  }
  readr::write_csv(anno, fn_out, progress = F)
  return(invisible(anno))
}

#' @title update_feature
#' @name update_feature using this to avoid function name confiction.
#' @param sce data frame, a singcellexperiment object.
#' @param cell_anno a data frame.
#' @export
update_feature <- function(sce, cell_anno, feature = "cell_type"){
  if(!all(cell_anno$cell_id %in% colData(sce)$cell_id)){
    stop("The provided annotation data frame is wrong. Please check the cell id.")
  }else{
    j <- match(cell_anno$cell_id, colData(sce)$cell_id)
    j <- j[!is.na(j)]
    colData(sce)[j, feature] <- cell_anno[[feature]]
    return(sce)
  }
}

#' @title merge_cell_type
#' @name merge_cell_type using this to avoid function name confiction.
#' @param sce data frame, a singcellexperiment object.
#' @param cell_anno a data frame.
#' @export
merge_cell_type <- function(sce, from = "cell_type", to = "cell_type10"){
  if( to == "cell_type"){
    stop("Cannot merge to cell_type")
  }else{
    coldat <- colData(sce) %>% as.data.frame()
    coldat <- coldat %>% mutate({{to}} := case_when(
      str_detect(!!sym(from), "^Epithelial_tumor") ~ "Epithelial_tumor",
      TRUE ~ !!sym(from)
    ))
    colData(sce)[[to]] <- coldat[[to]]
    return(sce)
  }
}
