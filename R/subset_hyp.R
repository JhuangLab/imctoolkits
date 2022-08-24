#' @title subset_hyp
#' @name subset_hyp using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- subset_hyp(sce, cell_number = 2000)
#' sce <- subset_hyp(sce, cells = cells)
#'
subset_hyp <- function(sce, cell_number = 300000, force = TRUE){
  if(length(metadata(sce)) > 0 && force == FALSE){
    message("WARNING: We suggested run subset_hyp at very begin of analysis.\n
            All cells are returned.")
    return(sce)
  }else{
    if(!is.null(cell_number)){
      if(cell_number > ncol(sce)){
        message("WARNING: The input cell number is larger than all cell numbers. \n
                All cells will be returned.")
      }else{
        i <- sample(ncol(sce), size = cell_number, replace = FALSE)
        sce <- sce[, i]
      }
    }
    return(sce)
  }
}

#' @title subset_by_cell_ids
#' @name subset_by_cell_ids using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param cell_ids integer vector, c(45, 48, 52, 7).
#' @param include in or not in.
#' @return a singcellexperiment object
#' @export
#' @examples
#'
subset_by_cell_ids <- function(sce, cell_ids, include = TRUE){
  fil <- colData(sce)$cell_id %in% cell_ids
  if(include){
    sce_fil <- sce[, fil]
  }else{
    sce_fil <- sce[, !fil]
  }
  return(sce_fil)
}

#' @title subset_by_cluster_id
#' @name subset_by_cluster_id using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param fastpg data frame, a singcellexperiment object.
#' @param cluster_ids integer vector, c(45, 48, 52, 7).
#' @return a singcellexperiment object
#' @export
#'
#' @examples
subset_by_cluster_id <- function(sce, fastpg, cluster_ids = c(45, 48, 52, 7)){
  fil <- fastpg$cluster_id %in% cluster_ids
  sce_fil <- sce[, fastpg$cell_id[fil]]
  return(sce_fil)
}

#' @title subset_by_feature
#' @name subset_by_feature using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param feture_name data frame, a singcellexperiment object.
#' @param features integer vector, c(45, 48, 52, 7).
#' @return a singcellexperiment object
#' @export
#'
#' @examples
subset_by_feature <- function(sce, feature_name = "cell_type", features = c("Unknown"), include = TRUE){
  fil <- colData(sce)[, feature_name] %in% features
  if(include){
    sce_fil <- sce[, fil]
  }else{
    sce_fil <- sce[, !fil]
  }
  return(sce_fil)
}

#' @title subset_by_marker_heatmap_exp
#' @name subset_by_marker_heatmap_exp using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param marker string, "CD45".
#' @param fastpg data frame, a singcellexperiment object.
#' @param dat data frame, for heatmap plot.
#' @param cutoff double, a singcellexperiment object.
#' @param direction string, up or down.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- subset_by_marker_heatmap_exp(sce, cell_number = 2000)
#' sce <- subset_by_marker_heatmap_exp(sce, cells = cells)
#'
subset_by_marker_heatmap_exp <- function(sce, marker = "CD45", fastpg, dat, cutoff, direction = "up"){
  cd45p_fil <- dat[[marker]] >= cutoff
  if(direction == "up"){
    cd45p_cluster <- rownames(dat)[cd45p_fil] %>% str_replace("c", "")
    fil <- fastpg$cluster_id %in% cd45p_cluster
    sce_cd45p <- sce[, fastpg$cell_id[fil]]
    return(sce_cd45p)
  }
  if(direction == "down"){
    cd45n_cluster <- rownames(dat)[!cd45p_fil] %>% str_replace("c", "")
    fil <- fastpg$cluster_id %in% cd45n_cluster
    sce_cd45n <- sce[, fastpg$cell_id[fil]]
    return(sce_cd45n)
  }
}
