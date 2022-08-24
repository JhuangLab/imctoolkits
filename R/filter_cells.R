#' @title filter_cells
#' @name filter_cells using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'
#' sce <- filter_cells(sce, cluster.id = 1)
#' sce <- filter_cells(sce, cluster.id = c(1,2))
#' a <- filter_cells(sce, cluster.id = 35)

filter_cells <- function(sce, assay = "compcounts", cluster.id = NULL, phenograph.id = NULL, som.id = NULL, cell_id =NULL){
  stopifnot("input must be a SingleCellExperiment object" = "SingleCellExperiment" %in% class(sce))
  if(!is.null(cluster.id)){
    fil <- ! colData(sce)$cluster.id  %in% cluster.id
    msg <- glue("{Sys.time()}: Total {sum({!fil})} cells are filtered by cluster.id: {cluster.id}.")
    message(msg)
    sce <- sce[, fil]
  }
  if(!is.null(phenograph.id)){
    fil <- ! colData(sce)$phenograph.id  %in% phenograph.id
    msg <- glue("{Sys.time()}: Total {sum({!fil})} cells are filtered by phenograph.id: {phenograph.id}.")
    message(msg)
    sce <- sce[, fil]
  }
  if(!is.null(som.id)){
    fil <- ! colData(sce)$som.id  %in% som.id
    msg <- glue("{Sys.time()}: Total {sum({!fil})} cells are filtered by som.id: {som.id}.")
    message(msg)
    sce <- sce[, fil]
  }
  if(!is.null(cell_id)){
    fil <- ! colData(sce)$cell_id  %in% cell_id
    msg <- glue("{Sys.time()}: Total {sum({!fil})} cells are filtered by cell_id: {cell_id}.")
    message(msg)
    sce <- sce[, fil]
  }

  # Any cells with eccentricity larger than the threshold will be removed.
  # 'area_threshold' represents for the z-score of cell size in each ROI.
  # Increase the value to filter more cells, decrease to retain more cells.
  if (is.null(assay(sce, assay))){
    message("Please perform spillover correction first!")
  }else{
    norm <- assay(sce, assay)
    norm <- as.data.frame(colSums(log10(norm +1)))
    colnames(norm)[1] <- "Value"
    norm$cut <- cut(norm$Value, breaks = c(0, 2, 16, 100),
                    labels = c("too small", "OK", "too large"))
    norm$cell_id <- rownames(norm)
    log_sig_removed <- norm[norm$cut != "OK", ]
    c_for_filter <- as.data.frame(colData(sce))
    c_for_filter <- c_for_filter %>% group_by(sample_tiff_id) %>%
      mutate(zscore_area = (area - mean(area)) / sd(area))
    dim(c_for_filter[c_for_filter$eccentricity >= eccentricity_threshold, ])
    dim(c_for_filter[c_for_filter$zscore_area < area_threshold, ])
    dim(c_for_filter[c_for_filter$zscore_area < area_threshold, ])
    c_removed <- rbind(c_for_filter[c_for_filter$eccentricity >= eccentricity_threshold, ],
                       c_for_filter[c_for_filter$zscore_area < area_threshold, ])
    remove_celllist <- unique(c(c_removed$cell_id, log_sig_removed$cell_id))
    sce <- sce[, colData(sce)$cell_id %notin% remove_celllist]
  }
  return(sce)
}
