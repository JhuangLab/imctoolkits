#' @title cal_distance
#' @name cal_distance using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param tocell string, CD8T.
#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- cal_distance(sce, tocell = "CD8T")
#' df <- cal_distance(sce)
#'}
cal_distance <- function(sce, tocell = "CD8T") {
  cli::cli_alert_info("Working on: {tocell}")
  coldata_dt <- colData(sce) %>% as.data.table
  cells <- setdiff(unique(coldata_dt$cell_type10), tocell)
  # min distance between immune cells and CD8T, cr: ly
  dt_to <- coldata_dt[cell_type10 %in% tocell,
                      .(sample_id, Pos_X, Pos_Y)]
  setnames(dt_to, c("Pos_X", "Pos_Y"), c("Pos_X_to", "Pos_Y_to"))
  dt_im <- coldata_dt[cell_type10 %in% cells,
                      .(sample_id, cell_id, Pos_X, Pos_Y, stype2, cell_type10)]
  n_thred <- 50000
  if(nrow(dt_to) > n_thred) {
    cli::cli_alert_info("There are {nrow(dt_to)} rows. Split it by {n_thred}.")
    split <- dt_to %>% group_split(group_id = row_number() %/% n_thred)
    pt_wide_list <- list()
    for (i in 1:length(split)) {
      dt <- merge(dt_im, split[[i]], all.x = TRUE, by = "sample_id", allow.cartesian = TRUE)
      dt[, dist := sqrt((Pos_X - Pos_X_to)^2 + (Pos_Y - Pos_Y_to)^2)]
      pt_wide_list[[i]] <- dt[, min(dist), by = .(sample_id, cell_type10, stype2)]
    }
    pt_wide <- rbindlist(pt_wide_list)
    pt_wide <- pt_wide[, min(V1, na.rm = TRUE), by = .(sample_id, cell_type10, stype2)]
    pt_wide$V1[pt_wide$V1 == Inf] <- NA
  } else {
    dt <- merge(dt_im, dt_to, all.x = TRUE, by = "sample_id", allow.cartesian = TRUE)
    dt[, dist := sqrt((Pos_X - Pos_X_to)^2 + (Pos_Y - Pos_Y_to)^2)]
    pt_wide <- dt[, min(dist, na.rm = TRUE), by = .(sample_id, cell_type10, stype2)]
    pt_wide$V1[pt_wide$V1 == Inf] <- NA
  }
  setnames(pt_wide, "V1", "distance")
  # construct wide data for comparison
  pt_wide$sample_id <- sub("\\_.*", "", pt_wide$sample_id)
  pt_wide <- pt_wide[!duplicated(pt_wide[, c("sample_id", "cell_type10")]), ]
  pt_wide <- pt_wide[, c("sample_id", "cell_type10", "distance")] %>%
    pivot_wider(names_from = sample_id, values_from = distance, values_fill = 0) %>%
    as.data.frame %>% column_to_rownames(var = "cell_type10") %>% as.data.frame()
  return(pt_wide)
}

#' @title cal_distance
#' @name cal_distance using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param tocell string, CD8T.
#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- cal_distance(sce, tocell = "CD8T")
#' a <- cal_distance_all(sce)
#'}
cal_distance_all <- function(sce) {
  dist_lst <- list()
  coldata_dt <- colData(sce) %>% as.data.frame()
  cells <- unique(coldata_dt$cell_type10)
  for(cell in cells){
    cli::cli_h1("Processing distance to: {cell}")
    dist_lst[[cell]] <- cal_distance(sce, tocell = cell)
  }
  metadata(sce)$distance <- dist_lst
  return(sce)
}

#' @title dist_min
#' @name dist_min using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param tocell string, CD8T.
#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' # usage: calculate the distance to the nearest Epithelial_normal cell for each CD8T cell
#' sce <- readRDS("../05.anno/all_anno.rds")
#' coldata <- colData(sce) %>% as.data.table
#' res <- coldata[, dist_min(.SD, "CD8T", "Epithelial_normal"), by = sample_tiff_id]
#' }
dist_min <- function(dt, q_cell, cell){
  dt <- dt |> as.data.table()
  cells <- dt$cell_type
  if((q_cell %in% cells) & (cell %in% cells)){
    data <- dt[cell_type == cell, .(Pos_X, Pos_Y)] %>% as.matrix
    query <- dt[cell_type == q_cell, .(Pos_X, Pos_Y)] %>% as.matrix
    d_min <- BiocNeighbors::queryKNN(data, query, k = 1)$distance %>% as.vector
    return(d_min)
  } else {
    return(NA_real_)
  }
}
