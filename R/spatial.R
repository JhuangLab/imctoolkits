#' @title cal_closest_dist
#' @name cal_closest_dist
#' @param roi object, a singcellexperiment object
#' @param cld object, a singcellexperiment object
#' @param query_cells object, a singcellexperiment object
#' @param taget_cells object, a singcellexperiment object
#' @param k object, a singcellexperiment object
#' @export
#'
cal_closest_dist <- function(roi, cld, query_cells, taget_cells, k = 1){
  query <- cld %>% dplyr::filter(cell_type10 %in% query_cells) %>%
    dplyr::filter(sample_tiff_id == roi) %>% dplyr::select("Pos_X", "Pos_Y", "cell_id") %>%
    remove_rownames() %>% column_to_rownames(var = "cell_id") %>% as.matrix()
  data <- cld %>% dplyr::filter(cell_type10 %in% taget_cells) %>%
    dplyr::filter(sample_tiff_id == roi) %>% dplyr::select("Pos_X", "Pos_Y", "cell_id") %>%
    remove_rownames() %>% column_to_rownames(var = "cell_id") %>% as.matrix()
  if (nrow(query) != 0 & nrow(data) != 0) {
    d_min <- BiocNeighbors::queryKNN(data, query, k = k)
    res <- data.frame(from_cell = rownames(query),
                      to_cell = rownames(data)[as.vector(d_min$index)],
                      distance = as.vector(d_min$distance))
    return(res)
  }
}
