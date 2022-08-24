#' @title filter_makers
#' @name filter_makers using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @param makers string vector, makers will be removed.
#' @description  Filter ion makers. We think most of markers will be kept. Therefore, we
#'               feed a makers vector will be filtered.
#' @return a singcellexperiment object
#' @export
#'
#' @examples
#'\dontrun{
#'   sce <- filter_makers(sce, makers = 2000)
#'   sce <- filter_makers(sce, makers = cells)
#'}
filter_makers <- function(sce, makers = c("80ArAr", "120Sn")){
  message(Sys.time(), " filter markers: ", glue::glue_collapse(makers, sep=" "))
  fil <- rownames(sce) %in% makers
  sce <- sce[!fil, ]
  return(sce)
}
