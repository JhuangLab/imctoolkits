#' #' @title qc_hyp
#' #' @name qc_hyp using this to avoid function name confiction.
#' #' @param sce object, a singcellexperiment object.
#' #' @param sd_num integer, set the number of range standard deviation. cell aera out of this size will be removed.
#' #' @description  Filter ion makers. We think most of markers will be kept. Therefore, we
#' #'               feed a makers vector will be filtered.
#' #' @return a singcellexperiment object
#' #' @export
#' #' @examples
#' #'\dontrun{
#' #'   sce <- qc_hyp(sce, sd_num = 3)
#' #'   sce <- filter_makers(sce, makers = cells)
#' #'}
#' qc_hyp <- function(sce, sd_num = 4, iqr_num = 7){
#'   sce <- sce[, colSums(counts(sce)) > 0]
#'   message(Sys.time(), " filter cell outline: ", sd_num , " standard deviation.")
#'   area <- colData(sce)$area
#'   fil <- abs(area - mean(area)) < sd_num * sd(area)
#'   message(Sys.time(), " based on cell area, ", sum(!fil), " cells are filtered by sd: ", sd_num )
#'   sce <- sce[, fil]
#'   dat <- assay(sce, "counts") %>% t()
#'   for(marker in colnames(dat)){
#'     sce <- marker_IQR_delete(sce, marker = marker, iqr_num = iqr_num)
#'   }
#'   return(sce)
#' }


#' @title remove_duplicate
#' @name remove_duplicate using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @description  Filter ion makers. We think most of markers will be kept. Therefore, we
#'               feed a makers vector will be filtered.
#' @return a singcellexperiment object
#' @export
#'
remove_duplicate <- function(sce) {
  message(Sys.time(), ": remove_duplicate")
  suppressMessages(features_per_ROI <- as.data.frame(colData(sce)) %>%
                     group_by(sample_tiff_id, sample_id) %>%   ## All intrinsic data are tibbles
                     dplyr::summarise(                                ## to avoid potential errors.
                       amount = n(), area = mean(area),
                       eccentricity = mean(eccentricity)
                     ))                                        ## Use the metrics to find duplicate
  features_per_ROI <- features_per_ROI %>%
    group_by(sample_id) %>%
    dplyr::mutate(ROIs = n_distinct(sample_tiff_id)) ## Calculate ROI counts for each sample
  dup <- features_per_ROI %>%
    group_by(amount, area, eccentricity) %>%
    dplyr::filter(n() > 1)
  if (dim(dup)[1] > 0) {
    dup_type_1 <- dup %>%
      group_by(sample_tiff_id) %>%
      dplyr::filter(n() > 1)
    if (dim(dup_type_1)[1] > 0) {
      message(paste0(
        "WARNING: Found duplicated samples: \n",
        paste(unique(dup_type_1$sample_id), collapse = ", ")
      ))
    }
    suppressMessages(dup_type_2 <- anti_join(dup, dup_type_1))
    dup_type_1 <- dup_type_1 %>% dplyr::arrange(ROIs, sample_id) ## Decides which duplicate
    dup_type_2 <- dup_type_2 %>% dplyr::arrange(ROIs, sample_id) ## to remove by ROI counts
    remain_1 <- dup_type_1 %>%
      dplyr::distinct(amount, area, .keep_all = TRUE)
    remain_2 <- dup_type_2 %>%
      dplyr::distinct(amount, area, .keep_all = TRUE)
    remains <- rbind(remain_1, remain_2)
    remains2 <- rbind(dup_type_1, remain_2)
    # Please switch the next line of code, if type 1 duplicates are to be removed
    # suppressMessages(toberemoved <- anti_join(dup, remains2))
    suppressMessages(toberemoved <- anti_join(dup, remains))
    if (dim(toberemoved)[1] > 0) {
      message(paste0(
        "WARNING: ", dim(toberemoved)[1], " ROIs are removed"
      ))
      message(paste0(
        "Removed ROIs: \n",
        paste(unique(toberemoved$sample_tiff_id), collapse = ", ")
      ))
    }
    scematrix <- as.matrix(cbind(sce$sample_id, sce$sample_tiff_id))
    removematrix <- as.matrix(cbind(toberemoved$sample_id, toberemoved$sample_tiff_id))
    removeindex <- fastercheck(scematrix, removematrix)
    sce <- sce[, !removeindex]
    sce <- sce[, as.vector(!duplicated(t(assay(sce))))]
  }
  return(sce)
}
#' @title fastercheck
#' @name x using this to avoid function name confiction.
#' @param matrix object, a singcellexperiment object.
#' @description  not export
#' @return a singcellexperiment object
fastercheck <- function(x, matrix){
  nc <- ncol(matrix)
  rec.check <- function(r,i,id){
    id[id] <- matrix[id,i] %in% r[i]
    if(i<nc & any(id)) rec.check(r,i+1,id) else any(id)
  }
  apply(x,1,rec.check,1,rep(TRUE,nrow(matrix)))
}
#' @title filter_rare_marker_exp
#' @name filter_rare_marker_exp using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @description  Filter ion makers. We think most of markers will be kept. Therefore, we
#'               feed a makers vector will be filtered.
#' @return a singcellexperiment object
#' @export
#'
filter_rare_marker_exp <- function(sce, exps_markers_num = 5){
  message(Sys.time(), ": filter_rare_marker_exp")
  fil1 <- colSums(counts(sce)) / dim(sce)[1] < exps_markers_num / dim(sce)[1]
  fil2 <- colSums(counts(sce) > 0) < exps_markers_num
  sce <- sce[, !(fil1 & fil2)]
  dat <- t(counts(sce))
  return(sce)
}

#' @title filter_cells_by_area
#' @name filter_cells_by_area using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.
#' @description  Filter ion makers. We think most of markers will be kept. Therefore, we
#'               feed a makers vector will be filtered.
#' @return a singcellexperiment object
#' @export
#'
filter_cells_by_area <- function(sce, sd_num = 4){
  message(Sys.time(), ": filter_cells_by_area")
  sce <- sce[, colSums(counts(sce)) > 0]
  message(Sys.time(), " filter cell outline: ", sd_num , " standard deviation.")
  area <- colData(sce)$area
  fil <- abs(area - mean(area)) < sd_num * sd(area)
  message(Sys.time(), " based on cell area, ", sum(!fil), " cells are filtered by sd: ", sd_num )
  sce <- sce[, fil]
  return(sce)
}
# sce <- rds_fil
# dat <- t(counts(sce))
# head(dat[, 1:2])
# dim(dat)
# dat_dist <- dplyr::distinct(dat %>% as.data.frame())
# dim(dat_dist)

