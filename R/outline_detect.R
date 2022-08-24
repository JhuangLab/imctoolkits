#' @title marker_outlier_delete
#' @name marker_outlier_delete remove outlier.
#' @param sce single cell experiment object.
#' @description
#'     http://r-statistics.co/Outlier-Treatment-With-R.html
#'     https://statsandr.com/blog/outliers-detection-in-r/
#'     https://rstudio-pubs-static.s3.amazonaws.com/300795_d5c82be1d7c548d881b78216e0b51c5d.html
#'     https://stackoverflow.com/questions/41462073/multivariate-outlier-detection-using-r-with-probability
#'     https://stackoverflow.com/questions/45289225/removing-multivariate-outliers-with-mvoutlier
#' @export
#'
marker_outlier_delete <- function(sce, marker = "CD45", k = 20){
  dat <- assay(sce, "counts") %>% t()
  marker_test <- EnvStats::rosnerTest(dat[,marker], k = k)
  fil <- marker_test$all.stats$Outlier
  index <- marker_test$all.stats$Obs.Num[fil]
  msg <- glue("{Sys.time()} based on {marker}: total {length(index)} cells are removed.")
  message(msg)
  sce <- sce[, - index]
  return(sce)
}

#' @title marker_IQR_delete
#' @name marker_IQR_delete remove outlier.
#' @param sce single cell experiment object.
#' @export

marker_IQR_delete <- function(sce, marker = "CD45", iqr_num = 5 ){
  dat <- assay(sce, "counts") %>% t()
  dat[,marker] <- log(dat[,marker]+1, 2)
  Q1 <- quantile(dat[,marker], .25)
  Q3 <- quantile(dat[,marker], .75)
  IQR <- IQR(dat[,marker])
  fil <- dat[ ,marker] > (Q1 - iqr_num * IQR) & dat[ ,marker] < (Q3 + iqr_num * IQR)
  msg <- glue("{Sys.time()} based on {marker}: total {sum(!fil)} cells are removed.")
  message(msg)
  sce <- sce[, fil]
  return(sce)
}


