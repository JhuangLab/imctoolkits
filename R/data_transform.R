#' @title imcAsinh_trans
#' @name imcAsinh_trans
#' @param sce object, a singcellexperiment object
#' @export
#' @examples
#' \dontrun{
#' sce <- imcAsinh_trans(sce)
#' assayNames(sce)
#'}
#'
imcAsinh_trans <- function(sce, assay = "compcounts"){
  exp <- assay(sce, assay)
  dat <- t(exp)
  assay(sce, "imcAsinh") <- apply(dat, 2, imcAsinh) %>% as.data.frame() %>% t()
  return(sce)
}

#' @title simpleAsinh_trans
#' @name simpleAsinh_trans
#' @param sce object, a singcellexperiment object
#' @export
#'
#'
simpleAsinh_trans <- function(sce, assay = "compcounts"){
  exp <- assay(sce, assay)
  dat <- t(exp)
  assay(sce, "simpleAsinh") <- apply(dat, 2, simpleAsinh) %>% as.data.frame() %>% t()
  return(sce)
}

#' @title normalize_hyp
#' @name normalize_hyp
#' @param sce object, a singcellexperiment object
#' @export
normalize_hyp <- function(sce, method = c("autoLgcl", "lognormcounts", "scater", "log2counts"), assay = "compcounts"){
  if(is.null(method)){method = "log2counts"}
  method = match.arg(method)
  exp <- assay(sce, assay)
  if(sum(colSums(exp) == 0) > 0){
    message("You are suggested to run: sce <- sce[, colSums(exp) > 0]")
    stop("There are some columns with no singling, \n
         or colSums(exp) == 0 .\n
         running this code: sce <- sce[, colSums(exp) > 0] ? ")
  }
  message(Sys.time(), " normalize data with method: ", method)
  switch(method,
         scater = {
           sce <- scater::logNormCounts(sce, transform = "log", assay.type = assay,
                                        BPPARAM = BiocParallel::SerialParam())
           # dat_tmp <- assay(sce, assay) %>% t() %>% as.matrix()
           # dat_scale <- apply(dat_tmp, 1, scale)
           # rownames(dat_scale) <- row.names(assay(sce, assay))
           # assay(sce, assay) <- dat_scale
         },
         autoLgcl = {
           sce <- lgcl_transform(sce)
         },
         lognormcounts = {
           sce <- lognormcounts(sce)
           # sce <- scuttle::logNormCounts(sce)
           # assay(sce, assay) <- assay(sce, assay)
         },
          log2counts = {
           sce <- log2counts(sce)
         })
  return(sce)
}

#' @title lgcl_transform
#' @name lgcl_transform
#' @param sce object, a singcellexperiment object
lgcl_transform <- function(sce, assay_name = "counts"){
  channels <- rownames(sce)
  exprs <- t(assay(sce, assay_name))
  ff <- new("flowFrame", exprs = exprs)
  lgcl <- flowCore::estimateLogicle(ff, channels = channels)
  after <- flowCore::transform(ff, lgcl)
  exprs_n <- t(after@exprs)
  assay(sce, glue("autoLgcl")) <- exprs_n
  return(sce)
}

#' @title lognormcounts
#' @name lognormcounts
#' @param sce object, a singcellexperiment object
lognormcounts <- function(sce, assay = "compcounts", verbose = F){
  if (verbose) message(paste0(Sys.time(), " Determining normalization factors"))
  all.log.data <- abs(t(assay(sce, assay)))
  cs <- apply(all.log.data, 2, sum)
  norm_factors <- (10**ceiling(log10(median(cs))))/cs
  norm_factors_idx <- which(!is.na(norm_factors))
  if (length(which(is.na(norm_factors))) > 0) {
    warning(paste0(Sys.time(), " Unavailable log data column, please check your data"))
  }
  if (verbose) message(paste0(Sys.time(), " Normalization and log-transformation."))
  all.log.data[, norm_factors_idx] <- round(
    log10(sweep(all.log.data[, norm_factors_idx], 2,
                norm_factors[norm_factors_idx], "*")+1), digits=3)
  assay(sce, "lognormcounts") <- t(all.log.data)
  return(sce)
}

#' @title lognormcounts_cell
#' @name lognormcounts_cell
#' @param sce object, a singcellexperiment object
lognormcounts_cell <- function(sce, assay = "compcounts"){
  message(Sys.time(), " normalize data with method: logcounts")
  counts <- assay(sce, assay)
  libsizes <- colSums(counts)
  size.factors <- libsizes/mean(libsizes)
  assay(sce, "lognormcounts_cell") <- log2(t(t(counts)/size.factors) + 1)
  return(sce)
}

#' @title logcounts
#' @name logcounts
#' @param sce object, a singcellexperiment object
log2counts <- function(sce, assay = "compcounts"){
  message(Sys.time(), " normalize data with method: logcounts")
  counts <- assay(sce, assay)
  assay(sce, "log2counts") <- t(log2(t(counts) + 1))
  # libsizes <- colSums(counts)
  # size.factors <- libsizes/mean(libsizes)
  # assay(sce, assay) <- log2(t(t(counts)/size.factors) + 1)
  return(sce)
}

#' @title range_scale
#' @name range_scale
#' @param sce object, a singcellexperiment object
range_scale <- function(sce, assay = "log2counts_censored"){
  message(Sys.time(), " normalize data with method: scales::rescale")
  mat <- assay(sce, assay) %>% t() %>% as.data.frame()
  mat$roi <- colData(sce)$sample_tiff_id
  mat_censor <- mat %>% group_by(roi) %>%
    mutate(across(where(is.numeric), ~ scales::rescale(.x) )) %>%
    ungroup() %>% dplyr::select(-roi) %>% as.data.frame()
  assy_censored <- glue("{assay}_scaled")
  rownames(mat_censor) <- rownames(mat)
  assay(sce, assy_censored) <- t(mat_censor)
  return(sce)
}


#' @name simpleAsinh
#' @title simpleAsinh
#' @param value object, a singcellexperiment object
simpleAsinh <- function(value, cofactor = 5){
  value <- value / cofactor
  value <- asinh(value)
  return(value)
}

#' @name simpleAsinh_Zero_Max
#' @title simpleAsinh_Zero_Max
#' @param value object, a singcellexperiment object
simpleAsinh_Zero_Max <- function(value, cofactor = 5, b=0) {
  value <- value / cofactor
  value <- asinh(value)
  if(!all(value==0)){
    value<-value/max(value,b)
  }
  return(value)
}

#' @title imcAsinh
#' @name imcAsinh
#' @param value object, a singcellexperiment object
imcAsinh <- function (value, cofactor = 5){
  value <- value - 0.05
  loID <- which(value <= 0)
  if (length(loID) > 0) {
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  }
  value <- value/cofactor
  value <- asinh(value)
  return(value)
}
#' @title percentile_hyp
#' @name percentile_hyp
#' @param sce object, a singcellexperiment object
#' @description Data transformation and normalization
#'   The presented data were not transformed, and all analyses were
#'   based on raw IMC measurements. Single-cell marker expressions are
#'   summarized by mean pixel values for each channel. The single-cell
#'   data were censored at the 99th percentile to remove outliers, and
#'   z-scored cluster means were visualized in heat maps. For t-SNE and
#'   PhenoGraph the data were normalized to the 99th percentile, as is
#'   suggested for these algorithms23,31. To visualize the number of cells
#'   per image or patient and for survival modelling, the cell counts were
#'   normalized by the image area (total number of pixels) and displayed
#'   as cell density. For Cox proportional hazards survival modelling, these
#'   densities were multiplied by a factor of 107 to yield values larger than 1
#'   and then log-transformed. ref: https://doi.org/10.1038/s41586-019-1876-x
#' @export

percentile_hyp <- function(sce, assay = "compcounts", percentile = 0.99){
  message(Sys.time(), ": percentile censor")
  dat <- t(assay(sce, assay))
  # dat_per <- purrr::map_df(1:ncol(dat), ~ scales::oob_squish(dat[, .x], quantile(dat[, .x], c(0, percentile)))) %>%
  #     as.data.frame()
  for(i in 1:ncol(dat)){
    dat[, i] <- scales::oob_squish(dat[, i], quantile(dat[, i], c(0, percentile)))
  }
  #rownames(dat_per) <- colnames(dat)
  assay(sce, assay) <- t(dat)
  return(sce)
}
#' @title zscore_censor
#' @name zscore_censor
#' @param sce object, a singcellexperiment object.
#' @param range zscore sd range.
#' @import data.table
#' @description Data transformation and normalization
#'   The presented data were not transformed, and all analyses were
#'   based on raw IMC measurements. Single-cell marker expressions are
#'   summarized by mean pixel values for each channel. The single-cell
#'   data were censored at the 99th percentile to remove outliers, and
#'   z-scored cluster means were visualized in heat maps. For t-SNE and
#'   PhenoGraph the data were normalized to the 99th percentile, as is
#'   suggested for these algorithms23,31. To visualize the number of cells
#'   per image or patient and for survival modelling, the cell counts were
#'   normalized by the image area (total number of pixels) and displayed
#'   as cell density. For Cox proportional hazards survival modelling, these
#'   densities were multiplied by a factor of 107 to yield values larger than 1
#'   and then log-transformed. ref: https://doi.org/10.1038/s41586-019-1876-x
#' @export
#'
zscore_censor <- function(obj, assay = "compcounts", range = c(-2.5, 2.5)) {
  message(Sys.time(), ": zscore censor")
  #mat <- counts(obj)
  mat <- assay(obj, assay)
  roi <- data.table::data.table(cell_id = colnames(obj), roi = factor(colData(obj)$sample_tiff_id))
  dt <- as.data.table(mat, keep.rownames = "channel") %>%
    data.table::melt(id.vars = "channel", variable.name = "cell_id", value.name = "intensity") %>%
    merge(roi)
  dt[, zscore := scale(intensity), by = .(channel, roi)]
  dt[is.na(zscore), zscore := 0]
  dt[, zscore_censored := scales::oob_squish(zscore, range = range)]
  dt_cast <- dt %>% dcast(channel ~ cell_id, value.var = "zscore_censored")
  mat_z <- dt_cast[, -1] %>% as.matrix
  rownames(mat_z) <- dt_cast$channel
  mat_z <- mat_z[match(rownames(obj), rownames(mat_z)), match(colnames(obj), colnames(mat_z))]
  assay(obj, "zscore") <- mat_z
  return(obj)
}

#' @title sd_censor
#' @name sd_censor
#' @param x value vector
#' @param range zscore sd range.
#' @export
sd_censor <- function(x, range){
  sd_range_down <- mean(x, na.rm = T) + sd(x, na.rm = T) * range[1]
  sd_range_up <- mean(x, na.rm = T) + sd(x, na.rm = T) * range[2]
  sd_range <- c(sd_range_down, sd_range_up)
  scales::oob_squish(x, range = sd_range)
}

#' @title sd_censor
#' @name sd_censor
#' @param x value vector
#' @param range zscore sd range.
#' @export
sd_censor_max <- function(x, range){
  sd_range_down <- min(x, na.rm = T)
  sd_range_up <- mean(x, na.rm = T) + sd(x, na.rm = T) * range[2]
  sd_range <- c(sd_range_down, sd_range_up)
  scales::oob_squish(x, range = sd_range)
}

#' @title value_censor
#' @name value_censor
#' @param sce object, a singcellexperiment object.
#' @param range zscore sd range.
#' @import data.table
#' @description Data transformation and normalization
#'   The presented data were not transformed, and all analyses were
#'   based on raw IMC measurements. Single-cell marker expressions are
#'   summarized by mean pixel values for each channel. The single-cell
#'   data were censored at the 99th percentile to remove outliers, and
#'   z-scored cluster means were visualized in heat maps. For t-SNE and
#'   PhenoGraph the data were normalized to the 99th percentile, as is
#'   suggested for these algorithms23,31. To visualize the number of cells
#'   per image or patient and for survival modelling, the cell counts were
#'   normalized by the image area (total number of pixels) and displayed
#'   as cell density. For Cox proportional hazards survival modelling, these
#'   densities were multiplied by a factor of 107 to yield values larger than 1
#'   and then log-transformed. ref: https://doi.org/10.1038/s41586-019-1876-x
#' @export
#'
value_censor <- function(obj, assay = "compcounts", range = c(-2.5, 2.5)) {
  message(Sys.time(), ": value censor on ", assay)
  mat <- assay(obj, assay) %>% t() %>% as.data.frame()
  mat$roi <- colData(obj)$sample_tiff_id
  mat_censor <- mat %>% group_by(roi) %>%
                mutate(across(where(is.numeric), ~ sd_censor_max(.x, range) )) %>%
                ungroup() %>% dplyr::select(-roi) %>% as.data.frame()
  assy_censored <- glue("{assay}_censored")
  rownames(mat_censor) <- rownames(mat)
  assay(obj, assy_censored) <- t(mat_censor)
  return(obj)
}

#' @title spillover_correction
#' @name spillover_correction
#' @param sce object, a singcellexperiment object.
#' @param co_factor zscore sd range.
#' @description Data transformation and normalization
#'   The presented data were not transformed, and all analyses were
#'   based on raw IMC measurements. Single-cell marker expressions are
#'   summarized by mean pixel values for each channel. The single-cell
#'   data were censored at the 99th percentile to remove outliers, and
#'   z-scored cluster means were visualized in heat maps. For t-SNE and
#'   PhenoGraph the data were normalized to the 99th percentile, as is
#'   suggested for these algorithms23,31. To visualize the number of cells
#'   per image or patient and for survival modelling, the cell counts were
#'   normalized by the image area (total number of pixels) and displayed
#'   as cell density. For Cox proportional hazards survival modelling, these
#'   densities were multiplied by a factor of 107 to yield values larger than 1
#'   and then log-transformed. ref: https://doi.org/10.1038/s41586-019-1876-x
#' @export
#'
spillover_correction <- function(sce, assayname = "counts", methodname = "nnls", co_factor = 5) {
  message(Sys.time(), ": spillover correction")
  fn <- fs::path_package("jhuanglabHyperion", "extdata/data/Isotope_Purity_Matrix.csv")
  compensation_matrix <- read.csv(fn)
  rownames(compensation_matrix) <- compensation_matrix[, 1]
  compensation_matrix <- compensation_matrix[, -1]
  compensation_matrix[is.na(compensation_matrix)] <- 0
  compensation_matrix <- compensation_matrix * 0.01
  stopifnot("input must be a SingleCellExperiment object" = "SingleCellExperiment" %in% class(sce))
  if (is.null(sce@rowRanges@elementMetadata@listData$channel) | is.null(sce@rowRanges@elementMetadata@listData$name)) {
    message("Please make sure [channel] and [name] are under \n sce@rowRanges@elementMetadata@listData")
  } else {
    sce@rowRanges@elementMetadata@listData$channel_name <- sce@rowRanges@elementMetadata@listData$channel
    sce@rowRanges@elementMetadata@listData$marker_name <- sce@rowRanges@elementMetadata@listData$name
    sce@rowRanges@elementMetadata@listData$marker_class <- "none"
    sce@rowRanges@elementMetadata@listData$use_channel <- sce@rowRanges@elementMetadata@listData$keep == 1
    sce <- CATALYST::compCytof(sce, as.matrix(compensation_matrix),
                     assay = assayname, method = methodname,
                     cofactor = co_factor, overwrite = FALSE
    )
  }
  #corrected intensities are in assay(sce, "compcounts")
  return(sce)
}

#' @title .normalize_quantile
#' @name .normalize_quantile
#' @param x value vector
#' @param range zscore sd range.
.normalize_quantile <- function(x, method = "qt_norm", range = c(-2.5, 2.5), quantile = 0.99){
  switch(method,
         quantile = {
           q <- quantile(x, quantile)
           x[which(x > q)] <- q
         },
         qt_norm = {
           q <- quantile(x, quantile)
           x <- x/q
           x[x > 1] <- 1
         },
         sd = { # censored by zscore and return the original value
           zscore <- scale(x)
           zscore[is.na(zscore)] <- 0
           y1 <- x[which(zscore > range[1])]
           v1 <- y1[which.min(y1)]
           y2 <- x[which(zscore < range[2])]
           v2 <- y2[which.max(y2)]
           x[which(x < v1)] <- v1
           x[which(x > v2)] <- v2
         },
         zscore = {
           zscore <- scale(x)
           zscore[is.na(zscore)] <- 0
           x <- scales::oob_squish(zscore, range = range) %>% as.vector
         }
  )
  return(x)
}

#' @title normalize_hyp
#' @name normalize_quantile
#' @param obj value vector
#' @param assay zscore sd range.
#' @export
#' @example
#' sce <- normalize_quantile(sce)
normalize_quantile <- function(obj, assay = "compcounts", ...){
  assay_new <- paste0(assay, "_censored")
  m <- assay(obj, assay)
  roi <- data.table(cell_id = colnames(obj),
                    roi = factor(colData(obj)$sample_tiff_id))
  dt <- as.data.table(m, keep.rownames = "channel") %>%
        melt(id.vars = "channel", variable.name = "cell_id", value.name = "intensity") %>%
        merge(roi)
  dt[, intensity_censor := .normalize_quantile(intensity, ...), by = .(channel, roi)]
  dt_cast <- dcast(dt, channel ~ cell_id, value.var = "intensity_censor")
  m_cs <- dt_cast[, -1] %>% as.matrix
  rownames(m_cs) <- dt_cast$channel
  m_cs <- m_cs[match(rownames(obj), rownames(m_cs)), match(colnames(obj), colnames(m_cs))]
  assay(obj, assay_new) <- m_cs
  return(obj)
}
