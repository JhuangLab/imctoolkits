#' @title pca_hyp
#' @name pca_hyp
#' @param sce object, a singcellexperiment object
#' @export
pca_hyp <- function(sce, assay = "log2counts_censored", center = TRUE, scale. = FALSE,
                    verbose = FALSE, ...){
  message(Sys.time(), " run pca analysis on sce object.")
  #pca_data <- prcomp(t(assay(sce, assay)), rank=dims)
  pca_data <- gmodels::fast.prcomp(t(assay(sce, assay)), retx = TRUE,
                                   center = center, scale. = scale., ...)
  colnames(pca_data$x) <- paste0("PC_", seq_len(ncol(pca_data$x)))
  reducedDims(sce)$pca <- pca_data$x[, 1:2]
  return(sce)
}

#' @title tsne_hyp
#' @name tsne_hyp using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object. For more details, see \code{\link[Rtsne]{Rtsne}}.
#' @import Rtsne
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' @return a singcellexperiment object
#' @references
#'    Maaten, L. Van Der, 2014. Accelerating t-SNE using Tree-Based
#'    Algorithms. Journal of Machine Learning Research, 15, p.3221-3245.
#'
#'    van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional
#'    Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
#' @export
#'
#' @examples
#'
#' sce <- normalize_hyp(sce)
#' sce <- tsne_hyp(sce)
#'

tsne_hyp <- function(sce, assay = "log2counts_censored", dims = 40){
  message(Sys.time(), " run Rtsne analysis on sce object.")
  tsne_data <- Rtsne::Rtsne(t(assay(sce, assay)), check_duplicates = FALSE, num_threads = 0)
  colnames(tsne_data$Y) <- c("tSNE_1", "tSNE_2")
  reducedDims(sce)$tsne <- tsne_data$Y
  return(sce)
}

#' @title fast_tsne
#' @name fast_tsne
#' @param sce object, a singcellexperiment object
#' @export
fast_tsne <- function(sce, assay_name = "log2counts_censored"){
  source("/cluster/apps/FIt-SNE/1.2.1/fast_tsne.R", chdir = T)
  mat <- t(assay(sce, assay_name))
  tsne_result1 <- fftRtsne(mat)
  tsne_result1 <- as.data.frame(tsne_result1)
  colnames(tsne_result1) <- c("fasttsne_1","fasttsne_2")
  rownames(tsne_result1) <- rownames(mat)
  reducedDims(sce)$fast_tsne <- tsne_result1
  return(sce)
}

#' @title umap_hyp
#' @name umap_hyp
#' @param sce object, a singcellexperiment object
#' @export
umap_hyp <- function(sce, assay = "log2counts_censored", dims = 40, umap.method = "uwot"){
  message(Sys.time(), " run umap::umap analysis on sce object.")
  cpu_num <- parallel::detectCores()
  if(umap.method == "uwot"){
    umap_data <- uwot::umap(t(assay(sce, assay)), n_threads = cpu_num)
    colnames(umap_data) <- c("UMAP_1", "UMAP_2")
    rownames(umap_data) <- colnames(sce)
    reducedDims(sce)$umap <- umap_data
  }else{
    umap_data <- umap::umap(t(assay(sce, assay)), n_threads = cpu_num)
    colnames(umap_data$layout) <- c("UMAP_1", "UMAP_2")
    rownames(umap_data) <- colnames(sce)
    reducedDims(sce)$umap <- umap_data$layout
  }
  return(sce)
}

#' @title phate_hyp
#' @name phate_hyp
#' @param sce object, a singcellexperiment object
#' @export
phate_hyp <- function(sce, assay = "log2counts_censored", dims = 40){
  message(Sys.time(), " run phateR::phate analysis on sce object.")
  dat <- t(assay(sce, assay)) %>% as.matrix()
  phater_data <- phateR::phate(dat, n.jobs = -1)
  colnames(phater_data$embedding) <- c("PHATE_1", "PHATE_2")
  reducedDims(sce)$phate <- phater_data$embedding
  return(sce)
}

#'
#' Calculate diffusion map in IMC
#' @title diffusionmap_hyp
#' @name diffusionmap_hyp
#'
#' @param object a SCE object
#' @param sigma.use numeric. Diffusion scale parameter of the Gaussian kernel.
#'     One of '\code{local}',
#'     '\code{global}', a \code{\link[base]{numeric}} global sigma or a Sigmas object.
#'     When choosing '\code{global}', a global sigma will be calculated using find_sigmas
#'     (See \code{destiny}). A larger sigma might be necessary if the eigenvalues can not
#'    be found because of a singularity in the matrix. See \code{destiny}.
#' @param distance Distance measurement method applied to data or a distance matrix/dist.
#'    For the allowed values, see \code{destiny}
#' @param k numeric. By default is 30. \code{destiny} can be used to specify k.
#' @param density.norm logical. If TRUE, use density normalisation. See \code{destiny}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... options to pass on to the \code{destiny}.
#'
#' @seealso \code{destiny}
#'
#' @export
#' @return A SCE object
#'
#' @examples
#'\dontrun{
#' IMC.file <- system.file("extdata/IMC.rds", package = "IMCoTree")
#' IMC <- readRDS(file = IMC.file)
#'
#' IMC <- diffusionmap_hyp(IMC, verbose = TRUE)
#'}
#'
diffusionmap_hyp <- function(sce, assay = "log2counts_censored", sigma.use = NULL,
                            distance = c("euclidean", "cosine", "rankcor"),
                            k = 30,
                            density.norm = TRUE,  verbose = T,
                            ...) {
  #dm.data <- as.matrix(sce@log.data[which(sce@meta.data$dowsample == 1), sce@markers.idx])
  dm.data <- as.matrix(t(assay(sce, assay)))
  if (verbose) message(Sys.time(), " Calculating Diffusion Map.")
  # Figure out sigma
  # this function refered to URD calcDM function.
  if (is.null(sigma.use)) {
    sigma.use <- destiny::find_sigmas(dm.data, verbose = FALSE)@optimal_sigma
    if (verbose) message(Sys.time(), " Destiny determined an optimal global sigma: ", round(sigma.use, digits=3))
  } else if (is.numeric(sigma.use)) {
    if (verbose) message(Sys.time(), " Using provided global sigma: ", round(sigma.use, digits=3))
  } else if (sigma.use == "local") {
    if (verbose) message(Sys.time(), " Using local sigma ")
  } else {
    sigma.use <- destiny::find_sigmas(dm.data, verbose = FALSE)@optimal_sigma
    warning(Sys.time(), " Invalid sigma value. Using an optimal global sigma instead.")
  }
  # Calculate the Diffusion Map
  distance <- match.arg(distance)
  dm.obj <- destiny::DiffusionMap(dm.data, sigma=sigma.use, k=k, density_norm = density.norm, distance=distance, ...)

  #rownames(dm.obj@eigenvectors) <- rownames(dm.data)
  colnames(dm.obj@eigenvectors) <- paste0("DC_", seq_len(ncol(dm.obj@eigenvectors)))
  # rownames(dm.obj@transitions) <- rownames(dm.data)
  # colnames(dm.obj@transitions) <- rownames(dm.data)

  reducedDims(sce)$dfm_eign <- dm.obj@eigenvectors[ ,1:2]
  reducedDims(sce)$dfm_trans <- dm.obj@transitions
  if (verbose) message(Sys.time(), " Calculating Diffusion Map completed")
  return(sce)
}

