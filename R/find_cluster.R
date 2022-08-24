#'
#' Specific Clustering Method Toolkits
#' Adopt these methods from sceoTree.
#'
#' @name run_cluster
#'
#' @description Compute a specific clustering using the combined flow
#'    sceometry data. "som" \code{\link[FlowSOM]{SOM}}, "hclust" \code{\link[stats]{hclust}},
#'    "clara" \code{\link[cluster]{clara}}, "phenograph", "kmeans" \code{\link[stats]{kmeans}} are
#'    provided.
#'
#' @param object a sce object
#' @param cluster.method character. Four clustering method are provided: som, clara, kmeans and phenograph.
#'    Clustering method "hclust" and "mclust" are not recommended because of long computing time.
#' @param verbose logic. Whether to print calculation progress.
#' @param ... options to pass on to the clustering functions.
#'
#' @seealso \code{\link[FlowSOM]{SOM}}, \code{\link[stats]{hclust}},
#'    \code{\link[cluster]{clara}}, \code{\link[stats]{kmeans}}.
#'    You can use \code{run_SOM}, \code{run_clara},
#'    \code{runPhenotype}, \code{run_kmeans}, \code{run_mclust} and
#'    \code{run_hclust} to run clustering respectively.
#'
#' @export
#' @return A sce object with cluster
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' # After building an sce object
#' # Set random seed to make results reproducible
#'
#' set.seed(1)
#' sce <- run_cluster(sce, cluster.method = "som", xdim = 3, ydim = 3, verbose = TRUE)
#'
#' # K-means clustering
#' sce <- run_cluster(sce, cluster.method = "kmeans", k = 9, verbose = TRUE)
#'
#' # Clara clustering
#' sce <- run_cluster(sce, cluster.method = "clara", k = 9, verbose = TRUE)
#'
#' # phenoGraph clustering
#' sce <- run_cluster(sce, cluster.method = "phenograph", verbose = TRUE)
#'
#' # hclust clustering
#' # not recommended for large cell size
#' sce <- run_cluster(sce, cluster.method = "hclust", k = 9, verbose = TRUE)
#'
#' # mclust clustering
#' # not recommended for large cell size
#' sce <- run_cluster(sce, cluster.method = "mclust", verbose = TRUE)
#'
#'
run_cluster <- function(object, assay = "log2counts_censored", cluster.method = c("fastpg" ,"som", "kmeans", "clara", "phenograph", "hclust", "mclust"),
                       verbose = FALSE, ...) {

  if (missing(object)) {
    stop(Sys.time(), " sce object is missing ")
  }
  cluster.method <- match.arg(cluster.method)
  if (cluster.method == "som") {
    object <- run_SOM(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$som.id
  } else if (cluster.method == "hclust") {
    object <- run_hclust(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$hclust.id
  } else if (cluster.method == "mclust") {
    object <- run_mclust(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$mclust.id
  } else if (cluster.method == "clara") {
    object <- run_clara(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$clara.id
  } else if (cluster.method == "kmeans") {
    object <- run_kmeans(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$kmeans.id
  } else if (cluster.method == "phenograph") {
    object <- run_phenograph(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$phenograph.id
  } else if (cluster.method == "fastpg") {
    object <- run_fastpg(object, assay = assay, verbose = verbose, ...)
    colData(object)@listData$cluster.id <- colData(object)@listData$phenograph.id
  } else {
    warning(Sys.time(), " Invalid cluster.method parameter ")
  }

  # Initialization for root cells
  object@colData@metadata$network <- list()
  object@colData@metadata$network$is.root.cells <- 0
  object@colData@metadata$network$is.leaf.cells <- 0
  return(object)
}

#'
#' processingCluster
#'
#' @name processingCluster
#'
#' @description Calculate Principal Components Analysis (PCA), t-Distributed
#'    Stochastic Neighbor Embedding (tSNE), Diffusion Map and Uniform Manifold
#'    Approximation and Projection (UMAP) of clusters calculated by run_cluster.
#'
#' @param object a sce object
#' @param perplexity numeric. Perplexity parameter (should not be bigger than 3 *
#'    perplexity < nrow(X) - 1, see details for interpretation). See \code{\link[Rtsne]{Rtsne}}
#'    for more information.
#' @param k numeric. The parameter k in k-Nearest Neighbor.
#' @param downsampling.size numeric. Percentage of sample size of downsampling.
#'    This parameter is from 0 to 1. by default is 1.
#' @param force.resample logical. Whether to do resample if downsampling.size < 1
#' @param random.cluster logical. Whether to perfrom random downsampling. If FALSE,
#'    an uniform downsampling will be processed.
#' @param umap.config object of class umap.config. See \code{\link[umap]{umap}}.
#' @param verbose logic. Whether to print calculation progress.
#' @param ... options to pass on to the dimensionality reduction functions.
#'
#' @seealso \code{\link[umap]{umap}}, \code{\link[gmodels]{fast.prcomp}},
#'    \code{\link[Rtsne]{Rtsne}}, \code{destiny}
#'
#' @importFrom stats cutree
#'
#' @export
#' @return A sce object with dimensionality reduction of clusters
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' # After running clustering
#' set.seed(1)
#' sce <- run_cluster(sce, cluster.method = "som", xdim = 3, ydim = 3, verbose = TRUE)
#'
#' # Do not perfrom downsampling
#' sce <- processingCluster(sce, perplexity = 2)
#'
#' # Perform cluster based downsampling
#' # Only keep 50% cells
#' sce <- processingCluster(sce, perplexity = 2, downsampling.size = 0.5)
#'
#' # Processing clusters without downsampling step
#' sce <- processingCluster(sce, perplexity = 2, force.resample = FALSE)
#'
#'
#'
processingCluster <- function(object, assay = "log2counts_censored", perplexity = 5, k = 5,
                              downsampling.size = 1,
                              force.resample = TRUE,
                              random.cluster = FALSE,
                              umap.config = umap.defaults, verbose = FALSE,
                              ...) {

  if (missing(object)) {
    stop(Sys.time(), " sce object is missing ")
  }

  if (!"cluster.id" %in% colnames(colData(object))) {
    stop(Sys.time(), " cluster.id is not in colnames of sce object, please run run_cluster first ")
  }

  # checking index of markers in cluster
  cluster.meta <- fetchClustMeta(object, verbose = FALSE)
  cluster.mat <- cluster.meta[, match(object@markers, colnames(cluster.meta))]

  # run PCA
  if (verbose) message(Sys.time(), " Calculating PCA")
  pca.info <- fast.prcomp( t(cluster.mat), ...)
  colnames(pca.info$rotation) <- paste0("PC_", seq_len(ncol(pca.info$rotation)))
  if (verbose) message(Sys.time(), " Calculating tSNE")
  tsne.info <- Rtsne(as.matrix(cluster.mat), perplexity = perplexity, ...)
  colnames(tsne.info$Y) <- paste0("tSNE_", seq_len(ncol(tsne.info$Y)))
  if (verbose) message(Sys.time(), " Calculating Diffusion Map")
  dm.info <- DiffusionMap(cluster.mat, k=5, ...)
  colnames(dm.info@eigenvectors) <- paste0("DC_", seq_len(ncol(dm.info@eigenvectors)))
  if (verbose) message(Sys.time(), " Calculating UMAP")
  umap.config$n_neighbors <- k
  umap.info <- umap(cluster.mat, config = umap.config, ...)
  colnames(umap.info$layout) <- paste0("UMAP_", seq_len(ncol(umap.info$layout)))

  object@cluster <- data.frame(pca.info$rotation, tsne.info$Y, dm.info@eigenvectors, umap.info$layout)
  rownames(object@cluster) <- rownames(object@tree.meta$cluster)

  if (force.resample) {
    # Initialization
    object@network <- list()
    colData(object)$is.root.cells <- 0
    colData(object)$is.leaf.cells <- 0
    colData(object)$dowsample <- 0
    colData(object)$pseudotime <- 0
    colData(object)$traj.value <- 0
    colData(object)$traj.value.log <- 0
    colData(object)$is.root.cells <- 0
    colData(object)$is.leaf.cells <- 0
    colData(object)$branch.id <- 0
    object@pca.sdev <- vector()
    object@umap.value <- object@tsne.value <- object@pca.scores <- object@pca.value <- matrix()
    object@dm <- new("DiffusionMap")

    cell.sub <- NULL
    if (downsampling.size >= 1) {
      if (verbose) message(Sys.time(), " No downsampling performed")
      cell.name <- colData(object)$cell
    } else if ( downsampling.size <= 0) {
      warning(Sys.time(), " The value of downsampling.size must be larger than 0 ")
      cell.name <- colData(object)$cell
    } else {
      if (random.cluster) {
        cell.name <- sapply(unique(colData(object)$cluster.id), function(x) sample(colData(object)$cell[which(colData(object)$cluster.id == x)], ceiling(length(which(colData(object)$cluster.id == x)) * downsampling.size )) )
        cell.name <- unlist(cell.name)
      } else {
        cell.name <- sapply(unique(colData(object)$cluster.id), function(x) {
          cell.sub <- as.character(colData(object)$cell[which(colData(object)$cluster.id == x)])
          cell.sub <- cell.sub[seq(1, length(cell.sub), by = 1/downsampling.size)]
        } )
        cell.name <- unlist(cell.name)
      }

    }

    object@cell.name <- as.character(cell.name)
    colData(object)$dowsample[match(cell.name, colData(object)$cell)] <- 1

  }

  return(object)

}


#'
#' run_hclust
#'
#' @name run_hclust
#'
#' @description Hierarchical cluster analysis on a set of dissimilarities
#'    and methods for analyzing it.
#'
#' @param object a sce object
#' @param hclust.method character or a function. The agglomeration method to be used.
#'    This should be one of "ward.D", "ward.D2", "single", "complete", "average",
#'    "mcquitty", "median" or "centroid". Or you can specify an equation as input, for example
#'    \code{function(x) hclust(x,method = 'ward.D2')}.
#' @param dist.method character or a function. The distance measure to be used.
#'    This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary"
#'    or "minkowski". Or you can specify an equation as input, for example
#'    \code{function(x) as.dist((1-cor(t(x)))/2)}.
#' @param k numeric. The number of clusters.
#' @param verbose logical. Whether to print calculation progress.
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{dist}}
#'
#' @importFrom stats hclust dist
#'
#' @export
#' @return A sce object with cluster
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_hclust(sce, k = 9, verbose = TRUE)
#'
#'
#'
run_hclust <- function(object, assay = "log2counts_censored", k = 25,
                      hclust.method = "complete", dist.method = "euclidean",
                      verbose = FALSE) {

  if (verbose) message(Sys.time(), " Calculating Hclust.")

  # check dist parameters
  if (is.character(dist.method)) {
    d <- stats::dist(t(assay(sce, assay)), method = dist.method)
  } else if (is.function(dist.method)) {
    d <- dist.method(t(assay(sce, assay)))
  } else {
    warning(Sys.time(), " Invalid dist.method parameter.")
    d <- stats::dist(t(assay(sce, assay)))
  }

  # check hclust parameters
  if (is.character(hclust.method)) {
    hc <- stats::hclust(d, method = hclust.method)
  } else if (is.function()) {
    hc <- dist.method(d)
  } else {
    warning(Sys.time(), " Invalid hclust.method parameter.")
    hc <- stats::hclust(d)
  }


  hc.tree <- cutree(hc, k = k)

  colData(object)$hclust.id <- colData(object)$cluster.id <- hc.tree

  if (verbose) message(Sys.time(), " Calculating Hclust completed.")
  return(object)
}



#'
#' run_kmeans
#'
#' @name run_kmeans
#'
#' @description Perform k-means clustering on a data matrix.
#'
#' @param object  a sce object
#' @param k numeric. The number of clusters.
#' @param iter.max numeric. The maximum number of iterations allowed.
#' @param nstart numeric. If k is a number, how many random sets should be chosen.
#' @param algorithm character. Type of algorithm that will be choosen to calculate
#'    kmeans. Four algoritms are provided: Hartigan-Wong, Lloyd, Forgy, MacQueen.
#' @param trace logical or integer number.
#' @param scale logical. Whether to use scaled data in kmeans.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[stats]{kmeans}} function
#'
#' @return a sce object with kmeans.id in meta.data
#'
#' @seealso \code{\link[stats]{kmeans}}
#'
#' @importFrom stats kmeans
#' @export
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_kmeans(sce, k = 25, verbose = TRUE)
#'
#'
run_kmeans <- function(object, assay = "log2counts_censored", k = 25, iter.max = 10, nstart = 1,
                      algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                      trace=FALSE, scale = FALSE, verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " Calculating Kmeans.")

  if (scale) kmeans.data <- scale(t(assay(sce, assay))) else kmeans.data = t(assay(sce, assay))

  algorithm <- match.arg(algorithm)
  kmeans.info <- kmeans(kmeans.data, centers = k, iter.max = iter.max, nstart = nstart,
                        algorithm = algorithm, trace = FALSE)

  colData(object)$kmeans.id <- colData(object)$cluster.id  <- kmeans.info$cluster

  if (verbose) message(Sys.time(), " Calculating Kmeans completed.")
  return(object)
}


#'
#' run_clara
#'
#' @name run_clara
#'
#' @description Clustering a data matrix into k clusters
#'
#' @param object  a sce object
#' @param k numeric. The number of clusters. It is required that
#'    0 < k < n where n is the number of observations (i.e., n = nrow(x)).
#' @param metric character. string specifying the metric to be used for
#'    calculating dissimilarities between observations.
#' @param stand logical. Indicating if the measurements in x are
#'    standardized before calculating the dissimilarities.
#' @param samples numeric. Say N, the number of samples to be drawn from the dataset.
#'    The default is N = 5,
#' @param trace numberic. Indicating a trace level for diagnostic output during the algorithm
#' @param scale logical. Whether to use scaled data in kmeans.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[cluster]{clara}} function
#'
#' @return a sce object with clara.id in meta.data
#'
#' @seealso \code{\link[cluster]{clara}}
#'
#' @importFrom cluster clara
#' @export
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_clara(sce, k = 25, verbose = TRUE)
#'
#'
run_clara <- function(object, assay = "log2counts_censored", k = 25, metric = c("euclidean", "manhattan", "jaccard"),
                     stand = FALSE, samples = 5, scale = TRUE,
                     trace = 0, verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " Calculating Clara")

  if (scale) clara.data <- scale(t(assay(sce, assay))) else clara.data = t(assay(sce, assay))

  metric <- match.arg(metric)
  clara.info <- clara(clara.data, k = k, metric = metric, stand = stand, samples = samples,
                      trace = trace, ...)

  colData(object)$clara.id <- colData(object)$cluster.id  <- clara.info$clustering

  if (verbose) message(Sys.time(), " Calculating Clara completed.")
  return(object)
}

#'
#' run_mclust
#'
#' @name run_mclust
#'
#' @description Model-based clustering based on parameterized finite Gaussian mixture models.
#'    This function is based on \code{\link[mclust]{Mclust}}.
#'
#' @param object  a sce object
#' @param scale logical. Whether to use scaled data in Mclust.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[mclust]{Mclust}} function
#'
#' @return a sce object with mclust.id in meta.data
#'
#' @seealso \code{\link[mclust]{Mclust}}
#'
#' @export
#'
#' @importFrom mclust Mclust mclustBIC
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_mclust(sce, verbose = TRUE)
#'
#'
run_mclust <- function(object, assay = "log2counts_censored", scale = FALSE,
                      verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " Calculating Mclust.")

  if (scale) mclust.data <- scale(t(assay(sce, assay))) else mclust.data = t(assay(sce, assay))

  mod <- Mclust(mclust.data, ...)

  colData(object)$mclust.id <- colData(object)$cluster.id  <- mod$classification

  if (verbose) message(Sys.time(), " Calculating Mclust completed.")
  return(object)
}



#'
#' calculation SOM in sce object
#'
#' @description Build a self-organizing map
#'
#' @param object  a sce object
#' @param xdim  Width of the grid.
#' @param ydim  Hight of the grid.
#' @param rlen  Number of times to loop over the training data for each MST
#' @param mst   Number of times to build an MST
#' @param alpha Start and end learning rate
#' @param radius Start and end radius
#' @param init  Initialize cluster centers in a non-random way
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param method the distance measure to be used. This must be one of "euclidean",
#'      "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'      Any unambiguous substring can be given. See \code{\link[stats]{dist}}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[FlowSOM]{SOM}} function
#'
#' @return a sce object with som.id in sce object
#' @seealso \code{\link{BuildSOM}}
#'
#' @references This code is strongly based on the \code{\link[FlowSOM]{SOM}} function.
#'             Which is developed by Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2018).
#'
#' @importFrom FlowSOM SOM
#'
#' @seealso \code{\link[FlowSOM]{SOM}}
#'
#' @export
#'
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_SOM(sce, xdim = 10, ydim = 10, verbose = TRUE)
#'
#'
run_SOM <- function(object, assay = "log2counts_censored", xdim = 6, ydim = 6, rlen = 8, mst = 1,
                   alpha = c(0.05,  0.01), radius = 1, init = FALSE,
                   distf = 2, codes = NULL, importance = NULL,
                   method = "euclidean", verbose= FALSE, ...) {

  if (verbose) message(Sys.time(), " Calculating FlowSOM.")
  # FlowSOM
  flowset <- as.matrix(t(assay(object, assay)))
  suppressMessages(
    flowsom <- FlowSOM::SOM(flowset,
                            xdim = xdim, ydim = ydim, rlen = rlen, mst = mst,
                            alpha = alpha[1], radius = radius,
                            init = init,
                            distf = distf, silent = verbose,
                            codes = codes, importance = importance,
                            ...))

  # generating som network
  colData(object)$som.id <- colData(object)$cluster.id  <- flowsom$mapping[, 1]
  colData(object)$som.value <- flowsom$mapping[, 2]
  #object@metadata$som <- flowsom
  #object@som.network <- buildSOMnet(flowsom, object, method = method)
  if (verbose) message(Sys.time(), " Calculating FlowSOM completed.")
  return(object)
}

#' RphenoGraph clustering
#'
#' @description
#'    A simple R implementation of the phenograph
#'    [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm,
#'    which is a clustering method designed for high-dimensional single-cell
#'    data analysis. It works by creating a graph ("network") representing
#'    phenotypic similarities between cells by calculating the Jaccard
#'    coefficient between nearest-neighbor sets, and then identifying communities
#'    using the well known [Louvain method](https://sites.google.com/site/findcommunities/)
#'    in this graph.
#'
#' @param object a sce object.
#' @param scale logical. Whether to scale the expression matrix
#' @param knn numeric. Number of nearest neighbours, default is 30.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{igraph} function
#'
#'
#' @importFrom igraph graph.adjacency simplify distances
#' @return A sce object with cluster
#'
#' @export
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_phenograph(sce, knn = 30, verbose = TRUE)
#'
#'
run_phenograph <- function(object, assay = "log2counts_censored", k = 30, scale = FALSE, verbose = FALSE, ...){

  if (verbose) message(Sys.time(), " Calculating phenoGraph")
  if (scale) phenograph.data <- scale(t(assay(object, assay))) else phenograph.data = t(assay(object, assay))
  #mod <- Rphenograph(phenograph.data, k = 30)
  mod <- cytofkit::Rphenograph(data = phenograph.data, k = k)
  colData(object)@listData$phenograph.id <- colData(object)@listData$cluster.id  <- as.numeric(igraph::membership(mod))
  if (verbose) message(Sys.time(), " Calculating phenoGraph completed.")
  return(object)
}


#' run_fastpg clustering
#'
#' @description
#'    A simple R implementation of the FastPG
#'    [FastPG](https://github.com/sararselitsky/FastPG) algorithm
#'
#' @param object a sce object.
#' @param k numeric. Number of nearest neighbours, default is 30.
#' @param ... Parameters passing to \code{igraph} function
#'
#' @return A sce object with cluster
#'
#' @export
#' @examples
#'
#' sce.file <- system.file("extdata/sce.rds", package = "sceoTree")
#' sce <- readRDS(file = sce.file)
#'
#' sce <- run_phenograph(sce, knn = 30, verbose = TRUE)
#'
#'
run_fastpg <- function(object, assay = "log2counts_censored", k = 30, verbose = FALSE, ...){
  if (verbose) message(Sys.time(), " Calculating FastPG")
  dat <- assay(object, assay)
  clusters <- FastPG::fastCluster(t(dat), k = k, num_threads = 90)
  colData(object)@listData$fastpg.id <- clusters$communities
  if (verbose) message(Sys.time(), " Calculating FastPG completed.")
  return(object)
}


#' @name meta_clustering
#' @title meta_clustering
#' @param object singlecellexperiment object;
#' @param assay string, posible value are: "counts" "imcAsinh" "simpleAsinh" "log2counts"
#' @param secclustername string character.
#' @param metaClustering_method string
#' @param k_value integer
#' @param elbow_test bool
#' @param xdim integer, this is for FlowSOM::SOM
#' @param ydim integer, this is for FlowSOM::SOM
#' @param seed integer, set seed.
#' @description  Uses a kd-tree to find the p number of
#'    This method is from xinlei cheng's package.
#' @export
#'
meta_clustering <- function(object,
                            assay = "log2counts_censored",
                            secclustername = "FlowSOM",
                            metaClustering_method = "metaClustering_PhenoGraph",
                            k_value=10,
                            elbow_test=F,
                            xdim=40,
                            ydim=40,
                            seed=123){
  dat_logcounts <- t(assay(object, assay))
  som_map <- FlowSOM::SOM(data= as.matrix(dat_logcounts),
                       xdim=xdim,
                       ydim=ydim,
                       silent = F)
  FlowSOM_combined <- data.frame(FlowSOM = som_map$mapping[,1], dat_logcounts)
  metacluster_result <- metaClustering(FlowSOM_combined,
                                        clustername = "FlowSOM",
                                        metaClustering_method = "metaClustering_PhenoGraph",
                                        k_value=10,
                                        elbow_test=F,
                                        seed=123)
  colData(object)$cluster.id <- colData(object)$metacluster.id <- metacluster_result
  return(object)
}





