#' @title sub_cluster
#' @name sub_cluster
#' @param sce
#' @param fastpg
#' @param cluster_ids
#' @param assay
#' @param k
#' @export
sub_cluster <- function(sce, fastpg, cluster_ids, assay = "logcounts_censored", k = 20){
  sce_sub <- subset_by_cluster_id(sce, fastpg, cluster_ids)
  sce_sub <- run_cluster(sce_sub, cluster.method = "fastpg", assay = assay)
  return(sce_sub)
}
