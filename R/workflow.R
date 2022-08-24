#' @title make_rds
#' @name make_rds
#' @param assay a string c("logcounts_censored", "logcounts", "zscore")
#' @param outputdir a string
#' @export

workflow <- function(sce, assay = "logcounts_censored", outputdir="~/projects/hyperion/analysis/qzhang/human/steinbock/rds/merge"){
  rds_fn <- glue("~/projects/hyperion/analysis/qzhang/human/steinbock/rds/all_tsne_umap_{assay}.rds")
  if(!fs::file_exists(rds_fn)){
    sce_anno <- sce %>% tsne_hyp(assay = assay) %>% umap_hyp(assay = assay)
    readr::write_rds(sce_anno, file = rds_fn)
  }else{
    message(system.time(), "file is existed: ", rds_fn)
  }
  rds_fn <- glue("~/projects/hyperion/analysis/qzhang/human/steinbock/rds/all_tsne_umap_phenograph_{assay}.rds")
  if(!fs::file_exists(rds_fn)){
    sce_anno <- sce_anno %>% run_cluster(cluster.method = "phenograph", assay = assay)
    readr::write_rds(sce_anno, file = rds_fn)
  }else{
    message(system.time(), "file is existed: ", rds_fn)
  }
  rds_fn <- glue("~/projects/hyperion/analysis/qzhang/human/steinbock/rds/all_all_{assay}.rds")
  if(!fs::file_exists(rds_fn)){
    sce_anno <- sce_anno %>% phate_hyp(assay = assay) %>% diffusionmap_hyp(assay = assay) %>% run_cluster(cluster.method = "som")
    readr::write_rds(sce_anno, file = rds_fn)
  }else{
    message(system.time(), "file is existed: ", rds_fn)
  }
}

#' @title update_sce
#' @name update_sce
#' @param sce a string c("logcounts_censored", "logcounts", "zscore")
#' @export
update_sce <- function(sce){
  cli::cli_alert_info("merge_cell_type")
  sce <- merge_cell_type(sce)
  cli::cli_alert_info("entropy_calc")
  sce <- entropy_calc(sce)
  cli::cli_alert_info("add_cell_ratio")
  sce <- add_cell_ratio(sce)
  cli::cli_alert_info("add_cell_per_mm")
  sce <- add_cell_per_mm(sce)
  return(sce)
}
