#' steinbock
#' conda activate steinbock
#' Wrapper for steinbock python package, steinbock need to be installed.
#'
#' @param stype character string, either c("normal", "paracancerous", "puncture", "tumor")
#' @param sample_id character string, for example: p155_pdac
#' @param panel_fn character string, path of panel file.
#' for example: panel_fn <- glue("~/projects/{project}/data/{dataset}/{species}/hyperion/panel.csv")
#' @param dmax integer, default 30
#' @param analysis_dir character string, a path where the analysis output.
#' @examples
#' \dontrun{
#' conda activate steinbock
#' steinbock(stype, sample_id, panel_fn, dmax = 30, analysis_dir)
#' We will setup run steinbock with singularity or docker.
#' }

steinbock <- function(stype, sample_id, panel_fn, dmax = 30, analysis_dir){
  steinbock_ss <- path_package("imctoolkits", "extdata/scripts", "steinbock_ss.R")
  classify_dir <- glue("{analysis_dir}/{stype}/measure/{sample_id}/classify") %>% checkdir()
  deepcell_dir <- glue("{analysis_dir}/{stype}/measure/{sample_id}/deepcell") %>% checkdir()
  intensities_dir <- glue("{analysis_dir}/{stype}/measure/{sample_id}/intensities") %>% checkdir()
  regionprops_dir <- glue("{analysis_dir}/{stype}/measure/{sample_id}/regionprops") %>% checkdir()
  neighbors_centroids_dir <- glue("{analysis_dir}/{stype}/measure/{sample_id}/neighbors") %>% checkdir()
  log_dir <- glue("logs/{stype}") %>% checkdir()
  log_fn <- glue("{log_dir}/{sample_id}.log")
  img_dir <- glue("{analysis_dir}/{stype}/img/{sample_id}")
  cmd <- glue("{steinbock_ss} --sample_id={sample_id} --classify_dir={classify_dir} --deepcell_dir={deepcell_dir} \
              --intensities_dir={intensities_dir} --regionprops_dir={regionprops_dir} --neighbors_centroids_dir={neighbors_centroids_dir} \
              --img_dir={img_dir} --dmax={dmax} --panel_fn={panel_fn} &>{log_fn} &") %>% trim_spaces()
  run_cmds(cmd)
}

steinbock_ilastik <- function(config){
  sample_id <- config$steinbock$samplename
  cropsize <- config$steinbock$cropsize
  steinbock_path <- config$steinbock$bin_path
  steinbock <- glue("{steinbock_path}/steinbock")
  panel_fn <- config$steinbock$panel_fn
  ilastik_dir <-
  classify_dir <-
  ilastik_crops <-
  ilp <-
  ilastik_probabilities <-
  cmd1 <- glue("{steinbock} classify ilastik prepare --cropsize {cropsize} --seed 123 --img {img_dir} \
              --panel {panel_fn} -o {ilastik_dir}/pixel_classifier.ilp --imgout {classify_dir} \
              --cropout {ilastik_crops} &>{ilastik_dir}/{sample_id}.log") %>% trim_spaces()
  cmd3 <- glue("{steinbock} classify ilastik run --ilp {ilp} --img {classify_dir} \
             -o {ilastik_probabilities} &>{ilastik_dir}/{sample_id}_ilastik_probabilities.log") %>% trim_spaces()
}

deepcell <- function(config){
  cmd2 <- glue("~/bin/tools/steinbock segment deepcell --minmax --img {img_dir} --panel {panel_fn} -o {deepcell_dir}") %>% trim_spaces()

}

cellprofiler <- function(config){
  cmd4 <- glue("~/bin/tools/steinbock segment cellprofiler run --pipe {pipe} --probabs {ilastik_probabilities} \
             -o {masks} &>{ilastik_dir}/{sample_id}_pipe.log") %>% trim_spaces()
}

measure <- function(config){
  cmd5 <- glue("~/bin/tools/steinbock measure intensities --img {img_dir} --masks {masks} --panel {panel_fn} -o {intensities_dir}")
  cmd6 <- glue("~/bin/tools/steinbock measure neighbors --type centroids --dmax {dmax} --masks {masks} -o {neighbors_centroids_dir}")
  cmd_graphs <- glue("steinbock export graphs --neighbors {neighbors_centroids_dir} --data {intensities_dir} --data {regionprops_dir} -o {graphs_outdir}")
  cmds <- c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd_graphs)
  cmds <- c(cmd1, cmd4)
  #run_cmds(cmds)
}



