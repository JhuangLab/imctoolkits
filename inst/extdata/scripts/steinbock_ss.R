#!/usr/bin/env Rscript
# conda activate steinbock
pkgs <- c("fs", "tidyverse", "jhtools","futile.logger", "configr", "stringr",
          "optparse", "glue")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-s","--sample_id"),
              help="sample_id"),
  make_option(c("-c","--classify_dir"),
              help="classify_dir"),
  make_option(c("-d","--deepcell_dir"),
              help="deepcell_dir"),
  make_option(c("-i","--intensities_dir"),
              help="intensities_dir"),
  make_option(c("-r","--regionprops_dir"),
              help="regionprops_dir"),
  make_option(c("-n","--neighbors_centroids_dir"),
              help="neighbors_centroids_dir"),
  make_option(c("-m","--img_dir"),
              help="img_dir"),
  make_option(c("-p","--panel_fn"),
              help="panel_fn"),
  make_option(c("-x","--dmax"),
              help="dmax")
)
opt <- parse_args(OptionParser(option_list=option_list))
sample_id <- opt$sample_id
classify_dir <- opt$classify_dir
deepcell_dir <- opt$deepcell_dir
intensities_dir <- opt$intensities_dir
regionprops_dir <- opt$regionprops_dir
neighbors_centroids_dir <- opt$neighbors_centroids_dir
img_dir <- opt$img_dir
dmax <- opt$dmax
panel_fn <- opt$panel_fn
graphs_outdir <- stringr::str_replace_all(intensities_dir, "intensities", "graphs") %>% checkdir()
if(is.na(sample_id)){quit("no")}
steinbockBinPath <- "~/bin/immu/hyperion/steinbock/steinbock_0.10.2.sif"
steinbock <- glue("singularity run -B /cluster:/cluster {steinbockBinPath}")
cmd1 <- glue("{steinbock} classify ilastik prepare --cropsize 50 --seed 123 --img {img_dir} \
              --panel {panel_fn} -o {classify_dir}/pixel_classifier.ilp --imgout {classify_dir} \
              --cropout {classify_dir}/ilastik_crops &>{classify_dir}/{sample_id}.log") %>% trim_spaces()
cmd2 <- glue("{steinbock} segment deepcell --minmax --img {img_dir} --panel {panel_fn} -o {deepcell_dir}") %>% trim_spaces()
cmd3 <- glue("{steinbock} measure regionprops --img {img_dir} --mask {deepcell_dir} -o {regionprops_dir}")
cmd4 <- glue("{steinbock} measure intensities --img {img_dir} --mask {deepcell_dir} --panel {panel_fn} -o {intensities_dir}")
cmd5 <- glue("{steinbock} measure neighbors --type centroids --dmax {dmax} --masks {deepcell_dir} -o {neighbors_centroids_dir}")
cmd_graphs <- glue("{steinbock} export graphs --neighbors {neighbors_centroids_dir} --data {intensities_dir} --data {regionprops_dir} -o {graphs_outdir}")
cmds <- c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd_graphs)
run_cmds(cmds)
