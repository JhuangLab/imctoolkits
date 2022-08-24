#' @title make_rds
#' @name make_rds
#' @param stype a string c("paracancerous", "chemotherapy", "normal", "tumor", "puncture")
#' @param outputdir a string
#' @export
make_rds <- function(stype = "normal",
                     outputdir = "~/projects/hyperion/analysis/qzhang/human/steinbock/rds",
                     dataset = "qzhang",
                     overwrite = F,
                     filter_makers = c("80ArAr", "120Sn", "127I", "134Xe", "138Ba", "DNA1", "DNA2", "Histone3", "208Pb"),
                                       glob = NULL){
  if(is.null(glob)){
    rds_fn <- glue("{outputdir}/{stype}.rds")
  }else{
    suffix <- str_replace(glob, fixed("*"), "")
    rds_fn <- glue("{outputdir}/{stype}{suffix}.rds")
  }
  analysis_dir <- glue("~/projects/hyperion/analysis/{dataset}/human/steinbock")
  message("Processing: ", stype)
  if(fs::file_exists(rds_fn) && !overwrite){
    message(Sys.time(), "File ", rds_fn, " is existed. If you want to rerun it, please set overwrite = T")
  }else{
    if(is.null(glob)){
      para_samples <- dir_ls(glue("~/projects/hyperion/analysis/{dataset}/human/steinbock/{stype}/measure")) %>% path_file()
    }else{
      para_samples <- dir_ls(glue("~/projects/hyperion/analysis/{dataset}/human/steinbock/{stype}/measure"), glob = glob) %>% path_file()
    }
    dat_dir <- glue("{analysis_dir}/{stype}/measure")
    punc_sce <- seq_along(para_samples) %>%
                map( ~ imcRtools::read_steinbock(glue("{dat_dir}/{para_samples[.x]}"), return_as = "sce"))
    sce_lst <- punc_sce %>% map( ~ filter_makers(.x,filter_makers))
    sce_lst <- seq_along(para_samples) %>% map( ~ init_metadata(sce_lst[[.x]], para_samples[.x], stype = stype))
    sce <- purrr::reduce(sce_lst, BiocGenerics::cbind)
    sce <- sce[, colSums(counts(sce)) > 0]
    sce <- sce %>% remove_duplicate() %>% filter_cells_by_area() %>% filter_rare_marker_exp() %>%
           spillover_correction() %>% percentile_hyp() %>% normalize_hyp(method = "log2counts") %>%
           normalize_hyp(method ="lognormcounts") %>% imcAsinh_trans() %>% simpleAsinh_trans() %>%
           value_censor(assay = "log2counts") %>% zscore_censor()
    message(Sys.time(), "Writing rds file: ", rds_fn)
    readr::write_rds(sce, file = rds_fn)
    return(invisible(sce))
  }
}
