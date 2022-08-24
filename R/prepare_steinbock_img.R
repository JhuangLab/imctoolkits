#' @title prepare_steinbock_imgcsv
#' @name prepare_steinbock_imgcsv
#' @description this is a quick an dirty function.
#' @export
prepare_steinbock_imgcsv <- function(dataset = "qzhang",
                                     stypes = c("paracancerous", "chemotherapy", "normal", "tumor", "puncture")){
  panel_fn <- glue("~/projects/hyperion/data/{dataset}/human/hyperion/panel.csv")
  for (stype in stypes){
    dat_dir <- glue("~/projects/hyperion/analysis/{dataset}/human/steinbock/{stype}/measure")
    dirs <- dir_ls(dat_dir) %>%
      purrr::walk(~link_create(panel_fn, glue("{.x}/panel.csv")))
  }
  ##### combind all images.csv file
  img_path <- glue("~/projects/hyperion/data/{dataset}/human/hyperion/array/img")
  dat_img_files <- dir_ls(img_path, recurse = TRUE, glob = "*images.csv")
  dat_img <- read_csv(dat_img_files, show_col_types = FALSE, progress = F) %>%
             dplyr::distinct(image, .keep_all = TRUE)
  #stypes <- c("paracancerous", "chemotherapy", "normal", "tumor", "puncture")
  for (stype in stypes){
    dat_dir <- glue("~/projects/hyperion/analysis/{dataset}/human/steinbock/{stype}/measure")
    dir_name <- "ilastik/ilastik_probabilities"
    samples <- dir_ls(dat_dir, recurse = F) %>% map_chr(~glue("{.x}/{dir_name}"))
    for(sample_id in samples){
      message("Processing file: ", sample_id)
      csv_files <- dir_ls(sample_id, glob = "*.tiff") %>% path_file() %>% str_replace(".csv", ".tiff")
      fil <- dat_img$image  %in% csv_files
      dat_img_fil <- dat_img[fil, ] %>% dplyr::distinct(image, .keep_all = TRUE)
      fil <- ! csv_files %in% dat_img$image
      if(sum(fil) > 0){
        csv_files_rm <- csv_files[fil] %>% str_replace(".tiff", ".csv")
        for(d in c("intensities", "neighbors", "regionprops")){
          sample_id_tmp <- sample_id %>% str_replace(dir_name, d)
          rm_fn <- glue("{sample_id_tmp}/{csv_files_rm}")
          if( fs::file_exists(rm_fn) ){
            message(rm_fn)
            file_delete(rm_fn)
          }
        }
      }
      sample_name <- path_dir(sample_id) %>% path_dir() %>% path_file()
      fn_out <- glue("~/projects/hyperion/analysis/{dataset}/human/steinbock/{stype}/measure/{sample_name}/images.csv")
      write_csv(dat_img_fil, fn_out)
    }
  }
}

