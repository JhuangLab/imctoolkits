#' @title cytomapper_plot
#' @name cytomapper_plot
#' @param dat_dir
#' @param pdf_fn
#' @param sample_id
#' @param colour_by colour_by = c("Collagen1", "Vimentin", "DNA1", "E_Cadherin", "a_smooth")
#' @param bcg bcg = list(Collagen1 = c(0, 20, 1), Vimentin = c(0, 20, 1), DNA1 = c(0, 2, 1), E_Cadherin = c(0, 20, 1), a_smooth = c(0, 20, 1))
#' @export
cytomapper_plot <- function(dat_dir, pdf_fn, sample_id, colour_by = c("CD3", "CD20"),
                            bcg = list(CD3 = c(0, 20, 1), CD20 = c(0, 20, 1)),
                            display = "single",
                            panel_fn = "~/projects/hyperion/data/qzhang/human/hyperion/panel.csv",
                            image_title = NULL){
  # img_path <- file.path(glue("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/{type}"), "img", x)
  # mask_path <- file.path(glue("/cluster/home/jhuang/projects/hyperion/analysis/qzhang/human/steinbock/{type}/measure"), x, "masks")
  panel <- fread(panel_fn)
  img_path <- glue("{dat_dir}/img/{sample_id}")
  mask_path <- glue("{dat_dir}/measure/{sample_id}/masks")
  if(file.info(mask_path)$size != 0){
    # img
    img <- cytomapper::loadImages(img_path, BPPARAM = BiocParallel::MulticoreParam(80), pattern = "tiff")
    cytomapper::channelNames(img) <- panel$name
    mcols(img) <- data.frame(sample_id = names(img))
    # mask
    mask <- cytomapper::loadImages(mask_path, BPPARAM = BiocParallel::MulticoreParam(80)) %>%
            cytomapper::scaleImages(2^16-1)
    mcols(mask) <- data.frame(sample_id = names(mask))
    # plot
    pdf(pdf_fn, width = 10, height = 8)
      plotPixels(img,
                 mask = mask,
                 img_id = "sample_id",
                 # colour_by = c("Collagen1", "Vimentin", "DNA1", "E_Cadherin", "a_smooth"),
                 colour_by = colour_by,
                 bcg = bcg,
                 display = display,
                 image_title = image_title)
    dev.off()
  }
}
#' @title cytomapper_cell_plot
#' @name cytomapper_cell_plot
#' @param dat_dir
#' @param pdf_fn
#' @param sample_id
#' @param colour_by colour_by = c("Collagen1", "Vimentin", "DNA1", "E_Cadherin", "a_smooth")
#' @param bcg bcg = list(Collagen1 = c(0, 20, 1), Vimentin = c(0, 20, 1), DNA1 = c(0, 2, 1), E_Cadherin = c(0, 20, 1), a_smooth = c(0, 20, 1))
#' @export
cytomapper_cell_plot <- function(sce, sample_id,
                                 dat_dir, colour_by = "cluster",
                                 pdf_fn,
                                 panel_fn = "~/projects/hyperion/data/qzhang/human/hyperion/panel.csv"){
  img_path <- glue("{dat_dir}/img/{sample_id}")
  mask_path <- glue("{dat_dir}/measure/{sample_id}/masks")
  # img
  panel <- fread(panel_fn)
  img <- loadImages(img_path,
                    BPPARAM = BiocParallel::MulticoreParam(80))
  cytomapper::channelNames(img) <- panel$name
  mcols(img) <- data.frame(sample_id = names(img))

  # mask
  mask <- loadImages(mask_path, BPPARAM = BiocParallel::MulticoreParam(80)) %>%
    scaleImages(2^16-1)
  mcols(mask) <- data.frame(sample_tiff_id = names(mask))
  # subset sce
  sce_sub <- sce[, colData(sce)$sample_tiff_id %in% names(img)]
  # plot
  pdf(pdf_fn)
    plotCells(mask, object = sce_sub,
              cell_id = "ObjectNumber",
              img_id = "sample_tiff_id",
              colour_by = "cluster",
              display = "single",
              image_title = NULL,
              legend = list(colour_by.labels.cex = 0.7,
                            colour_by.legend.cex = 0.7))
  dev.off()
}
