#' @title pairwise_comparison
#' @name pairwise_comparison using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- pairwise_comparison(sce)
#' sce <- pairwise_comparison(sce)
#'}

pairwise_comparison <- function(feature_dat, class_info, method = "t") {
  # feature_dat is m x n matrix or data.frame, each row indicates one feature
  # class_info is 2 columns data.frame with matched samples (colnames of feature_dat) and classifications
  # method could be one in c("t", "f", "pairt", "blockf", "t.equalvar", "wilcoxon")
  match_index <- na.omit(match(unlist(class_info[, 1]), colnames(feature_dat)))
  stopifnot("The inputs has no overlap observations" = length(match_index) > 0)
  matched_feature_dat <- feature_dat[, match_index, drop = FALSE]
  match_index2 <- na.omit(match(colnames(feature_dat), unlist(class_info[, 1])))
  matched_class_info <- class_info[match_index2, ]
  ordered_class_info <- matched_class_info[match(colnames(matched_feature_dat), matched_class_info[, 1]), ]
  sample_class <- unlist(ordered_class_info[, 2])

  all_clusters <- unique(sample_class)
  comb <- combn(all_clusters, 2)
  k <- ncol(comb)
  pvalues <- data.frame(matrix(
    ncol = nrow(feature_dat), nrow = k,
    dimnames = list(NULL, rownames(feature_dat))
  ))
  rownames_comp <- character()
  for (i in 1:k) {
    x_sub <- matched_feature_dat[, sample_class %in% comb[, i], drop = FALSE]
    y_sub <- sample_class[sample_class %in% comb[, i]]
    comp_name <- paste(comb[, i], collapse = "_vs_")
    y_sub <- as.numeric(as.factor(y_sub)) - 1
    pvals <- mt.teststat(x_sub, y_sub, test = method)
    names(pvals) <- rownames(x_sub)
    rawp1 <- 2 * (1 - pnorm(abs(pvals)))
    rownames_comp[i] <- comp_name
    pvalues[i, ] <- rawp1
  }
  rownames(pvalues) <- rownames_comp
  return(pvalues)
}

#' @title getdiff_hyp
#' @name getdiff_hyp using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- getdiff_hyp(sce)
#' sce <- getdiff_hyp(sce)
#'}

getdiff_hyp <- function(clinical_dat, keycol, feature_dat, method = "t") {
  for_comp <- feature_dat[, colnames(feature_dat) %in% clinical_dat$patient_id]
  res <- pairwise_comparison(for_comp, clinical_dat[, c("patient_id", keycol)], method)
  return(res)
}

#' @title entropy_calc
#' @name entropy_calc using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- entropy_calc(sce)
#' sce <- entropy_calc(sce)
#'}
entropy_calc <- function(sce, keycol = "cell_type10"){
  stopifnot("input must be a SingleCellExperiment object" = "SingleCellExperiment" %in% class(sce))
  suppressMessages(
    frac_cells <- colData(sce) %>% as.data.frame() %>% group_by(sample_tiff_id) %>%
      mutate(allcount = n()) %>% group_by(across(c(keycol, "sample_tiff_id", "allcount"))) %>%
      summarise(count = n()) %>% na.omit() %>% mutate(frac = count / allcount))
  suppressMessages(
    frac_wide <- frac_cells %>% select(any_of(c(keycol, "sample_tiff_id", "frac"))) %>%
      pivot_wider(names_from = sample_tiff_id, values_from = frac, values_fill = 0) %>%
      as.data.frame() %>% column_to_rownames(var = keycol) %>% t())
  h_index <- vegan::diversity(frac_wide)
  temp <- data.frame(sample_tiff_id = names(h_index), H = h_index)
  temp2 <- as.data.frame(colData(sce))
  temp2 <- left_join(temp2, temp)
  colData(sce)[, paste0("entropy_", keycol)] <- temp2$H
  return(sce)
}


#' @title quick_violin
#' @name quick_violin using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- quick_violin(sce)
#' sce <- quick_violin(sce)
#'}
quick_violin <- function(data, x, y, z = NULL, ymin = -0.2, ymax = 0.9, step = 0.2, title) {
  ll <- combn(unique(na.omit(unlist(data[, x]))), 2, simplify = FALSE)
  plt <- ggviolin(data, x = x, y = y, facet.by = z,
                  color = x, add = "jitter", draw_quantiles = 0.5,
                  title = title, xlab = F, ylab = F, size = 1,
                  ggtheme = theme(panel.background = element_blank(),
                                  panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",
                                                                                               size = rel(1)), legend.key = element_blank(),
                                  strip.background = element_rect(fill = "white", colour = "black",
                                                                  size = rel(2)), complete = TRUE,
                                  text = element_text(size = 13),
                                  plot.title = element_text(size = 15, hjust = 0.02),
                                  legend.position = "none")) +
    # scale_y_continuous(expand = expansion(mult = c(0, z))) +
    stat_compare_means(comparisons = ll, na.rm = T, method = "t.test", size = 3, vjust = 0.1,
                       label.y = seq(ymax, length.out = length(ll), by = step)) +
    ylim(ymin, ymax + step * length(ll) + 0.1)
  return(plt)
}

#' @title quick_violin2
#' @name quick_violin2 using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- quick_violin2(sce)
#' sce <- quick_violin2(sce)
#'}
quick_violin2 <- function(data, x, y = "sparse_score", z = "cell_type10") {
  ll <- combn(unique(na.omit(unlist(data[, x]))), 2, simplify = FALSE)
  plt <- ggviolin(data, x = x, y = y, facet.by = z,
                  color = x, add = "jitter", draw_quantiles = 0.5,
                  title = "Sparse score of cells among groups", xlab = F, ylab = F, size = 1,
                  ggtheme = theme(panel.background = element_blank(),
                                  panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",
                                                                                               size = rel(1)), legend.key = element_blank(),
                                  strip.background = element_rect(fill = "white", colour = "black",
                                                                  size = rel(2)), complete = TRUE,
                                  text = element_text(size = 13),
                                  plot.title = element_text(size = 15, hjust = 0.02),
                                  legend.position = "none")) +
    # scale_y_continuous(expand = expansion(mult = c(0, z))) +
    stat_compare_means(comparisons = ll, na.rm = T, method = "t.test", size = 3, vjust = 0.1,
                       label.y = seq(0.9, length.out = length(ll), by = 0.2)) +
    ylim(-0.2, 1 + 0.2 * length(ll))
  return(plt)
}


#' @title temp_plot
#' @name temp_plot using this to avoid function name confiction.
#' @param sce object, a singcellexperiment object.

#' @export
#' @return A data frame.
#' @examples
#'\dontrun{
#' sce <- temp_plot(sce)
#' sce <- temp_plot(sce)
#'}
temp_plot <- function(data, list, u, group_var) {
  data2 <- list[[u]]
  names <- names(list)
  entropy <- left_join(data, data2[, c("patient_id", group_var)])
  entropy <- entropy[!is.na(entropy[, group_var]), ]
  # en_p <- entropy[, c("value", "value2")]
  # en_p <- as.data.frame(t(en_p))
  # colnames(en_p) <- entropy$sample_tiff_id
  # xx <- pairwise_comparison(en_p, entropy[, c("sample_tiff_id", group_var)])
  # return(xx)

  h_vio <- quick_violin(entropy, group_var, "value")
  h_vio2 <- quick_violin(entropy, group_var, "value2")

  ggsave(glue("{out_dir}/entropy_{names[u]}_{group_var}_2.png"),
         h_vio2, dpi=300, width = 5, height = 5, units = "in")
  ggsave(glue("{out_dir}/entropy_{names[u]}_{group_var}.png"),
         h_vio, dpi=300, width = 5, height = 5, units = "in")
}

