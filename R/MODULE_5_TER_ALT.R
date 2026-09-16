#' @importFrom data.table data.table setnames setkey setcolorder setorder tstrsplit CJ .N .SD .EACHI :=
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>% count rename filter_all all_vars
#' @import corrplot
#' @import rlist
#' @importFrom nlme gls corCompSymm
#' @importFrom gridExtra grid.arrange
#' @import EnhancedVolcano
#' @import fitdistrplus


#' @title generate_ENZ
#' @description Function to generate empirical null distribution.
#' @param x A sample-by-sample list of RNA and RPF count data and sample attributes produced by \code{\link{partition_to_uniques}}.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' @param outcome The variable determining whether a vector of read counts is RNA or RPF.
#' This is usually the name of the response variable in the TER test logistic regression performed through \code{\link{logit_seq}}. Default: \code{"read_type"}.
#' @param uniqueID A variable (column) of the design matrix defining unique experimental preparations
#' from each of which one RNA sample and one RPF sample was derived. It corresponds to the highest resolution
#' (lowest level) of classification of samples in the data set apart from the RNA/RPF distinction
#' and is usually equal to replicate name in biological experiments.
#' @param groupID A variable (column) of the design matrix indicating which replicates should be grouped together.
#' All experimental units having the same \code{groupID} will be considered replicates of the same biological sample
#' (or members of the same group of samples).
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return
#' A vector of z statistics constituting empirical null.
#' @details
#' In large scale hypothesis testing e.g. genomic data sets, it may be possible to observe the null distribution, instead
#' of relying the theoretically assumed distribution (standard normal for regression). Ribolog compares replicates of
#' each biological sample (items with the same \code{groupID}) and pools the z values from those regressions to produce the
#' empirical null.
#' @examples
#' rr.v2_dummy.split <- Ribolog::partition_to_uniques(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, "replicate_name")
#' rr.v2_dummy.enz <- generate_ENZ(x = rr.v2_dummy.split, design = Ribolog::sample_attributes_dummy, outcome = "read_type",
#'                            uniqueID = "replicate_name", groupID = "cell_line", adj_method = "none")
#' @export
generate_ENZ <- function(x, design, outcome = "read_type", uniqueID, groupID, adj_method){

  homo_pair_z_vectors <- list()
  n <- length(x)
  n_design_cols <- dim(design)[2]

  for (i in c(2:n)){
    for (j in c(1:(i-1))){
      list_ij <- list()
      list_ij[["uniqueIDs"]] <- sort(c(names(x)[i], names(x)[j]))

      list_ij[["groupIDs"]] <- sort(c(as.character(x[[i]][, groupID][1]), as.character(x[[j]][, groupID][1])))
      if (identical(list_ij[["groupIDs"]][1], list_ij[["groupIDs"]][2])) {
        list_ij[["pair_type"]] = "homo"
        model1 <- as.formula(paste(outcome, uniqueID, sep = "~"))
        data_ij <- rbind(x[[i]], x[[j]])[, -c(1:n_design_cols)]
        design_ij <- rbind(x[[i]], x[[j]])[, c(1:n_design_cols)]
        list_ij[["z"]] <- Ribolog::logit_seq(t(data_ij), design_ij, model1, adj_method=adj_method, long_output = TRUE)[,7]
        name_ij <- paste(list_ij[["uniqueIDs"]], collapse = "_vs_")
        homo_pair_z_vectors[[name_ij]] <- list_ij
      }
    }
  }
  # extract the z vectors and unlist into a single vector
  homo_pair_z_vectors_merged <- as.vector(unlist(lapply(homo_pair_z_vectors, function(x) x[["z"]])))
  return(homo_pair_z_vectors_merged)
}


#' @title build_empirical_null
#' @description Convenience wrapper that generates the empirical null (EN) distribution of TER test z
#' statistics in one call, from a combined RNA+RPF count matrix and its design matrix. Combines
#' \code{\link{partition_to_uniques}} (partitioning the data into unique replicate preparations) and
#' \code{\link{generate_ENZ}} (pooling z statistics from replicate-vs-replicate "homo" comparisons) -
#' see their documentation for details on the underlying mechanics.
#' @param x A data matrix of RNA and RPF read counts (transcript ID column excluded), as used by \code{\link{logit_seq}}.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' @param uniqueID A variable (column) of \code{design} defining unique experimental preparations, from each
#' of which one RNA sample and one RPF sample was derived (usually a replicate name).
#' @param groupID A variable (column) of \code{design} indicating which replicates should be grouped together
#' as replicates of the same biological sample. All experimental units sharing the same \code{groupID} are
#' compared to each other to generate the empirical null.
#' @param outcome The variable determining whether a vector of read counts is RNA or RPF. Default: "read_type".
#' @param adj_method P-value adjustment method passed to the underlying \code{\link{logit_seq}} calls.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default: "none".
#' @return A vector of z statistics constituting the empirical null (same as \code{\link{generate_ENZ}}'s output).
#' @examples
#' rr.v2_dummy.enz <- build_empirical_null(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, uniqueID = "replicate_name", groupID = "cell_line")
#' @export

build_empirical_null <- function(x, design, uniqueID, groupID, outcome = "read_type", adj_method = "none"){
  x_split <- Ribolog::partition_to_uniques(x = x, design = design, uniqueID = uniqueID)
  enz <- Ribolog::generate_ENZ(x = x_split, design = design, outcome = outcome,
                                uniqueID = uniqueID, groupID = groupID, adj_method = adj_method)
  return(enz)
}


#' @title test_ENZ
#' @description Function to calculate p-values by empirical null hypothesis testing.
#' @param x Regression fit object produced by \code{\link{logit_seq}} with option \code{long_output = TRUE},
#' or any data frame or matrix containing a column
#' or columns of test statistics (z).
#' @param zcols Column(s) in \code{x} containing test z statistics.
#' @param enz A vector of z values constituting empirical null e.g. produced by \code{\link{generate_ENZ}}.
#' @param method Can be "direct" or "fit". "direct": two-sided p-value of a z is calculated as
#' the percentage of ENZ larger than abs(z) or smaller than -abs(z). "fit": a normal distribution is
#' fitted to the empirical null and is used to obtain two-sided p-values. Default: \code{"fit"}.
#' @param keep_data A logical variable indicating the original data set should be included in the output. If it is set to
#' FALSE, only the empirical null p is reported back. Default: TRUE.
#' @details
#' Function \code{\link{logit_seq}} outputs regression coefficients and p-values by defaults. Option
#' \code{long_output = TRUE} adds SD(beta) and z columns to the output matrix.
#' Output column names are derived from \code{colnames(x)[zcols]} by replacing "zvalue" with "ep"
#' (e.g. \code{zvalue_lung_metastasisY} -> \code{ep_lung_metastasisY}, matching \code{\link{logit_seq}}'s
#' naming convention); columns not following that convention are returned unrenamed.
#' @examples
#' fit1 <- Ribolog::logit_seq(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                             as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' rr.v2_dummy.enz <- Ribolog::build_empirical_null(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy,
#'                                             uniqueID = "replicate_name", groupID = "cell_line")
#' fit1_ENZ_p <- test_ENZ(x = fit1, enz = rr.v2_dummy.enz, zcols = 7)
#' @return
#' Vector(s) of p-values calculated from an empirical null, reported alone or appended to the original data set.
#' @export
test_ENZ <- function(x, zcols, enz, method = "fit", keep_data = TRUE){

  fitz <- x[, zcols, drop = FALSE]
  x_ecdf <- ecdf(enz)

  if (method == "direct"){
    fitz_enz_p <- apply(fitz, 1:2, function(y) ( 1 - x_ecdf(abs(y)) + x_ecdf((-1)*abs(y))))
  }
  if (method == "fit"){
    enz_normal_fit <- fitdistrplus::fitdist(enz, "norm")
    mean_fitz <- as.vector(enz_normal_fit$estimate[1])
    sd_fitz <- as.vector(enz_normal_fit$estimate[2])
    fitz_std <- apply(fitz, 1:2, function(y) ((y - mean_fitz) / sd_fitz))
    fitz_enz_p <- apply(fitz_std, 1:2, function(y) ( 2 * (1 - pnorm(abs(y)))))
  }
  colnames(fitz_enz_p) <- gsub("zvalue", "ep", colnames(fitz_enz_p))
  if (keep_data == TRUE){
    output <- data.frame(x, fitz_enz_p)
  } else if (keep_data == FALSE) {
    output <- fitz_enz_p
  }
  return(output)
}


#' @title visualize_pvalue_comparison
#' @description Scatter plot comparing theoretical vs empirical p-values (e.g. the output of
#' \code{\link{test_ENZ}}) on a log-log scale, with a y = x reference line - the same comparison shown in
#' the empirical null testing vignette, formalized into a reusable, themed plot.
#' @param x_data Data frame containing both a theoretical and an empirical p-value column, such as the
#' output of \code{\link{test_ENZ}} run with \code{keep_data = TRUE}.
#' @param covariate (Optional) Name (or unambiguous prefix) of the model term to plot, e.g.
#' "lung_metastasis" - matched against the \code{pvalue_<term>} columns of \code{x_data} (same matching
#' logic as \code{\link{volcano_plot}}'s \code{covariate} argument). When given, sets \code{x_col} to
#' \code{pvalue_<term>} and \code{y_col} to \code{ep_<term>}, overriding them. Leave \code{NULL} to specify
#' \code{x_col}/\code{y_col} manually instead.
#' @param x_col Column of \code{x_data} holding the theoretical p-value. Ignored if \code{covariate} is given.
#' @param y_col Column of \code{x_data} holding the empirical p-value. Ignored if \code{covariate} is given.
#' @param title Plot title. Default: "Theoretical vs Empirical P-values".
#' @param subtitle Plot subtitle. Default: \code{NULL}.
#' @param xlab,ylab Axis labels. Default: "Theoretical P" / "Empirical P".
#' @param xlim,ylim (Optional) Length-2 numeric vectors giving the axis display range, in raw
#' (untransformed) p-value units, e.g. \code{xlim = c(1e-40, 1)}. Passed straight through to
#' \code{ggplot2::coord_cartesian}, which applies the plot's own log10 scale transform internally - so,
#' like \code{\link{visualise_empirical_null}}'s \code{xlim} - this only zooms the view; it never drops data
#' or distorts the plotted points, unlike \code{ggplot2::xlim}/\code{scale_x_log10(limits =)}.
#' @param title_size Font size for the title. Default: 22.
#' @param axis_title_size Font size for the (bold) axis titles. Default: 18.
#' @param axis_text_size Font size for the (bold) axis tick labels. Default: 14.
#' @param point_color Colour of the p-value points. Default: "#1F3864".
#' @param point_size,point_alpha Size/transparency of the p-value points. Defaults: 0.6, 0.5.
#' @param ref_color Colour of the y = x reference line. Default: "#C0392B".
#' @param legend_position Where to place the legend. Default: "top".
#' @param legend_text_size Font size for the (bold) legend labels. Default: 14.
#' @return A ggplot2 object.
#' @details
#' Both axes are on a log10 scale, since p-values span many orders of magnitude. Exact-zero p-values
#' (common for very significant hits, from floating point underflow) can't be log-transformed and would
#' otherwise be silently dropped with a "values <= 0 omitted" warning; instead they are floored to
#' one-tenth of the smallest non-zero value in their column before plotting (the same convention
#' EnhancedVolcano uses for exact-zero p-values), so every point stays visible.
#' @examples
#' ter_testing_results <- Ribolog::logit_seq(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                                            as.vector(rr.v2_dummy$transcript), adj_method = 'fdr')
#' rr.v2_dummy.enz <- Ribolog::build_empirical_null(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy,
#'                                             uniqueID = "replicate_name", groupID = "cell_line")
#' ter_testing_ENZ_p <- Ribolog::test_ENZ(x = ter_testing_results, enz = rr.v2_dummy.enz, zcols = 8)
#' Ribolog::visualize_pvalue_comparison(ter_testing_ENZ_p, covariate = "lung_metastasis")
#' Ribolog::visualize_pvalue_comparison(ter_testing_ENZ_p, covariate = "lung_metastasis",
#'                                       title = "TER test calibration", subtitle = "LMCN dataset",
#'                                       xlab = "Theoretical p-value", ylab = "Empirical p-value",
#'                                       xlim = c(1e-40, 1), ylim = c(1e-40, 1),
#'                                       axis_title_size = 20, axis_text_size = 16,
#'                                       legend_text_size = 16)
#' @export

visualize_pvalue_comparison <- function(x_data, covariate = NULL,
                                         x_col = "pvalue_lung_metastasisY", y_col = "ep_lung_metastasisY",
                                         title = "Theoretical vs Empirical P-values", subtitle = NULL,
                                         xlab = "Theoretical P", ylab = "Empirical P",
                                         xlim = NULL, ylim = NULL,
                                         title_size = 22, axis_title_size = 18, axis_text_size = 14,
                                         point_color = "#1F3864", point_size = 0.6, point_alpha = 0.5,
                                         ref_color = "#C0392B", legend_position = "top",
                                         legend_text_size = 14){

  if (!is.null(covariate)) {
    pcols <- grep("^pvalue_", names(x_data), value = TRUE)
    terms <- sub("^pvalue_", "", pcols)
    match_idx <- which(startsWith(terms, covariate))

    if (length(match_idx) == 0) {
      stop(paste0("No pvalue column found matching covariate '", covariate, "'. Available terms: ",
                  paste(terms, collapse = ", ")))
    }
    if (length(match_idx) > 1) {
      stop(paste0("covariate '", covariate, "' matches more than one term: ",
                  paste(terms[match_idx], collapse = ", "), ". Use a more specific value."))
    }

    term <- terms[match_idx]
    x_col <- paste0("pvalue_", term)
    y_col <- paste0("ep_", term)
  }

  if (! x_col %in% names(x_data)) {
    stop(paste0("The column '", x_col, "' does not exist in x_data. Available columns: ",
                paste(names(x_data), collapse = ", ")))
  }
  if (! y_col %in% names(x_data)) {
    stop(paste0("The column '", y_col, "' does not exist in x_data. Available columns: ",
                paste(names(x_data), collapse = ", ")))
  }

  floor_zero <- function(v){
    v <- as.numeric(v)
    nz <- v[v > 0 & !is.na(v)]
    if (length(nz) > 0) {
      v[!is.na(v) & v <= 0] <- min(nz) / 10
    }
    v
  }

  plot_df <- data.frame(
    theoretical = floor_zero(x_data[[x_col]]),
    empirical = floor_zero(x_data[[y_col]])
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = theoretical, y = empirical)) +
    ggplot2::geom_point(ggplot2::aes(colour = "P-values"), size = point_size, alpha = point_alpha) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1, linetype = "y = x"), colour = ref_color) +
    ggplot2::scale_colour_manual(name = NULL, values = c("P-values" = point_color)) +
    ggplot2::scale_linetype_manual(name = NULL, values = c("y = x" = "dashed")) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1,
                      override.aes = list(size = 3.5, alpha = 1)),
                    linetype = ggplot2::guide_legend(order = 2,
                      override.aes = list(linewidth = 1.3))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = axis_title_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5,
                                          margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = title_size * 0.6, colour = "grey40", hjust = 0.5,
                                             margin = ggplot2::margin(b = 12)),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text = ggplot2::element_text(size = axis_text_size, colour = "grey20", face = "bold"),
      axis.line = ggplot2::element_line(colour = "grey50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text_size, face = "bold"),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1.4, "lines")
    )

  # xlim/ylim are given in raw p-value units (matching the un-transformed data); coord_cartesian applies
  # the log10 scale's own transform internally, so they're passed through as-is here (pre-transforming
  # them too, e.g. with log10(), would double-transform and produce NaNs for any value < 1). This still
  # only zooms the view - it never drops data or re-triggers the "values <= 0 omitted" warning this
  # function already floors away.
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }

  p
}


#' @title comb_cor_p
#' @description Function to combine correlated p-values for one feature (row). Not exported; used only
#' internally by \code{\link{meta_test}}.
#' @param p Vector of p-values (one per test being combined) for a single feature (gene/transcript).
#' @param w Vector of weights (one per test being combined) for the same feature.
#' @param cm Correlation matrix of z-values or any other test statistics of the tests being combined (mxm,
#' m = number of tests).
#' @details
#' This is based on the the thoery expounded in a 2003 paper by Kepher Makambi.
#' @references
#' Makambi, K. 2003. Weighted inverse chi-square method for correlated significance tests.
#' Journal of Applied Statistics, 30(2): 225-234.
#' @return
#' A one-row data frame with three columns: \code{M} (combined test statistic), \code{nu} (degrees of
#' freedom of the chi-square distribution), \code{p_com} (the combined/meta-analytical p-value).
#' @noRd

comb_cor_p <- function(p, w, cm){
  M <- (-2) * sum(w * log(p))
  ww <- w %*% t(w)
  nu <- 8 / sum(3.25 * ww * cm + 0.75 * ww * cm^2)
  p_out <- pchisq((nu * M)/2, df = nu, lower.tail = FALSE)
  out <- data.frame(M, nu, p_out)
  colnames(out) <- c("M", "nu", "p_com")
  return(out)
}


#' @title meta_test
#' @description Function to combine the output of related tests in a meta-analytical framework.
#' @param x A list of test results to be combined. Each element of the list is a four-column matrix:
#' [1] beta (regression coefficient or any other measures of effect size) [2] sd.beta (standard deviation of beta),
#' [3] z (the test statistic from which original test p-values were calculated) [4] p (original test p-value).
#' In the Ribolog pipeline, it is produced by the \code{\link{logit_seq}} function with the \code{long_output = TRUE} option.
#' @param features_equivar Method used to calculate the correlation matrix among test statistics of the tests
#' that are being combined. FALSE: The correlation matrix on the provided vectors of test statistics is calculated
#' as is. Features (genes, transcripts, etc) with higher variability of results among tests will
#' the correlation matrix more. TRUE: row-normalize the test statistics first, so that the feature-wise e.g. gene-wise
#' variance of z's for all features = 1. All features contribute equally to the correlation matrix. Deafult: FALSE.
#' @param comp_sym Compound symmetry. TRUE: All pairs of tests are assumed to be equally correlated. If more than two
#' tests are being combined, all pairwise correlation coefficients are replaced by a common value (more in Details).
#' FALSE: Correlation matrix is calculated by the R default i.e. separately between each pair of z vectors.
#' This would be equivalent to unstructured covariance matrix in repreated measure ANOVA.
#' Default: FALSE.
#' @param effect_wt Indicator variable. TRUE: Combine effect sizes weighted by their inverse of variance.
#' FALSE: Combine effect sizes with equal weights. Default: TRUE.
#' @param feature_list (Optional) A vector containing IDs of genes/transcripts.
#' Must have the same length as the row number of input list data frames.
#' @param long_output Indicates whether the parameters M and nu of the Makambi method should be included in the output.
#' Default: FALSE.
#' @details
#' The columns in \code{x} do not need to be named exactly as mentioned in parameter description, but they must be
#' in that exact order. They must all be given, too. The combination of test results cannot be done with p-values alone,
#' because p-values do not contain information on the direction of the effect, and thus, cannot be used to calculate correlation
#' between tests. Correlation among tests is calculated from z scores, not p-values. The vector of tests statistics should
#' preferably come directly from the original test. But, many software tools report only a measure of
#' effect size and the p-value, and omit SDs or test statistics from the output for the sake of brevity. If the p-value was calculated
#' using a Wald test (estimated parameter divided by its SD and compared against the standard normal distribution), use the function
#' \code{\link{add_sd_z}} to reproduce the SD and z score columns, and then feed into \code{\link{meta_test}}. SD is used by
#' \code{\link{meta_test}} to calculate weights which is made to proprtional to the inverse of variance. If other measures of weights
#' are desired, manualy substitute the SD column in \code{x} with inverse square roots of relative weights. If you want to give all
#' samples equal weight regardless of estimation uncertainty (not recommended), replace the SD column with a vector of 1s.
#'
#' Compound symmetry is a covariance structure encountered in repeated measure ANOVA i.e. when a quantity
#' is measured multiple times on each subject. In meta-analysis, we are measuring the association test statistic, e.g.
#' z score from regression, multiple times for each transcript. Compound symmetry assumes a common error variance for all measurements,
#' and a common covariance between all pairs of measurements on the same subject. This leads to equal correlation
#' between all pairs of measurements on the same subject. In meta-analysis, this translates to equal correlation of all
#' the tests that are being combined. If this option is selected, a common correlation coefficient will be estimated
#' that best describes the average similarity between all pairs of z vectors.
#'
#' Note that the \code{\link{meta_test}} function is designed to combined tests on one predictor variable at a time. Ribolog and
#' some other bioanalytic tools can incorporate several predictors in one model. Also, models usually include an intercept term.
#' Results (the four previously described columns) for each predictor have to be extracted from the \code{\link{logit_seq}} output
#' and combined separately using \code{\link{meta_test}}.
#' @return
#' A two or four column matrix. By default, a two column matrix: [1] meta_beta: combined effect size [2] meta_p: combined p-value.
#' If \code{long_output = TRUE} is specified, a four columns matrix: [1] meta_beta: combined effect size [2] M [3] nu
#' [4] meta_p: combined p-value. \code{M} is the meta-analysis test statistic. A derivative of \code{M} follows a chi-squared distribution
#' with \code{nu} degrees of freedom.
#' @examples
#' # The effect of lung metastasis is tested in two genetic backgrounds separately, then combined:
#' fit1a <- Ribolog::logit_seq(rr.v2_dummy[,c(2:5,10:13)], Ribolog::sample_attributes_dummy[c(1:4,9:12), ],
#'                              read_type ~ lung_metastasis, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' fit1b <- Ribolog::logit_seq(rr.v2_dummy[,c(6:9,14:17)], Ribolog::sample_attributes_dummy[c(5:8,13:16), ],
#'                              read_type ~ lung_metastasis, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # We are not interested in the intercept terms, so keep only the last four columns of the outputs which pertain to lung metastasis:
#' fit1ab_list <- list("fit1a" = fit1a[, c(5:8)], "fit1b" = fit1b[, c(5:8)])
#' # And now, combine:
#' fit1ab_meta <- meta_test(fit1ab_list, feature_list = rr.v2_dummy$transcript)
#' @references
#' Makambi, K. 2003. Weighted inverse chi-square method for correlated significance tests.
#' Journal of Applied Statistics, 30(2): 225-234.
#' @export
meta_test <- function(x, features_equivar = FALSE, comp_sym = FALSE, effect_wt = TRUE,
                      feature_list = NULL, long_output = FALSE){
  n <- length(x)
  beta_mat <- do.call(cbind,lapply(x, function(y) y[,1]))
  weight_mat <- do.call(cbind,lapply(x, function(y) y[,2]))^-2
  weight_df <- t(data.frame(apply(weight_mat, 1, function(y) y/sum(y))))
  z_df <- data.frame(do.call(cbind,lapply(x, function(y) y[,3])))

  if (features_equivar == TRUE){
    z_df <- apply(z_df, 1, function(y) y/sd(y))
  }
  z_cor <- cor(z_df)
  print(z_cor)

  if (comp_sym == TRUE){
    z_df2 <- z_df
    z_df2$transcript <- factor(rownames(z_df2))
    z_df2l <- tidyr::gather(z_df2, z_source, z, -transcript, factor_key = TRUE)
    z_df2lg <- nlme::groupedData(z ~ 1| transcript, data = z_df2l)
    fit_cs_z <- gls(z ~ 1, data = z_df2lg, correlation = corCompSymm(, form= ~ 1 | transcript))
    csx <- fit_cs_z$modelStruct$corStruct[1]
    cs_rho <- (exp(csx) - 1 ) / (exp(csx) + 1)
    print(paste("Compound symmetry Rho:", cs_rho))
    z_cor <- matrix(cs_rho, nrow = nrow(z_cor), ncol = ncol(z_cor))
    diag(z_cor) <- 1
  }

  p_df <- data.frame(do.call(cbind,lapply(x, function(y) y[,4])))
  p_w_df <- data.frame(p_df, weight_df)

  if (effect_wt == FALSE){
    beta_vec <- rowMeans(beta_mat)
  } else if (effect_wt == TRUE){
    beta_vec <- rowSums(beta_mat * weight_df) / rowSums(weight_df)
  }

  M_nu_p <- do.call(rbind, apply(p_w_df, 1, function(y) {comb_cor_p(y[1:n], y[(n+1):(2*n)], z_cor)} ))
  result <- data.frame(beta_vec, M_nu_p)
  rownames(result) <- feature_list
  colnames(result) <- c("meta_beta", "M", "Nu", "meta_p")
  if (long_output == "FALSE"){
    result <- result[, c(1,4)]
  }
  return(result)
}


#' @title add_sd_z
#' @description Function to add standard deviation and z statistic to test results
#' @param x A matrix representing output from a regression test or similar, with two columns per predictor (may include intercept).
#' The first column contains a measure of effect size e.g. regression coefficient. The second column contains p-values. The matrix
#' follows a [, beta1 p1 beta2 p2 ...] format.
#' @details
#' Regression fitting algorithms calculate a standard deviation (SD) for each coefficient (beta) and then calculate a z staistic by dividing
#' each coefficient by its SD. P-value is obtained by comapring the  statistic to the standard normal distribution. SD and z are usually
#' not reported in the regression output because they can be recovered easily from beta and p. Z is required for empirical null testing.
#' Both SD and z are required for the correlated tests meta-analysis. Function \code{\link{add_sd_z}} adds these two components to the
#' (beta & p) outputs so that they can be used as input to functions such as \code{\link{test_ENZ}} or \code{\link{meta_test}}.
#' If p-value was calculated from a distribution other than standard normal, the \code{\link{add_sd_z}} function should not be used.
#' The user must add the SD and z columns to the output manually using apprporiate operations.
#' @return
#' An augmented matrix where SD(beta) and z columns are added between beta and p for each predictor (and intercept, if applicable).
#' @examples
#' fit1 <- Ribolog::logit_seq(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                             as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' fit1_short <- fit1[, c(TRUE, FALSE, FALSE, TRUE)]
#' fit1_restored <- add_sd_z(fit1_short)
#' @export

add_sd_z <- function(x){
  b <- x[, c(TRUE, FALSE)]
  p <- x[, c(FALSE, TRUE)]
  z <- qnorm( p / 2) * sign(b) * (-1)
  sd <- b / z
  x_out <- t(gdata::interleave(t(b),t(sd),t(z),t(p)))
  colnames(x_out) <-c(as.vector(gdata::interleave(colnames(x)[c(TRUE, FALSE)],
                      gsub("Estimate.", "", paste0("sd.", colnames(x)[c(TRUE, FALSE)])),
                      gsub("Estimate.", "", paste0("z.", colnames(x)[c(TRUE, FALSE)])),
                      colnames(x)[c(FALSE, TRUE)])))
  return(x_out)
}


#' @title validate_meta_test
#' @description Test one predictor's effect separately within each level of a splitting variable, combine
#' the split-level results with \code{\link{meta_test}}, and compare the combined (meta) and per-split
#' results against a single test on the fully pooled data - the "gold standard" this alternative approach
#' is meant to approximate. Automates the workflow described in Ribolog's empirical null / meta-analysis
#' vignette (module 5, "same biological question, separate datasets").
#' @param x A data matrix of RNA and RPF read counts (transcript ID column excluded), as used by \code{\link{logit_seq}}.
#' @param design Design matrix of the experiment describing samples and their attributes. Must contain \code{split_var}.
#' @param model Regression model formula passed to \code{\link{logit_seq}} (e.g. \code{read_type ~ lung_metastasis}).
#' @param split_var Name of the column in \code{design} used to split \code{x}/\code{design} into independent
#' sub-datasets to be tested separately and then combined (e.g. two labs, protocols, or batches). Each level
#' of \code{split_var} is tested on its own with \code{\link{logit_seq}}.
#' @param term_name Name of the model term whose coefficient/p-value should be compared and combined (e.g.
#' \code{"lung_metastasisY"}); matched against the suffix of the relevant \code{\link{logit_seq}} long-output columns.
#' @param feature_list (Optional) A vector containing IDs of genes/transcripts, passed to \code{\link{logit_seq}} and \code{\link{meta_test}}.
#' @param adj_method P-value adjustment method passed to the underlying \code{\link{logit_seq}} calls. Must be
#' "none": \code{\link{meta_test}} needs the raw (unadjusted) z and p from each split-level test, and this
#' function does not currently support pre-adjusted p-values as input. Adjust \code{meta_p} yourself afterward
#' if desired. Default: "none".
#' @param cor_type Correlation method passed to \code{Hmisc::rcorr} when comparing pooled, split and meta
#' results: "pearson" or "spearman". Default: "spearman".
#' @param comp_sym Passed to \code{\link{meta_test}}. Default: FALSE.
#' @param effect_wt Passed to \code{\link{meta_test}} (its \code{wt_effects} argument). Default: TRUE.
#' @return A list with elements: \code{pooled_fit} (the gold-standard \code{\link{logit_seq}} fit on the full
#' data), \code{split_fits} (a named list of per-level \code{\link{logit_seq}} fits), \code{meta_fit} (the
#' \code{\link{meta_test}} combination of the split fits), \code{betas} and \code{ps} (data frames comparing
#' the pooled, per-split and meta effect sizes / p-values for \code{term_name}), and \code{beta_cor} / \code{p_cor}
#' (\code{Hmisc::rcorr} objects computed on \code{betas} / \code{ps}, ready to pass to \code{corrplot::corrplot}).
#' @examples
#' val <- validate_meta_test(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                            split_var = "cell_line_origin", term_name = "lung_metastasisY",
#'                            feature_list = as.vector(rr.v2_dummy$transcript))
#' corrplot::corrplot(val$beta_cor$r, method = "color", addCoef.col = "white")
#' @export

validate_meta_test <- function(x, design, model, split_var, term_name, feature_list = NULL,
                                adj_method = "none", cor_type = "spearman", comp_sym = FALSE, effect_wt = TRUE){

  pooled_fit <- Ribolog::logit_seq(x, design, model, feature_list = feature_list,
                                    long_output = TRUE, adj_method = adj_method)

  levels <- unique(as.character(design[[split_var]]))
  split_fits <- lapply(levels, function(lvl){
    idx <- which(design[[split_var]] == lvl)
    Ribolog::logit_seq(x[, idx], design[idx, ], model, feature_list = feature_list,
                        long_output = TRUE, adj_method = adj_method)
  })
  names(split_fits) <- levels

  term_cols <- grep(paste0("_", term_name, "$"), colnames(pooled_fit))
  if (length(term_cols) != 4) {
    stop(paste0("Expected exactly 4 columns (Estimate, SD, z, p) matching term '", term_name,
                "' in the logit_seq long output, found ", length(term_cols),
                ". Check that term_name matches a real model term and that adj_method = 'none'."))
  }

  meta_input <- lapply(split_fits, function(f) f[, term_cols])
  meta_fit <- Ribolog::meta_test(meta_input, comp_sym = comp_sym, effect_wt = effect_wt, feature_list = feature_list)

  beta_col <- term_cols[1]
  p_col <- term_cols[4]

  beta_mat <- cbind(pooled = pooled_fit[, beta_col],
                     do.call(cbind, lapply(split_fits, function(f) f[, beta_col])),
                     meta = meta_fit$meta_beta)
  colnames(beta_mat) <- c("pooled", names(split_fits), "meta")
  betas <- as.data.frame(beta_mat)

  p_mat <- cbind(pooled = pooled_fit[, p_col],
                  do.call(cbind, lapply(split_fits, function(f) f[, p_col])),
                  meta = meta_fit$meta_p)
  colnames(p_mat) <- c("pooled", names(split_fits), "meta")
  ps <- as.data.frame(p_mat)

  beta_cor <- Hmisc::rcorr(as.matrix(betas), type = cor_type)
  p_cor <- Hmisc::rcorr(as.matrix(ps), type = cor_type)

  list(pooled_fit = pooled_fit, split_fits = split_fits, meta_fit = meta_fit,
       betas = betas, ps = ps, beta_cor = beta_cor, p_cor = p_cor)
}


#' @title meta_analysis_hetero
#' @description Combine rep-by-rep pairwise TER tests between two groups of interest (e.g. two cell lines,
#' treatments, or batches) via meta-analysis, instead of testing the fully pooled data or splitting into
#' independent sub-datasets (c.f. \code{\link{validate_meta_test}}). Runs \code{\link{TER_all_pairs}} on every
#' replicate-vs-replicate pair among samples belonging to \code{group_A} or \code{group_B}, keeps only the pairs
#' that compare the two groups to each other (\code{pair_type == "hetero"}), and combines them with
#' \code{\link{meta_test}} using \code{comp_sym = TRUE} (compound symmetry is appropriate since these are all
#' repeated measurements of the same underlying contrast). Automates the workflow described in Ribolog's
#' empirical null / meta-analysis vignette (module 5, "merging rep-by-rep test results").
#' @param x A data matrix of RNA and RPF read counts (transcript ID column excluded), as used by \code{\link{logit_seq}}.
#' @param design Design matrix of the experiment describing samples and their attributes. Must contain \code{groupID}
#' and \code{uniqueID}.
#' @param uniqueID A variable (column) of \code{design} defining unique experimental preparations, from each
#' of which one RNA sample and one RPF sample was derived (usually a replicate name).
#' @param groupID A variable (column) of \code{design} defining the groups being compared (e.g. cell line,
#' treatment, batch). Passed to \code{\link{TER_all_pairs}}, which uses it to classify each replicate pair as
#' "homo" (same group) or "hetero" (different groups).
#' @param group_A,group_B The two levels of \code{groupID} to compare. Samples outside these two levels are
#' excluded before pairwise testing.
#' @param outcome The variable determining whether a vector of read counts is RNA or RPF. Default: "read_type".
#' @param adj_method P-value adjustment method passed to the underlying \code{\link{logit_seq}} calls (via
#' \code{\link{TER_all_pairs}}). Default: "none".
#' @param feature_list (Optional but strongly recommended) A vector containing IDs of genes/transcripts,
#' one per row of \code{x}. \code{\link{logit_seq}} (called once per replicate pair, via
#' \code{\link{TER_all_pairs}}) drops any transcript with 0 counts across all RPF (or RNA) samples
#' \emph{within that specific pair}, and different pairs can drop different transcripts. Supplying
#' \code{feature_list} lets \code{meta_analysis_hetero} track transcript identity through each pair and
#' correctly align them before combining; without it, rows are matched up positionally, which is only
#' correct if no pair drops any transcripts.
#' @return A list with two elements: \code{rep_by_rep_meta}, the output of \code{\link{meta_test}}
#' restricted to the transcripts that survived filtering in every rep-by-rep pair (a data frame with
#' \code{meta_beta} and \code{meta_p} columns, one row per surviving feature; a warning reports how many,
#' if any, transcripts were dropped for this reason); and \code{pooled_fit_pair}, the single pooled
#' \code{\link{logit_seq}} fit of \code{outcome ~ groupID} on the same \code{group_A}/\code{group_B}
#' samples - the "gold standard" that \code{rep_by_rep_meta} is meant to approximate (c.f.
#' \code{\link{validate_meta_test}}).
#' @examples
#' result <- meta_analysis_hetero(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, uniqueID = "replicate_name",
#'                                 groupID = "cell_line", group_A = "CN34", group_B = "LM1a",
#'                                 feature_list = as.vector(rr.v2_dummy$transcript))
#' # rep_by_rep_meta only contains transcripts common to every rep-by-rep pair, so align by name first:
#' common <- intersect(rownames(result$pooled_fit_pair), rownames(result$rep_by_rep_meta))
#' cor.test(result$pooled_fit_pair[common, "logFC_cell_lineLM1a"],
#'          result$rep_by_rep_meta[common, "meta_beta"], method = "spearman")
#' @export
meta_analysis_hetero <- function(x, design, uniqueID, groupID, group_A, group_B,
                                  outcome = "read_type", adj_method = "none", feature_list = NULL){

  pair_mask <- design[[groupID]] %in% c(group_A, group_B)
  pair_design <- droplevels(design[pair_mask, ])
  pair_x <- x[, pair_mask]
  if (!is.null(feature_list)) {
    rownames(pair_x) <- feature_list
  }

  pooled_model <- as.formula(paste(outcome, groupID, sep = "~"))
  pooled_fit_pair <- Ribolog::logit_seq(pair_x, pair_design, pooled_model,
                                         feature_list = feature_list, adj_method = adj_method)

  pair_split <- Ribolog::partition_to_uniques(x = pair_x, design = pair_design, uniqueID = uniqueID)

  rep_by_rep_fits <- Ribolog::TER_all_pairs(x = pair_split, design = pair_design, outcome = outcome,
                                             uniqueID = uniqueID, groupID = groupID, adj_method = adj_method)

  hetero_fits <- Filter(function(p) p$pair_type == "hetero", rep_by_rep_fits)

  extract_term_block <- function(fit_df){
    stats <- c("logFC", "std", "zvalue", "pvalue")
    cn <- colnames(fit_df)
    term <- setdiff(unique(sub(paste0("^(", paste(stats, collapse = "|"), ")_"), "", cn)), "intercept")
    fit_df[, paste0(stats, "_", term)]
  }
  hetero_fit_blocks <- lapply(hetero_fits, function(p) extract_term_block(p[["fit"]]))

  # Each pair's fit can have silently dropped a different subset of transcripts (e.g. 0 counts across
  # RPF samples within that specific pair), so align all blocks to the transcripts common to every pair
  # before combining - meta_test assumes identical row order across its input list.
  common_features <- Reduce(intersect, lapply(hetero_fit_blocks, rownames))
  if (length(common_features) == 0) {
    stop("No transcripts survived filtering in every rep-by-rep pair; cannot combine via meta_test.")
  }
  all_features <- unique(unlist(lapply(hetero_fit_blocks, rownames)))
  n_dropped <- length(all_features) - length(common_features)
  if (n_dropped > 0) {
    warning(sprintf(
      "%s of %s transcripts were dropped from the meta-analysis because they were filtered out (e.g. 0 counts in one read type) in at least one rep-by-rep pair.",
      n_dropped, length(all_features)))
  }
  hetero_fit_blocks <- lapply(hetero_fit_blocks, function(b) b[common_features, , drop = FALSE])

  rep_by_rep_meta <- Ribolog::meta_test(hetero_fit_blocks, comp_sym = TRUE, feature_list = common_features)

  list(rep_by_rep_meta = rep_by_rep_meta, pooled_fit_pair = pooled_fit_pair)
}


#' @title visualize_meta_analysis_hetero
#' @description Plot the two comparisons produced by \code{\link{meta_analysis_hetero}} (effect size and
#' p-value agreement between the pooled fit and the rep-by-rep meta-analysis) side by side as a pair of
#' themed \code{ggplot2} scatter plots with a bold y = x reference line, instead of two separate base R
#' plots.
#' @param result Output list of \code{\link{meta_analysis_hetero}}, containing \code{pooled_fit_pair} and
#' \code{rep_by_rep_meta}. Only transcripts common to both (i.e. present in \code{rep_by_rep_meta}) are
#' plotted.
#' @param beta_title,p_title Titles for the effect size and p-value panels, in bold.
#' Default: "Effect Size Agreement" / "P-value Agreement".
#' @param subtitle Subtitle shown under both panel titles (e.g. a dataset name). Default: \code{NULL}.
#' @param xlab_beta,ylab_beta Axis labels for the effect size panel. Default: "Pooled log TER" / "Meta log TER".
#' @param xlab_p,ylab_p Axis labels for the p-value panel. Default: "Pooled TER p-value" / "Meta TER p-value".
#' @param title_size Font size for both titles. Default: 22.
#' @param axis_title_size Font size for the (bold) axis titles. Default: 18.
#' @param axis_text_size Font size for the (bold) axis tick labels. Default: 14.
#' @param point_color Colour of the transcript points. Default: "#1F3864".
#' @param point_size,point_alpha Size/transparency of the points. Defaults: 0.6, 0.5.
#' @param ref_color Colour of the y = x reference line. Default: "#C0392B".
#' @param ref_linewidth Line width of the y = x reference line. Default: 1.1 (thicker than ggplot2's
#' default so it stays visible against a dense point cloud, e.g. near the origin of the p-value panel).
#' @param legend_position Where to place the legend (identical for both panels). Default: "top".
#' @param legend_text_size Font size for the (bold) legend labels. Default: 14.
#' @param show_cor Whether to compute and display a correlation coefficient in each panel's subtitle
#' (e.g. "Spearman rho = 0.93"). Default: TRUE.
#' @param cor_method Correlation method used for \code{show_cor} and the returned \code{cor.test} objects:
#' "pearson", "spearman", or "kendall". Default: "spearman".
#' @return Invisibly, a list with the two \code{ggplot2} objects (\code{beta_plot}, \code{p_plot}) and the
#' underlying \code{cor.test} results (\code{beta_cor}, \code{p_cor}). Called primarily for its side effect
#' of drawing both plots side by side in the current graphics device via \code{gridExtra::grid.arrange}.
#' @examples
#' result <- Ribolog::meta_analysis_hetero(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, uniqueID = "replicate_name",
#'                                          groupID = "cell_line", group_A = "CN34", group_B = "LM1a",
#'                                          feature_list = as.vector(rr.v2_dummy$transcript))
#' Ribolog::visualize_meta_analysis_hetero(result)
#' Ribolog::visualize_meta_analysis_hetero(result, subtitle = "LMCN dataset",
#'                                          title_size = 24, axis_title_size = 20)
#' @export
visualize_meta_analysis_hetero <- function(result,
                                            beta_title = "Effect Size Agreement", p_title = "P-value Agreement",
                                            subtitle = NULL,
                                            xlab_beta = "Pooled log TER", ylab_beta = "Meta log TER",
                                            xlab_p = "Pooled TER p-value", ylab_p = "Meta TER p-value",
                                            title_size = 22, axis_title_size = 18, axis_text_size = 14,
                                            point_color = "#1F3864", point_size = 0.6, point_alpha = 0.5,
                                            ref_color = "#C0392B", ref_linewidth = 1.1,
                                            legend_position = "top", legend_text_size = 14,
                                            show_cor = TRUE, cor_method = "spearman"){

  pooled_fit_pair <- result$pooled_fit_pair
  rep_by_rep_meta <- result$rep_by_rep_meta

  term_col <- setdiff(grep("^logFC_", colnames(pooled_fit_pair), value = TRUE), "logFC_intercept")
  if (length(term_col) != 1) {
    stop(paste0("Expected exactly one non-intercept logFC_ column in pooled_fit_pair, found ",
                length(term_col), ": ", paste(term_col, collapse = ", ")))
  }
  pval_col <- sub("^logFC_", "pvalue_", term_col)

  common <- intersect(rownames(pooled_fit_pair), rownames(rep_by_rep_meta))
  if (length(common) == 0) {
    stop("No transcripts are common to pooled_fit_pair and rep_by_rep_meta; nothing to plot.")
  }

  make_panel <- function(x, y, xlab, ylab, panel_title, panel_cor){
    plot_df <- data.frame(x = x, y = y)
    panel_subtitle <- subtitle
    if (show_cor) {
      cor_label <- paste0(toupper(substring(cor_method, 1, 1)), substring(cor_method, 2))
      cor_line <- sprintf("%s rho = %.3f", cor_label, unname(panel_cor$estimate))
      panel_subtitle <- if (is.null(subtitle)) cor_line else paste(subtitle, cor_line, sep = " | ")
    }
    ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(ggplot2::aes(colour = "Transcripts"), size = point_size, alpha = point_alpha) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1, linetype = "y = x"),
                            colour = ref_color, linewidth = ref_linewidth) +
      ggplot2::scale_colour_manual(name = NULL, values = c("Transcripts" = point_color)) +
      ggplot2::scale_linetype_manual(name = NULL, values = c("y = x" = "dashed")) +
      ggplot2::guides(colour = ggplot2::guide_legend(order = 1,
                        override.aes = list(size = 3.5, alpha = 1)),
                      linetype = ggplot2::guide_legend(order = 2,
                        override.aes = list(linewidth = 1.3))) +
      ggplot2::labs(title = panel_title, subtitle = panel_subtitle, x = xlab, y = ylab) +
      ggplot2::theme_minimal(base_size = axis_title_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5,
                                            margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_text(size = title_size * 0.6, colour = "grey40", hjust = 0.5,
                                               margin = ggplot2::margin(b = 12)),
        plot.title.position = "plot",
        axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
        axis.text = ggplot2::element_text(size = axis_text_size, colour = "grey20", face = "bold"),
        axis.line = ggplot2::element_line(colour = "grey50"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey92"),
        legend.position = legend_position,
        legend.box = "horizontal",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = legend_text_size, face = "bold"),
        legend.key = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(1.4, "lines")
      )
  }

  beta_x <- pooled_fit_pair[common, term_col]
  beta_y <- rep_by_rep_meta[common, "meta_beta"]
  p_x <- pooled_fit_pair[common, pval_col]
  p_y <- rep_by_rep_meta[common, "meta_p"]

  beta_cor <- stats::cor.test(beta_x, beta_y, method = cor_method)
  p_cor <- stats::cor.test(p_x, p_y, method = cor_method)

  beta_plot <- make_panel(beta_x, beta_y, xlab_beta, ylab_beta, beta_title, beta_cor)
  p_plot <- make_panel(p_x, p_y, xlab_p, ylab_p, p_title, p_cor)

  gridExtra::grid.arrange(beta_plot, p_plot, ncol = 2)

  invisible(list(beta_plot = beta_plot, p_plot = p_plot, beta_cor = beta_cor, p_cor = p_cor))
}


#' @title visualize_meta_analysis
#' @description Plot the two correlation matrices produced by \code{\link{validate_meta_test}} (effect
#' size and p-value agreement between the pooled fit, each split-level fit, and the meta-analysis
#' combination) side by side as a pair of \code{corrplot::corrplot} heatmaps, instead of one above the other.
#' @param validation Output list of \code{\link{validate_meta_test}}, containing \code{beta_cor} and
#' \code{p_cor} (each an \code{Hmisc::rcorr} object).
#' @param beta_title Title for the effect size (logFC) correlogram, in bold. Default: "Effect Size Correlation".
#' @param p_title Title for the p-value correlogram, in bold. Default: "P-value Correlation".
#' @param title_size Font size (cex) for both titles. Default: 1.4.
#' @param label_color Colour of the row/column labels. Default: "black" (corrplot's own default is red).
#' @param label_size Font size (cex) for the row/column labels. Default: 1.3 (corrplot's own default, 1,
#' is quite small).
#' @param coef_color Colour of the correlation values written inside each cell. Default: "white".
#' @param coef_size Font size (cex) for those values. Default: 1.
#' @param method Cell rendering method passed to \code{corrplot::corrplot} (e.g. "color", "circle",
#' "number"). Default: "color".
#' @param col Colour palette passed to \code{corrplot::corrplot}. Default: \code{NULL} (corrplot's own
#' diverging blue-white-red palette).
#' @return Invisibly \code{NULL}; called for its side effect of drawing the two corrplots in the current
#' graphics device.
#' @details
#' Uses \code{graphics::par(mfrow = c(1, 2))} to place the two plots side by side (restoring the previous
#' \code{par()} settings afterward, so it doesn't affect later plots). The colour scale legend is only
#' drawn once, on the first (effect size) plot, since it's identical for both.
#' @examples
#' validation <- Ribolog::validate_meta_test(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                                            split_var = "cell_line_origin", term_name = "lung_metastasisY",
#'                                            feature_list = as.vector(rr.v2_dummy$transcript))
#' Ribolog::visualize_meta_analysis(validation)
#' Ribolog::visualize_meta_analysis(validation, beta_title = "logFC agreement",
#'                                   p_title = "P-value agreement", label_size = 1.5, title_size = 1.6)
#' @export

visualize_meta_analysis <- function(validation, beta_title = "Effect Size Correlation",
                                     p_title = "P-value Correlation", title_size = 1.4,
                                     label_color = "black", label_size = 1.3,
                                     coef_color = "white", coef_size = 1,
                                     method = "color", col = NULL){

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mfrow = c(1, 2))

  corrplot::corrplot(validation$beta_cor$r, method = method, col = col,
                      addCoef.col = coef_color, number.cex = coef_size,
                      tl.col = label_color, tl.cex = label_size,
                      cl.pos = "r", mar = c(0, 0, 3, 0))
  graphics::title(main = beta_title, font.main = 2, cex.main = title_size)

  corrplot::corrplot(validation$p_cor$r, method = method, col = col,
                      addCoef.col = coef_color, number.cex = coef_size,
                      tl.col = label_color, tl.cex = label_size,
                      cl.pos = "n", mar = c(0, 0, 3, 0))
  graphics::title(main = p_title, font.main = 2, cex.main = title_size)

  invisible(NULL)
}


#' @title visualise_empirical_null
#' @description Function to visualise the empirical null distribution
#' @param enz_dataframe A vector (or single-column data frame) of the empirical null distribution
#' generated by \code{\link{generate_ENZ}}/\code{\link{build_empirical_null}}.
#' @param plot_density A boolean to decide if the plot includes a density curve overlay or just a histogram.
#' @param title Plot title. Default: "ENZ distribution".
#' @param subtitle Plot subtitle. Default: "LMCN data".
#' @param xlab X-axis label. Default: "Empirical null Z".
#' @param ylab Y-axis label. Default: \code{NULL}, which lets ggplot2 auto-label it ("density" or "count"
#' depending on \code{plot_density}).
#' @param title_size Font size for the title. Default: 22.
#' @param axis_title_size Font size for the axis titles (\code{xlab}/\code{ylab}). Default: 16.
#' @param axis_text_size Font size for the axis tick labels. Default: 12.
#' @param binwidth Histogram bin width. Default: 0.25.
#' @param hist_fill Fill colour for the histogram bars. Default: "#5B9BD5".
#' @param density_color Line colour for the density curve (only used when \code{plot_density = TRUE}).
#' Default: "#1F3864".
#' @param legend_position Where to place the legend distinguishing the histogram from the density curve
#' (only shown when \code{plot_density = TRUE}, since a legend isn't useful for a single series).
#' Default: "top".
#' @param xlim (Optional) A length-2 numeric vector giving the x-axis display range manually, e.g.
#' \code{c(-5, 5)}. Zooms the view only (via \code{ggplot2::coord_cartesian}) - the histogram and density
#' curve are still computed from the full data, so this does not distort them or drop any values (unlike
#' \code{ggplot2::xlim}/\code{scale_x_continuous(limits=)}, which would). Overrides \code{auto_xlim}.
#' @param auto_xlim If \code{xlim} is not given, automatically crop the view to the central \code{cl}\%
#' of the data (trimming sparse outlier tails), the same idea used by \code{\link{rlength_distr_rW}}'s
#' \code{cl} argument. Default: TRUE. Set to FALSE to show the full range instead.
#' @param cl Central percentage of the data to keep in view when \code{auto_xlim = TRUE} and \code{xlim}
#' is not given. Default: 99.9.
#' @return A ggplot2 object of the histogram
#' @examples
#' rr.v2_dummy.enz <- Ribolog::build_empirical_null(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy,
#'                                             uniqueID = "replicate_name", groupID = "cell_line")
#' Ribolog::visualise_empirical_null(rr.v2_dummy.enz)
#' Ribolog::visualise_empirical_null(rr.v2_dummy.enz, title = "Empirical null", subtitle = "My dataset",
#'                                    axis_text_size = 14, axis_title_size = 18)
#' Ribolog::visualise_empirical_null(rr.v2_dummy.enz, xlim = c(-5, 5))  # crop manually
#' Ribolog::visualise_empirical_null(rr.v2_dummy.enz, auto_xlim = FALSE)  # show full tails
#' @export

visualise_empirical_null <- function(enz_dataframe, plot_density = TRUE,
                                      title = "ENZ distribution", subtitle = "LMCN data",
                                      xlab = "Empirical null Z", ylab = NULL,
                                      title_size = 22, axis_title_size = 16, axis_text_size = 12,
                                      binwidth = 0.25,
                                      hist_fill = "#5B9BD5", density_color = "#1F3864",
                                      legend_position = "top",
                                      xlim = NULL, auto_xlim = TRUE, cl = 99.9){

  enz_df <- data.frame(enz = as.numeric(unlist(enz_dataframe)))

  if (is.null(xlim) && auto_xlim) {
    trimmed <- unname(stats::quantile(enz_df$enz, c((1 - cl / 100) / 2, 1 - (1 - cl / 100) / 2), na.rm = TRUE))
    pad <- diff(trimmed) * 0.1
    xlim <- c(trimmed[1] - pad, trimmed[2] + pad)
  }

  p <- ggplot2::ggplot(enz_df, ggplot2::aes(x = enz))

  if (plot_density){
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density), fill = "Histogram"),
                               binwidth = binwidth, colour = "white", linewidth = 0.2, alpha = 0.9) +
      ggplot2::geom_density(ggplot2::aes(colour = "Kernel density"), fill = NA, linewidth = 1.1, bw = 0.5) +
      ggplot2::scale_fill_manual(name = NULL, values = c("Histogram" = hist_fill)) +
      ggplot2::scale_colour_manual(name = NULL, values = c("Kernel density" = density_color)) +
      ggplot2::guides(fill = ggplot2::guide_legend(order = 1), colour = ggplot2::guide_legend(order = 2))
  } else {
    p <- p +
      ggplot2::geom_histogram(binwidth = binwidth, colour = "white", fill = hist_fill, linewidth = 0.2, alpha = 0.9)
    legend_position <- "none"
  }

  p +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = axis_title_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5,
                                          margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = title_size * 0.6, colour = "grey40", hjust = 0.5,
                                             margin = ggplot2::margin(b = 12)),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text = ggplot2::element_text(size = axis_text_size, colour = "grey20"),
      axis.line = ggplot2::element_line(colour = "grey50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = axis_text_size),
      legend.key = ggplot2::element_blank()
    )
}


#' @title fit_empirical_null
#' @description Function to fit the empirical null distribution to a given standard distribution
#' @param enz_dataframe A dataframe of the empirical null distribution generated by generate_ENZ function
#' @param distribution A keyword for the distribution to be fitted to the empirical null
#' @return
#' Dataframe of statistical results from the fitting
#' @examples
#' rr.v2_dummy.enz <- Ribolog::build_empirical_null(rr.v2_dummy[, -1], Ribolog::sample_attributes_dummy,
#'                                             uniqueID = "replicate_name", groupID = "cell_line")
#' enz_fit <- Ribolog::fit_empirical_null(rr.v2_dummy.enz, 'norm')
#' @export

fit_empirical_null <- function(enz_dataframe, distribution){
  return(fitdistrplus::fitdist(enz_dataframe, distribution))
}



#' @title partition_to_uniques
#' @description Function to convert a RNA+RPF data frame to a sample-by-sample list.
#' @param x A data frame or matrix containing RNA+RPF count data where each row is a transcript and each column is RNA or RPF counts of one sample.
#' This object must contain only count data (and not, for example, a transcript ID column).
#' @param design Design matrix of the experiment describing samples and their attributes.
#' The i-th row in the design matrix describes the i-th column in the input data frame \code{x}.
#' @param uniqueID A variable (column) of the design matrix defining unique experimental preparations
#' from each of which one RNA sample and one RPF sample was derived. It corresponds to the highest resolution
#' (lowest level) of classification of samples in the data set apart from the RNA/RPF distinction
#' and is usually equal to replicate name in biological experiments.
#' @return A list where each element is a data frame containing RNA and RPF count of one replicate and its attributes from the design matrix.
#' @examples
#' rr.v2_dummy.split <- partition_to_uniques(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, "replicate_name")
#' # The first column of rr.v2_dummy contained transcript IDs and was thus excluded from input.
#' @export
partition_to_uniques <- function(x, design, uniqueID){
  xt <- t(x)
  xtd <- cbind(design, xt)
  xtdl <- split(xtd, xtd[,uniqueID])
  return(xtdl)
}


#' @title TER_all_pairs
#' @description Function to perform the logit TER test between all pairs of samples in a data set.
#' @param x A sample-by-sample list of RNA and RPF count data and sample attributes produced by \code{\link{partition_to_uniques}}.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' @param outcome The variable determining whether a vector of read counts is RNA or RPF.
#' This is usually the name of the response variable in the TER test logistic regression performed through \code{\link{logit_seq}}. Default: \code{"read_type"}.
#' @param uniqueID A variable (column) of the design matrix defining unique experimental preparations
#' from each of which one RNA sample and one RPF sample was derived. It corresponds to the highest resolution
#' (lowest level) of classification of samples in the data set apart from the RNA/RPF distinction
#' and is usually equal to replicate name in biological experiments.
#' @param groupID A variable (column) of the design matrix indicating which replicates should be grouped together.
#' All experimental units having the same \code{groupID} will be considered replicates of the same biological sample
#' (or members of the same group of samples).
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return
#' A list of lists containig the results of all pairwise TER tests. If there are n samples in the input list, the output list will consist of C(n,2) elements.
#' Each element of the list is in turn a list with four attributes:
#' - \code{uniqueID}s of the two samples compared
#' - \code{groupID}s of the two samples compared
#' - \code{pair_type} (\code{"homo"} if the two \code{groupID}s are equal and \code{"hetero"} otherwise)
#' - \code{fit} containing the output of the TER test in a data frame. See \code{\link{logit_seq}} for details.
#' @examples
#' # Subset to a handful of transcripts so the pairwise fits run quickly.
#' rr.v2_dummy.split <- partition_to_uniques(rr.v2_dummy[1:200, -1], Ribolog::sample_attributes_dummy, "replicate_name")
#' rr.v2_dummy.pairwise <- TER_all_pairs(rr.v2_dummy.split, Ribolog::sample_attributes_dummy, "read_type", "replicate_name",
#'                                  "cell_line", adj_method = "none")
#' @export
TER_all_pairs <- function(x, design, outcome = "read_type", uniqueID, groupID, adj_method){
  pair_results <- list()
  n <- length(x)
  n_design_cols <- dim(design)[2]

  for (i in c(2:n)){
    for (j in c(1:(i-1))){
      list_ij <- list()
      list_ij[["uniqueIDs"]] <- sort(c(names(x)[i], names(x)[j]))

      list_ij[["groupIDs"]] <- sort(c(as.character(x[[i]][, groupID][1]), as.character(x[[j]][, groupID][1])))
      if (identical(list_ij[["groupIDs"]][1], list_ij[["groupIDs"]][2])) list_ij[["pair_type"]] = "homo" else list_ij[["pair_type"]] = "hetero"

      model1 <- as.formula(paste(as.factor(outcome), as.factor(uniqueID), sep = "~"))
      data_ij <- rbind(x[[i]], x[[j]])[, -c(1:n_design_cols)]
      design_ij <- rbind(x[[i]], x[[j]])[, c(1:n_design_cols)]
      list_ij[["fit"]] <- Ribolog::logit_seq(t(data_ij), design_ij, model1, adj_method=adj_method, long_output = TRUE)
      name_ij <- paste(list_ij[["uniqueIDs"]], collapse = "_vs_")
      pair_results[[name_ij]] <- list_ij

    }
  }
  return(pair_results)
}