#' @import data.table
#' @import Biostrings
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import robustbase
#' @import qvalue
#' @import nortest
#' @import matrixStats
#' @import sm
#' @import corrplot
#' @import DescTools
#' @import GenomicAlignments
#' @import rlist
#' @import gdata
#' @import nlme
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
#' rr_LMCN.v2.enz <- generate_ENZ(x = rr_LMCN.v2.split, design = sample_attributes_LMCN, outcome = "read_type", uniqueID = "replicate_name", groupID = "cell_line")
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
#' @examples
#' fit1_LMCN_ENZ_p <- test_ENZ(x = fit1_LMCN, enz = rr_LMCN.v2.enz, zcols = 7)
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
  colnames(fitz_enz_p) <- gsub("z value", "ep", colnames(fitz_enz_p))
  if (keep_data == TRUE){
    output <- data.frame(x, fitz_enz_p)
  } else if (keep_data == FALSE) {
    output <- fitz_enz_p
  }
  return(output)
}


#' @title com_cor_p
#' @description Function to combine correlated p-values
#' @param p Array (nxm matrix) of p-values. n: number of features e.g. genes or transcripts, m: number of tests being combined.
#' @param w Array (nxm matrix) of weights. n: number of features e.g. genes or transcripts, m: number of tests being combined.
#' @param cm Correlation matrix of z-values or any other test statistics of individual tests that are being combined
#' (mxm). m: number of samples.
#' @details
#' This is based on the the thoery expounded in a 2003 paper by Kepher Makambi.
#' @references
#' Makambi, K. 2003. Weighted inverse chi-square method for correlated significance tests.
#' Journal of Applied Statistics, 30(2): 225-234.
#' @return
#' A matrix with three columns: [1] M (combined test statistic) [2] nu (degree of freedom of chi-sq distribution)
#' [3] p.com (the combined or meta-analytical p-value)

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
#' @param wt_effects Indicator variable. TRUE: Combine effect sizes weighted by their inverse of variance.
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
#' The effect of lung metastasis is tested in two genetic backgrounds sperately, then combined:
#' fit1a_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,c(2:5,10:13)], sample_attributes_LMCN[c(1:4,9:12), ], read_type ~ lung_metastasis, as.vector(rr_LMCN.v2$transcript))
#' fit1b_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,c(6:9,14:17)], sample_attributes_LMCN[c(5:8,13:16), ], read_type ~ lung_metastasis, as.vector(rr_LMCN.v2$transcript))
#' We are not interested in the intercept terms, so keep only the last four columns of the outputs which pertain to lung metastasis:
#' fit1ab_LMCN_list <- list("fit1a_LMCN" = fit1a_LMCN[, c(5:8)], "fit1b_LMCN" = fit1b_LMCN[, c(5:8)])
#' And now, combine:
#' fit1ab_LMCN_meta <- meta_test(fit1ab_LMCN_list, feature_list = rr_LMCN.v2$transcript)
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
    fit_cs_z <- gls(z ~ 1, data = z_df2lg, corr = corCompSymm(, form= ~ 1 | transcript))
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
#' @value
#' An augmented matrix where SD(beta) and z columns are added between beta and p for each predictor (and intercept, if applicable).
#' @examples
#' fit1_LMCN_short <- fit1_LMCN[, c(TRUE, FALSE, FALSE, TRUE)]
#' fit1_LMCN_restored <- add_sd_z(fit1_LMCN_short)
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



#' @title visualise_empirical_null
#' @description Function to visualise the empirical null distribution
#' @param enz_dataframe A dataframe of the empirical null distribution generated by generate_ENZ function
#' @param plot_density A boolean to decide if the plot includes a density plot or just a histogram
#' @value
#' A ggplot2 object of the histogram
#' @examples
#' stats <- Ribolog::visualise_empirical_null(rr_LMCN.v2.enz)
#' @export

visualise_empirical_null <- function(enz_dataframe, plot_density=TRUE){

  if (plot_density){
    ggplot(as.data.frame(enz_dataframe)) + aes(x=enz_dataframe,y=..density..) +
    geom_histogram(binwidth=0.25, colour="black", fill="white") +
    geom_density(aes(x=enz_dataframe), color="darkblue", fill="lightblue", alpha=0.5, bw=0.5) +
    ggtitle("ENZ distribution \n LMCN data") + xlab("Empirical null Z") +
    theme(plot.title = element_text(size=22, hjust = 0.5), axis.title=element_text(size=16, hjust=0.5))
  } else {
    ggplot(as.data.frame(enz_dataframe)) + aes(x=enz_dataframe) +
    geom_histogram(binwidth=0.25, colour="black", fill="white") +
    ggtitle("ENZ distribution \n LMCN data") + xlab("Empirical null Z") +
    theme(plot.title = element_text(size=22, hjust = 0.5), axis.title=element_text(size=16, hjust=0.5))
  }
}



#' @title fit_empirical_null
#' @description Function to fit the empirical null distribution to a given standard distribution
#' @param enz_dataframe A dataframe of the empirical null distribution generated by generate_ENZ function
#' @param distribution A keyword for the distribution to be fitted to the empirical null
#' @value
#' Dataframe of statistical results from the fitting
#' @examples
#' stats <- Ribolog::fit_empirical_null(rr_LMCN.v2.enz, 'norm')
#' @export

fit_empirical_null <- function(enz_dataframe, distribution){
  return(fitdistrplus::fitdist(enz_dataframe, distribution))
}
