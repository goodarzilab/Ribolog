#' @import data.table
#' @import Biostrings
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import robustbase
#' @import qvalue
#' @import nortests
#' @import fitdistrplus
#' @import matrixStats
#' @import sm
#' @import epiR
#' @import corrplot
#' @import mvmeta
#' @import DescTools
#' @import GenomicAlignments
#' @import rlists
#' @import gdata
#' @import nlme



#' @title logit_seq
#' @description Function to perform the logistic regression test for differential translational efficiency
#' @param x Input data frame where each column contains RNA or RPF counts from a single sample.
#' Rows are genes/transcripts.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' i-th row in the design matrix corresponds to the i-th column in the input data frame.
#' @param model Regression equation modeling the odds ratio of RPF/RNA counts against the selected design variables (sample attributes)
#' @param feature_list (Optional) A vector containing IDs of genes/transcripts.
#' Must have the same length as the row number of input data frame.
#' @return A matrix containing the output of the regression. If long output is requested, four
#' values are reported for each predictor in the regression 'model' including the intercept:
#' regression coefficient (beta), standard deviation of the estimated beta, z-score and Wald-test p-value.
#' If long output is not specified, only regression coefficients and p-values are reported.
#' Gene/transcript IDs provided by the 'feature_list' argument is added to the output matrix as row names.
#' @details The response variable is always the read type which should be a factor variable with two levels: "RNA" and "RPF".
#' Translational efficiency of each sample if formulated as the odds of RPF vs. RNA reads.
#' The change in translational efficiency between samples or based on unit values of any predictor
#' (design variable) is given by the exponentiated regression coefficient.
#' Exponentiated intercept gives the TE for the sample with all the attributes at the reference level.
#' If the predcitor variable is categorical, its levels are sorted alphabetically, the first level is set to reference
#' and all other levels are compared to it. The reference level can be manually changed using the \code{\link{relevel()}} function of R.
#' @examples
#' Test the effect of lung metastasis on translational efficiency:
#' fit1_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,-1], sample_attributes_LMCN, read_type ~ lung_metastasis, as.vector(rr_LMCN.v2$transcript))
#' Test the effects of lung metastasis and cell line origin on translational efficiency:
#' fit2_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,-1], sample_attributes_LMCN, read_type ~ lung_metastasis + cell_line_origin, as.vector(rna.rpf.combined.m5$transcript))
#' Test the effects of lung metastasis, cell line origin and their interaction on translational efficiency:
#' fit3_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,-1], sample_attributes_LMCN, read_type ~ lung_metastasis * cell_line_origin, as.vector(rna.rpf.combined.m5$transcript))
#' Test the effect of cell line on translational efficiency (cell line "CN34" is used as reference because it comes first alphabetically):
#' fit4_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,-1], sample_attributes_LMCN, read_type ~ cell_line, as.vector(rr_LMCN.v2$transcript))
#' Test the effect of cell line on translational efficiency with cell line "MDA" set as reference:
#' sample_attributes_LMCN$cell_line <- relevel(sample_attributes_LMCN$cell_line, ref = "MDA")
#' fit5_LMCN <- Ribolog::logit_seq(rr_LMCN.v2[,-1], sample_attributes_LMCN, read_type ~ cell_line, as.vector(rr_LMCN.v2$transcript))
#' @export

logit_seq <- function(x, design, model, feature_list=NULL, long_output = FALSE){
  logit_seq_gene <- function(m){
    prep <- data.frame(design,m)
    fit <- suppressWarnings(glm(model, data=prep, family="binomial"(link="logit"), weights = m))
    sfit <- summary(fit)
    return(c(t(sfit$coefficients)))
  }
  logit_x <- t(apply(x, 1, logit_seq_gene))
  prep1 <- cbind.data.frame(design, data.frame(counts1 = as.numeric(x[1,])))
  fit1 <- suppressWarnings(glm(model, data=prep1, family="binomial"(link="logit"), weights = counts1))
  sfit1 <- summary(fit1)

  colnames(logit_x) <- apply(expand.grid(colnames(sfit1$coefficients), rownames(sfit1$coefficients)), 1, paste, collapse="_")
  rownames(logit_x) <- feature_list
  if (long_output == FALSE){
    logit_x <- logit_x[,c(TRUE, FALSE, FALSE, TRUE)]
  }
  return(logit_x)
}



#' @title adj_TER_p
#' @description Function to adjust p-values from the TER tests
#' @param x Output object of a logit_seq function
#' @param pcols A vector specifying the positions of unadjusted p-value columns in x
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return A data frame containing all the original information with the adjusted p-values
#' appended to the end.
#' @details P-values from differential translational efficiency test on individual genes/transcripts
#' are correctd for multiple testing.
#' @examples
#' fit1_LMCN_qval <- adj_TER_p(fit1_LMCN, c(2,4), "qvalue")
#' fit3_LMCN_fdr <- adj_TER_p(fit3_LMCN, c(2,4,6,8), "fdr")
#' @export

adj_TER_p <- function(x, pcols, adj_method){
  y <- x[, pcols, drop = FALSE]
  if (adj_method == "qvalue"){
    y <- apply(y, 2, function(t) qvalue::qvalue(t)$qvalues)
  } else {
    y <- apply(y, 2, function(t) p.adjust(t, method = adj_method))
  }

  newnames <- paste0(adj_method, "_", colnames(y))
  z <- data.frame(x,y)
  colnames(z)[(NCOL(x)+1) : NCOL(z)] <- newnames
  return(z)
}
