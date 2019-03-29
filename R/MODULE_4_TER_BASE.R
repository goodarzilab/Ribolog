#' @title logit_seq
#' @description Function to perform the logistic regression test for differential translational efficiency
#' @param x Input data frame where each column contains RNA or RPF counts from a single sample.
#' Rows are genes/transcripts.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' i-th row in the design matrix corresponds to the i-th column in the input data frame.
#' @param model Regression equation modeling the odds ratio of RPF/RNA counts against the selected design variables (sample attributes)
#' @param n Number of estimated parameters (including intercept)
#' @param feature.list (Optional) A vector containing IDs of genes/transcripts.
#' Must have the same length as the row number of input data frame.
#' @return A matrix containing the output of the regression.
#' Four values are reported for each predictor in the regression 'model' including the intercept:
#' regression coefficient (beta), standard deviation of the estimated beta, z-score and Wald-test p-value.
#' Gene/transcript IDs provided by the 'feature.list' argument is added to the output matrix as row names.
#' @details The response variable is always the read type which should be a factor variable with two levels: "RNA" and "RPF".
#' Translational efficiency of each sample if formulated as the odds of RPF vs. RNA reads.
#' The change in translational efficiency between samples or based on unit values of any predictor
#' (design variable) is given by the exponentiated regression coefficient.
#' Exponentiated intercept gives the TE for the sample with all the attributes at the reference level.
#' @examples
#' fit.1 <- logit_seq(rna.rpf.combined.m5[,-1], sample.attributes, read.type ~ lung.metastasis, 2, as.vector(rna.rpf.combined.m5$transcript))
#' The n=2 parameters include the intercept and the effect of lung.metastasis.
#' fit.2 <- logit_seq(rna.rpf.combined.m5[,-1], sample.attributes, read.type ~ lung.metastasis*cell.line.origin, 4, as.vector(rna.rpf.combined.m5$transcript))
#' The n=4 parameters include the intercept, the main effects of lung.metastasis and cell.line.origin and their interaction.
#' @export
logit_seq <- function(x, design, model, n, feature.list=NULL){
  logit_seq_gene <- function(m){
    prep <- data.frame(design,m)
    fit <- glm(model, data=prep, family="binomial"(link="logit"), weights = m)
    sfit <<- summary(fit)
    return(c(t(sfit$coefficients)))
  }
  logit_seq_gene_try <- function(m){
    out<-tryCatch(
      {
        result_logistic <- logit_seq_gene(m)
      },
      error = function(cond){
        return(rep(NA, times=4*n))
      }
    )
    return(out)
  }
  logit.x <- t(apply(x, 1, logit_seq_gene_try))
  colnames(logit.x) <- apply(expand.grid(colnames(sfit$coefficients), rownames(sfit$coefficients)), 1, paste, collapse=".")
  rownames(logit.x) <- feature.list
  return(logit.x)
}



#' @title adj_TER_p
#' @description Function to adjust p-values from the TER tests
#' @param x Output object of a logit_seq function
#' @param p.columns A vector specifying the positions of unadjusted p-value columns in x
#' @param adj.method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @return A data frame containing all the original information with the adjusted p-values
#' appended to the end.
#' @details P-values from differential translational efficiency test on individual genes/transcripts
#' are correctd for multiple testing.
#' @examples
#' fit.1.qval <- adj_TER_p(fit.1, c(4,8), "qvalue")
#' fit.2.fdr <- adj_TER_p(fit.2, c(4,8,12,16), "fdr")
#' @export

adj_TER_p <- function(x, pcols, adj.method){
  y <- x[,pcols]
  if (adj.method == "qvalue"){
    y <- apply(y, 2, function(t) qvalue(t)$qvalues)
  } else {
    y <- apply(y, 2, function(t) p.adjust(t, method = adj.method))
  }
  colnames(y) <- paste0(adj.method, ".", colnames(y))
  return(data.frame(x,y))
}

