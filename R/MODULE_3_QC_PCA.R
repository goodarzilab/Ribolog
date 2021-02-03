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



#' @title min_count_filter
#' @description Function to filter out genes with counts below a minimum in one or more samples.
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data.frame may contain additional columns for gene/transcript ID or other metadata.
#' @param mincount A single number (float), minimum RNA or RPF count required for a gene/transcript to be retained.
#' @param columns A vector specifying the columns to be considered for minimum count filtering.
#' @param method The method of filtering. Options: \code{"all", "average"}. Default: \code{"all"}.
#' @details If \code{method="all"} is chosen, a gene passes the filtering only if all samples specified
#' by the \code{columns} argument have values \code{>= mincount}. If \code{method="average"} is chosen, a gene
#' passes the filtering if the average count among the specified \code{columns} is \code{>= mincount}.
#'
#' Use the \code{columns} argument to exclude gene/transcript ID and other metadata columns from the
#' calculations or to filter RNA and RPF counts separately while keeping them in the same data.frame.
#' @return Filtered data frame containing all the original columns but only the rows (genes)
#' that pass the filtering criterion.
#' @examples
#' rr_LMCN.v1 <- min_count_filter(rr_LMCN, mincount = 5, columns = c(2:9), method = "all")
#' rr_LMCN.v2 <- min_count_filter(rr_LMCN.v1, mincount = 2, columns = c(10:17), method = "average")
#' @export
min_count_filter <- function(x, mincount, columns, method="all"){
  if (method=="all"){
    x <- x[!rowSums(x[,columns] < mincount),]
  } else if (method=="average"){
    x <- x[rowMeans(x[,columns]) >= mincount,]
  }
  return(x)
}



#' @title create_te
#' @description Function to create a TE (translational efficiency) data frame from a combined RNA+RPF data frame.
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data.frame may contain additional columns for gene/transcript ID or other metadata.
#' @param idcolumns A vector specifying the columns to be excluded from the calculations.
#' These columns may contain the gene/transcript ID or any metadata that need to be preserved.
#' @param rnacolumns A vector specifying the columns containing RNA counts
#' @param rpfcolumns A vector specifying the columns containing RPF counts
#' @param allow_zero_rpf A boolean allowing the total RPF counts of a transcript to be 0
#' @details Translational efficiency is calculated as RPF/RNA.
#' The number and order of samples must be the same in RNA and RPF columns.
#' @return A data frame containing the original ID columns and the calculated TE columns.
#' @examples
#' te_LMCN <- create_te(rr_LMCN.v2, 1, c(2:9), c(10:17))
#' @export
create_te <- function(x, idcolumns=NULL, rnacolumns, rpfcolumns, allow_zero_rpf=FALSE){

  total_rpf_counts <- as.data.frame(rowSums(x[,rpfcolumns]))
  total_rna_counts <- as.data.frame(rowSums(x[,rnacolumns]))

  rownames(total_rpf_counts) <- x$transcript
  rownames(total_rna_counts) <- x$transcript

  empty_rpf_transcripts <- c()
  empty_rna_transcripts <- c()

  for (transcript in rownames(total_rpf_counts)) {
    if (total_rpf_counts[transcript,] == 0 && !allow_zero_rpf) {
      empty_rpf_transcripts <- c(empty_rpf_transcripts, transcript)
    }

    if (total_rna_counts[transcript,] == 0) {
      empty_rna_transcripts <- c(empty_rna_transcripts, transcript)
    }
  }

  if (length(empty_rpf_transcripts) > 0 && !allow_zero_rpf){

    warning(sprintf('There are ( %s ) transcripts that have 0 counts across all the RPF samples.You can allow RPF counts
    to be 0 using allow_zero_rpf=TRUE option or filter them using Ribolog::min_count_filter. These transcripts have been removed for now.', length(empty_rpf_transcripts)))
    x <- Ribolog::min_count_filter(x, mincount = 5, columns = rpfcolumns, method = "all")
  }

  if (length(empty_rna_transcripts) > 0){

    warning(sprintf('There are ( %s ) transcripts that have 0 counts across all the RNA samples.You can filter them using Ribolog::min_count_filter.
    These transcripts have been removed for now.', length(empty_rna_transcripts)))
    x <- Ribolog::min_count_filter(x, mincount = 2, columns = rnacolumns, method = "average")
  }


  y <- data.frame(x[,idcolumns], x[,rpfcolumns]/x[,rnacolumns])
  names(y)[idcolumns] <- names(x)[idcolumns]
  names(y) <- gsub("rpf", "te", names(y))
  return(y)
}



#' @title row_center
#' @description Function to center a selected block of a data frame on its row means.
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param columns A vector specifying the columns to be included for row-centering.
#' @return A data frame where the specified columns from the input are row-centered and
#' the rest is intact.
#' @examples
#' te_LMCN.v2.cent <- row_center(te_LMCN.v2, columns = c(2:9))
#' @export
row_center <- function(x, columns){
  x[,columns] <- t(apply(x[,columns], 1, function(y) y-mean(y)))
  return(x)
}



#' @title row_standardize
#' @description Function to standardize a selected block of a data frame row-wise.
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param columns A vector specifying the columns to be included for row-wise standardization
#' @return A data frame where the specified columns from the input are row-standardized and
#' the rest is intact.
#' @details Row mean is subtracted from each element in the row and the result is divided by row standard deviation.
#' @examples
#' te_LMCN.v2.stnd <- row_standardize(te_LMCN.v2, columns = c(2:9))
#' @export
row_standardize <- function(x, columns){
  x[,columns] <- t(apply(x[,columns], 1, function(y) (y-mean(y))/sd(y)))
  return(x)
}



#' @title pca_qc
#' @description Function to produce PCA output for data visualization and QC purposes.
#' @param x A numeric matrix containing RNA, RPF, TE or some other type of data.
#' Rows are genes/transcripts and columns are samples.
#' Gene/transcript IDs and sample names may be given as row names and column names,
#' respectively (but not as additional non-numeric columns or rows).
#' @param n The number of principal components to be plotted
#' @param outfile The path and name of the output pdf file containing the PCA plots (optional). Default: NULL.
#' @param ID A vector containing group IDs to color-code the samples on the PCA plot
#' (must correspond to the order of samples in \code{x}) (optional). Default: NULL.
#' @details A summary of the PC analysis is printed out to standard output.
#' If outfile is not specified, PCA plots will be printed to standard output i.e. the plots panel in Rstudio.
#' @examples
#' pca_qc(te_LMCN.v2.stnd[,-1], n = 4)
#' The first column of input data (transcript ID) had to be removed to create a full-numeric input dataset.
#' @export
pca_qc <- function(x, n, outfile = NULL, ID = NULL){
  x.pca <- prcomp(x[!rowSums(is.na(x)),])
  print(summary(x.pca))
  rotation.x.pca <- data.frame(x.pca$rotation)
  var.x.pca <- summary(x.pca)$importance[2,]

  if (is.null(outfile)){
    for (i in 1:dim(combn(n,2))[2]){
      x_i <- rotation.x.pca[,combn(n,2)[1,i]]
      y_i <- rotation.x.pca[,combn(n,2)[2,i]]
      margin_x <- (max(x_i) - min(x_i)) * 0.2
      margin_y <- (max(y_i) - min(y_i)) * 0.2
      print(ggplot(rotation.x.pca, aes(x_i, y_i, color = ID)) +
              geom_point(shape = 16) +
              geom_label_repel(aes(label=colnames(x)))+
              xlim(min(x_i)-margin_x, max(x_i)+margin_x)+
              ylim(min(y_i)-margin_y,max(y_i)+margin_y)+
              xlab(paste0(colnames(rotation.x.pca)[combn(n,2)[1,i]]," (", var.x.pca[combn(n,2)[1,i]]*100,"% of variance)"))+
              ylab(paste0(colnames(rotation.x.pca)[combn(n,2)[2,i]]," (", var.x.pca[combn(n,2)[2,i]]*100,"% of variance)"))+
              labs(title="PCA of samples"))
    }

  } else {
    pdf(outfile)
    for (i in 1:dim(combn(n,2))[2]){
      x_i <- rotation.x.pca[,combn(n,2)[1,i]]
      y_i <- rotation.x.pca[,combn(n,2)[2,i]]
      margin_x <- (max(x_i) - min(x_i)) * 0.2
      margin_y <- (max(y_i) - min(y_i)) * 0.2
      print(ggplot(rotation.x.pca, aes(x_i, y_i, color = ID)) +
              geom_point(shape = 16) +
              geom_label_repel(aes(label=colnames(x)))+
              xlim(min(x_i)-margin_x, max(x_i)+margin_x)+
              ylim(min(y_i)-margin_y,max(y_i)+margin_y)+
              xlab(paste0(colnames(rotation.x.pca)[combn(n,2)[1,i]]," (", var.x.pca[combn(n,2)[1,i]]*100,"% of variance)"))+
              ylab(paste0(colnames(rotation.x.pca)[combn(n,2)[2,i]]," (", var.x.pca[combn(n,2)[2,i]]*100,"% of variance)"))+
              labs(title="PCA of samples"))
    }
    dev.off()
  }
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
#' rr_LMCN.v2.split <- partition_to_uniques(rr_LMCN.v2[,-1], sample_attributes_LMCN, "replicate_name")
#' The first column of the rr_LMCN.v2 contained transcript IDs and was thus excluded from input.
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
#' rr_LMCN.v2.pairwise <- TER_all_pairs(rr_LMCN.v2.split, sample_attributes_LMCN, "read_type", "replicate_name", "cell_line")
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



#' @title pairs2pi0s
#' @description Function to estimate and plot the proportion of null features from pairwise TER tests.
#' @param x A list of the results of TER tests between all pairs of samples in a data set produced by \code{\link{TER_all_pairs}}.
#' @param outfile The path and name of the output pdf file. Default: NULL.
#' @return A data frame describing uniqueIDs and groupIDs of pairs of compared samples, the \code{pair_type ("Homo" or "Hetero")}
#' and the estimated \emph{pi0} (proportion of null features aka not differentially translated transcripts).
#' @details A histogram of \emph{pi0}s is created colored by pair type. If \code{outfile} is given, the histogram will be saved to the pdf file, too.
#' \code{"Homo"} pairs are expected to have a higher proportion of null features than \code{"Hetero"} pairs.
#' @examples
#' pi0df_LMCN <- pairs2pi0s(rr_LMCN.v2.pairwise)
#' @export
pairs2pi0s <- function(x, outfile = NULL){
  pi0df <- data.frame(t(sapply(x, function(y) c(y[[1]], y[[2]], y[[3]], qvalue::pi0est(y[[4]][, 8])$pi0))))
  names(pi0df) <- c("uniqueID1", "uniqueID2", "groupID1", "groupID2", "pair_type", "pi0")
  pi0df$pi0 <- as.numeric(as.character(pi0df$pi0))
  pi0hist <- ggplot(pi0df, aes(x=as.numeric(pi0), fill=as.factor(pair_type))) +
    geom_histogram(binwidth = 0.05) + guides(fill=guide_legend(title="Pair type")) +
    labs(x="Proportion of NULL features")
  print(pi0hist)
  if (!is.null(outfile)){
    pdf(outfile)
    print(pi0hist)
    dev.off()
  }
  return(pi0df)
}



#' @title generate_correlogram
#' @description Function to calculate and plot the correlation matrix of TER test z scores.
#' @param x A list of TER test outputs. Each element of the list is a data frame produced by the \code{\link{logit_seq}}
#' function comparing two samples.
#' @details
#' Columns 1-4 of each output data frame describe \emph{Estimate}, \emph{SD(Estimate)}, \emph{z score} and \emph{p-value}
#' of the logistic regression intercept. Columns 5-8 of the data frames describe \emph{Estimate}, \emph{SD(Estimate)},
#' \emph{z score} and \emph{p-value} of the independent variable (predictor) of the TER test. The \code{\link{generate_correlogram}}
#' function extracts the 7th column from all data frames and calculates and plots their correlation matrix. This function is an internal
#' component of the \code{\link{pairs2correlograms}} function.
#' @return Correlation matrix of z scores.
#' @export
generate_correlogram <- function(x){
  xz <- sapply(x, function(y) y$fit[,7])
  xz_cor <- cor(xz)
  print(corrplot::corrplot(xz_cor, method="color", addCoef.col = "white"))
  return(xz_cor)
}



#' @title pairs2correlograms
#' @description Function to calculate and plot correlograms from equivalent pairwise TER tests.
#' @param x A list of the results of TER tests between all pairs of samples in a data set produced by \code{\link{TER_all_pairs}}.
#' @param outfile The path and name of the output pdf file containing the correlograms (correlation matrix heat maps). Default: NULL.
#' @details
#' The Ribolog TER test can be performed on single replicates per biological sample. In a replicated experiment
#' such as (sample A: reps A1 and A2 + sample B: reps B1 and B2), correlation coefficients of regression z scores from equivalent tests
#' (A1 vs B1, A2 vs B1, A1 vs B2, A2 vs B2) are used to evaluate replicate homogeneity and help determine
#' the minimum advisable number of replicates to achieve reproducibility.
#' @return A list containing the correlation matrices of equivalent replicate-by-replicate TER tests in a data set.
#' @examples
#' rr_LMCN.v2.correlograms <- pairs2correlograms(rr_LMCN.v2.pairwise)
#' @export
pairs2correlograms <- function(x, outfile = NULL){
  xhets <- rlist::list.filter(x, pair_type == "hetero")
  xhets_grouped <- rlist::list.group(xhets, groupIDs)
  xzcors <- lapply(xhets_grouped, function(x) Ribolog::generate_correlogram(x))
  pdf(outfile)
  xzcors <- lapply(xhets_grouped, function(x) Ribolog::generate_correlogram(x))
  dev.off()
  return(xzcors)
}
