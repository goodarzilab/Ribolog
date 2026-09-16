#' @importFrom data.table data.table setnames setkey setcolorder setorder tstrsplit CJ .N .SD .EACHI :=
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>% count rename filter_all all_vars
#' @import corrplot
#' @import rlist
#' @import EnhancedVolcano
#' @import fitdistrplus


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
#' te.v2 <- Ribolog::create_te(rr.v2_dummy, idcolumns = 1, rnacolumns = c(2:9), rpfcolumns = c(10:17))
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
#' te.v2 <- Ribolog::create_te(rr.v2_dummy, idcolumns = 1, rnacolumns = c(2:9), rpfcolumns = c(10:17))
#' te.v2.cent <- Ribolog::row_center(te.v2, columns = c(2:9))
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
#' te.v2 <- Ribolog::create_te(rr.v2_dummy, idcolumns = 1, rnacolumns = c(2:9), rpfcolumns = c(10:17))
#' te.v2.stnd <- Ribolog::row_standardize(te.v2, columns = c(2:9))
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
#' te.v2 <- Ribolog::create_te(rr.v2_dummy, idcolumns = 1, rnacolumns = c(2:9), rpfcolumns = c(10:17))
#' te.v2.stnd <- Ribolog::row_standardize(te.v2, columns = c(2:9))
#' # The first column of input data (transcript ID) had to be removed to create a full-numeric input dataset.
#' pca_qc(te.v2.stnd[, -1], n = 2, ID = Ribolog::sample_attributes_dummy$cell_line[c(1:8)])
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


#' @title pairs2pi0s
#' @description Function to estimate and plot the proportion of null features from pairwise TER tests.
#' @param x A list of the results of TER tests between all pairs of samples in a data set produced by \code{\link{TER_all_pairs}}.
#' @param outfile The path and name of the output pdf file. Default: NULL.
#' @return A data frame describing uniqueIDs and groupIDs of pairs of compared samples, the \code{pair_type ("Homo" or "Hetero")}
#' and the estimated \emph{pi0} (proportion of null features aka not differentially translated transcripts).
#' @details A histogram of \emph{pi0}s is created colored by pair type. If \code{outfile} is given, the histogram will be saved to the pdf file, too.
#' \code{"Homo"} pairs are expected to have a higher proportion of null features than \code{"Hetero"} pairs.
#' @examples
#' # Subset to a handful of transcripts so the pairwise fits run quickly.
#' rr.v2_dummy.split <- Ribolog::partition_to_uniques(rr.v2_dummy[1:200, -1], Ribolog::sample_attributes_dummy, "replicate_name")
#' rr.v2_dummy.pairwise <- Ribolog::TER_all_pairs(rr.v2_dummy.split, Ribolog::sample_attributes_dummy,
#'                                           "read_type", "replicate_name", "cell_line", adj_method = "none")
#' pi0df <- pairs2pi0s(rr.v2_dummy.pairwise)
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
#' @details
#' Each pairwise fit in \code{x} can have silently dropped a different subset of transcripts (e.g. 0
#' counts in one read type within that specific pair; see \code{\link{logit_seq}}), so this function
#' aligns all fits to the transcripts common to every one of them before computing the correlation
#' matrix - otherwise the per-fit z score vectors would differ in length/order and either error out or
#' silently misalign.
#' @export
generate_correlogram <- function(x){
  common <- Reduce(intersect, lapply(x, function(y) rownames(y$fit)))
  if (length(common) == 0) {
    stop("No transcripts survived filtering in every fit being correlated; cannot compute a correlation matrix.")
  }
  xz <- sapply(x, function(y) y$fit[common, 7])
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
#' # Subset to a handful of transcripts so the pairwise fits run quickly.
#' rr.v2_dummy.split <- Ribolog::partition_to_uniques(rr.v2_dummy[1:200, -1], Ribolog::sample_attributes_dummy, "replicate_name")
#' rr.v2_dummy.pairwise <- Ribolog::TER_all_pairs(rr.v2_dummy.split, Ribolog::sample_attributes_dummy,
#'                                           "read_type", "replicate_name", "cell_line", adj_method = "none")
#' rr.v2_dummy.correlograms <- pairs2correlograms(rr.v2_dummy.pairwise)
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
