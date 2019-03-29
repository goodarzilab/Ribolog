#' @title gm_means_z
#' @description Function to calculate geometric mean of a vector of non-negative numbers
#' @details Zeros will remain in the calculation and render the geometric mean zero.
#' Compare to function: gm_means.
#' @param x Vector of non-negative numbers (may include zero)
#' @return Geometric mean
#' @export
gm_mean_z <- function(x){
  exp(sum(log(x)) / length(x))
}



#' @title gm_means
#' @description Function to calculate geometric mean of positive numbers in a vector
#' @details Zeros and negative numbers will be treated as 1:
#' They will not contribute to the sum of logs but will be included in calculation of vector length.
#' Compare to function: gm_mean_z.
#' @param x A vector of numbers (may include zero or negative numbers)
#' @return Geometric mean of the vector with the zeros and negative numbers replaced by 1
#' @export
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



#' @title normalize_median_of_ratios
#' @description Function to normalize a RNA-seq or ribo-seq dataset for library size variation
#' @details The original data columns are retained and normalized columns are added.
#' Compare to the method 'normalize_median_of_ratios2'.
#'
#' Use the 'data.columns' argument to exclude gene/transcript ID and other metadata columns from the
#' calculations or to normalize RNA and RPF counts separately while keeping them in the same data frame.
#' @param expression.data.frame A data frame containing RNA-seq, ribo-seq or similar data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param data.columns A vector of numbers specifying the columns to be normalized together.
#' @return An augmented data frame with all the original columns retained and normalized columns
#' appended at the end. New columns are named 'norm.<old_column>'.
#' @examples
#' rna.observed.normalized <- normalize_median_of_ratios(rna.observed, c(2:9))
#' rpf.observed.normalized <- normalize_median_of_ratios(rpf.observed, c(2:9))
#' ## The first column of each data frame (transcript ID) is excluded from normalization.
#' @export
normalize_median_of_ratios <- function (expression.data.frame, data.columns){
  gm_mean_z <- function(x){
    exp(sum(log(x)) / length(x))
  }
  edf <- expression.data.frame
  geo.mean.vec <- apply(edf[,data.columns], 1, function(x) gm_mean_z(x))
  ratios.df <- edf[,data.columns]/geo.mean.vec
  # Division by 0 gm_mean creates NAs here.
  normalization.factors <- apply(ratios.df, 2, function(x) median(x, na.rm=TRUE))
  # NAs are removed from calculation of median here.
  print("Normalization factors:")
  print(normalization.factors)
  normalized.edf <- t(t(edf[,data.columns])/normalization.factors)
  colnames(normalized.edf) <- paste0("norm.", colnames(normalized.edf))
  combined.df <- data.frame(edf,normalized.edf)
  return(combined.df)
}



#' @title normalize_median_of_ratios2
#' @description Function to normalize a RNA-seq or ribo-seq dataset for library size variation
#' @details The original columns of non-normalized counts are replaced in situ with normalized counts.
#' Compare to the method 'normalize_median_of_ratios'.
#'
#' Use the 'data.columns' argument to exclude gene/transcript ID and other metadata columns from the
#' calculations or to normalize RNA and RPF counts separately while keeping them in the same data frame.
#' @param expression.data.frame A data frame containing RNA-seq, ribo-seq or similar data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param data.columns A vector of numbers specifying the columns to be normalized together.
#' @return A data frame the same size as the input data frame where non-normalized counts in the 'data.columns'
#' are replaced in situ by normalized counts. Column names are not changed.
#' @examples
#' rna.observed.normalized2 <- normalize_median_of_ratios2(rna.observed, c(2:9))
#' rpf.observed.normalized2 <- normalize_median_of_ratios2(rpf.observed, c(2:9))
#' ## The first column of each data frame (transcript ID) is excluded from normalization.
#' @export
normalize_median_of_ratios2 <- function (expression.data.frame, data.columns){
  gm_mean_z <- function(x){
    exp(sum(log(x)) / length(x))
  }
  edf <- expression.data.frame
  id.names <- names(edf[-data.columns])
  geo.mean.vec <- apply(edf[,data.columns], 1, function(x) gm_mean_z(x))
  ratios.df <- edf[,data.columns]/geo.mean.vec
  # Division by 0 gm_mean will create NAs here.
  normalization.factors <- apply(ratios.df, 2, function(x) median(x, na.rm=TRUE))
  # NAs will be removed from calculation of median here.
  print("Normalization factors:")
  print(normalization.factors)
  normalized.edf <- t(t(edf[,data.columns])/normalization.factors)
  normalized.edf <- data.frame(edf[,-data.columns],normalized.edf)
  names(normalized.edf)[1:length(id.names)] <- id.names
  return(normalized.edf)
}
