#' @title min_count_filter
#' @description Function to filter out genes with counts below a minimum in one or more samples.
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data.frame may contain additional columns for gene/transcript ID or other metadata.
#' @param min.count A single number (float), minimum RNA or RPF count required for a gene/transcript to be retained.
#' @param columns A vector specifying the columns to be considered for minimum count filtering.
#' @param method The method of filtering. Options: "all", "average". Default: "all".
#' @details If method="all" is chosen, a gene passes the filtering only if all samples specified
#' by the "columns" argument have values >= "min.count". If method="average" is chosen, a gene
#' passes the filtering if the average count among the specified "columns" is >= "min.count".
#'
#' Use the columns argument to exclude gene/transcript ID and other metadata columns from the
#' calculations or to filter RNA and RPF counts separately while keeping them in the same data.frame.
#' @return Filtered data.frame containing all the original columns but only the rows (genes)
#' that pass the filtering criterion.
#' @examples
#' rna.rpf.combined.m5 <- min_count_filter(rna.rpf.combined, 5, c(2:9))
#' @export
min_count_filter <- function(x, min.count, columns, method="all"){
  if (method=="all"){
    x <- x[!rowSums(x[,columns] < min.count),]
  } else if (method=="average"){
    x <- x[rowMeans(x[,columns]) >= min.count,]
  }
  return(x)
}



#' @title create_te
#' @description Function to create a TE (translational efficiency) data frame from a combined RNA+RPF data frame
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data.frame may contain additional columns for gene/transcript ID or other metadata.
#' @param id.columns A vector specifying the columns to be excluded from the calculations.
#' These columns may contain the gene/transcript ID or any metadata that need to be preserved.
#' @param rna.columns A vector specifying the columns containing RNA counts
#' @param rpf.columns A vector specifying the columns containing RPF counts
#' @details Translational efficiency is calculated as RPF/RNA.
#' There must be an equal number of RNA and RPF columns.
#' @return A data frame containing the original ID columns and the calculated TE columns.
#' @examples
#' te.m5 <- create_te(rna.rpf.combined.m5, 1, c(2:9), c(10:17))
#' @export
create_te <- function(x, id.columns=NULL, rna.columns, rpf.columns){
  y <- data.frame(x[,id.columns], x[,rpf.columns]/x[,rna.columns])
  names(y)[id.columns] <- names(x)[id.columns]
  names(y) <- gsub("rpf", "te", names(y))
  return(y)
}



#' @title row_center
#' @description Function to center a selected block of a dataframe on its row means
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param columns A vector specifying the columns to be included for row-centering.
#' @return A data frame where the specified columns from the input are row-centered and
#' the rest is intact.
#' @examples
#' rna.rpf.combined.m5.centered <- row_center(rna.rpf.combined.m5, c(2:9))
#' rna.rpf.combined.m5.centered2 <- row_center(rna.rpf.combined.m5.centered, c(10:17))
#' ## RNA data (columns 2:9) and RPF data (columns 10:17) are centered separately.
#' @export
row_center <- function(x, columns){
  x[,columns] <- t(apply(x[,columns], 1, function(y) y-mean(y)))
  return(x)
}



#' @title row_standardize
#' @description Function to standardize a selected block of a dataframe row-wise
#' @param x Input data frame containing RNA-seq or Ribo-seq data.
#' Rows are genes/transcripts and columns are samples.
#' The data frame may contain additional columns for gene/transcript ID or other metadata.
#' @param columns A vector specifying the columns to be included for row-wise standardization
#' @return A data frame where the specified columns from the input are row-standardized and
#' the rest is intact.
#' @details Row mean is subtracted from each element in the row and the result is divided by row standard deviation.
#' @examples
#' te.m5.standardized <- row_standardize(te.m5, c(2:9))
#' @export
row_standardize <- function(x, columns){
  x[,columns] <- t(apply(x[,columns], 1, function(y) (y-mean(y))/sd(y)))
  return(x)
}



#' @title pca_qc
#' @description Function to produce PCA output for data visualization and QC purposes
#' @param x A numeric matrix containing RNA, RPF, TE or some other type of data.
#' Rows are genes/transcripts and columns are samples.
#' Gene/transcript IDs and sample names may be given as row names and column names,
#' respectively (but not as additional non-numeric columns or rows).
#' @param n The number of principal components to be plotted
#' @param plotsfile The path and name of the output pdf file containing the PCA plots.
#' @details A summary of the PC analysis is printed out to standard output in addition to creating the plots file.
#' @import ggplot2
#' @examples
#' pca_qc(te.m5.standardized[,-1], 4, "<file.path>/PCA.QC.te.standardized.pdf")
#' ## The first column of input data (transcript ID) had to be removed to create a full-numeric input dataset.
#' @export
pca_qc <- function(x, n, plotsfile){
  x.pca <- prcomp(x[!rowSums(is.na(x)),])
  print(summary(x.pca))
  rotation.x.pca <- data.frame(x.pca$rotation)
  var.x.pca <- summary(x.pca)$importance[2,]
  pdf(plotsfile)
  for (i in 1:dim(combn(n,2))[2]){
    x_i <- rotation.x.pca[,combn(n,2)[1,i]]
    y_i <- rotation.x.pca[,combn(n,2)[2,i]]
    margin_x <- (max(x_i) - min(x_i)) * 0.2
    margin_y <- (max(y_i) - min(y_i)) * 0.2
    print(ggplot(rotation.x.pca, aes(x_i,y_i)) + geom_point() +
            geom_text(aes(label=colnames(x)), vjust=2)+
            xlim(min(x_i)-margin_x, max(x_i)+margin_x)+
            ylim(min(y_i)-margin_y,max(y_i)+margin_y)+
            xlab(paste0(colnames(rotation.x.pca)[combn(n,2)[1,i]]," (", var.x.pca[combn(n,2)[1,i]]*100,"% of variance)"))+
            ylab(paste0(colnames(rotation.x.pca)[combn(n,2)[2,i]]," (", var.x.pca[combn(n,2)[2,i]]*100,"% of variance)"))+
            labs(title="PCA of samples"))
  }
  dev.off()
}

