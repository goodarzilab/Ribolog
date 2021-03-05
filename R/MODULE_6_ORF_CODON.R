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
#' @import rlists
#' @import gdata
#' @import nlme
#' @import EnhancedVolcano
#' @import fitdistrplus



#' @title get_tr_regions
#' @description Function to count the number of reads mapping to each of the 3 regions per transcript per sample.
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @return A dataframe containing the number of reads mapping to each region, for each transcript
#' @examples
#' tr_regions_df <- get_tr_regions(reads_psite_list)
#' @export

get_tr_regions <- function(reads_psite_list){
  rpl <- lapply(reads_psite_list, function (x) dplyr::select(x, transcript, psite_region))
  rpl.2 <- lapply( rpl , function(x) x %>% count(transcript, psite_region)  )
  rp.df <- data.table::rbindlist(rpl.2, idcol="sample")
  rp.df$psite_region <- relevel(as.factor(rp.df$psite_region), ref = "cds")
  names(rp.df)[4] <- "count"
  return(rp.df)
}



#' @title remove_0_1_region_transcripts
#' @description Function to remove transcripts with fewer than two active transcript sites from dataset
#' @param data Dataset containing site read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count").
#' Thus, each row in \code{data} contains the read count for one site - transcript - sample combination.
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @details
#' This function counts the number of pA sites with non-zero read counts for each transcripts and removes transripts
#' with fewer than two active pA sites. This is essential to avoid errors when running the regression models.
#' @return
#' A subset of the input dataset where all transcripts are guaranteed to have two or more non-zero read counts for transcript regions.
#' @export

remove_0_1_region_transcripts <- function(data){
  tt <- table(data$transcript, data$psite_region)
  ts <- rownames(tt[rowSums(tt>0) >= 2, , drop = FALSE])
  x2 <- subset(data, transcript %in% ts)
  x2 <- droplevels(x2)
  rate <- 100 * (1 - (length(ts)/dim(tt)[1]))
  print(paste0(rate, "% of transcripts had reads in <2 regions and were removed"))
  return(x2)
}



#' @title glm_deviance_test_p
#' @description Calculate the p-value from a deviance test comparing a model to its corresponding null.
#' @param x Output of a glm run
#' @return P-value calculated from the chisq test of deviance between the model and its corresponding null.
#' @export

glm_deviance_test_p <- function(x){
  p <- pchisq((x$null.deviance - x$deviance), df = (x$df.null - x$df.residual), lower.tail = FALSE)
  return(p)
}



#' @title normalize_with_ratios
#' @description Normalise the data with pre-provided normalisation ratios
#' @param edf Dataframe output of get_tr_regions()
#' @param num_samples Number of samples in the dataframe
#' @return data frame of the normalised transcript counts for each region
#' @export

normalize_with_ratios <- function(edf, num_samples, normalization_factors){
  edf <- edf %>% pivot_wider(names_from = sample, values_from = count)
  edf[is.na(edf)] = 0
  data_columns <- c(3:(2+num_samples))

  id_names <- names(edf)
  names(normalization_factors) <- names(edf[,data_columns])
  normalized_edf <- t(t(edf[,data_columns])/normalization_factors)
  edf[, data_columns] <- normalized_edf
  names(edf) <- id_names

  edf <- edf %>%  pivot_longer(all_of(data_columns), names_to = "sample", values_to = "count")
  return(edf)
}


#' @title tr_region_logit_dev
#' @description Function to evaluate the overall effect of predictors on transcript region count through a deviance test.
#' @param data Dataset containing transcript region read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count").
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @param model Regression model describing the dependence of region counts on sample attribute(s).
#' @param design (optional) Design matrix. A matrix describing sample attributes which can be used as predictors in the regression model.
#' @param sample_ID (optional) A key variable connecting the counts dataset (\code{data}) and the design matrix.
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return
#' Deviance test p-values (one per transcript).
#' @export

tr_region_logit_dev <- function(data, model, design = NULL, sample_ID = NULL, adj_method) {
  data <- Ribolog::remove_0_1_region_transcripts(data)
  xd <- merge(data, design, by = sample_ID)
  xd <- droplevels(xd)

  gfitx <- suppressWarnings(by(xd, xd$transcript, function(y) glm(as.formula(Reduce(paste,
    deparse(model)), env = new.env()), data = y, family = "binomial", weight = count)))

  dtest_gfitx <- lapply(gfitx, function(x) Ribolog::glm_deviance_test_p(x))
  dtest_gfitx.df <- tibble::rownames_to_column(as.data.frame(unlist(dtest_gfitx)), 'transcript') %>%
                        rename('p_devtest' = 'unlist(dtest_gfitx)' )

  if (adj_method != 'none'){
          dtest_gfitx.df <- Ribolog::adj_TER_p(dtest_gfitx.df, pcols = 2, adj_method = adj_method) }

  return(dtest_gfitx.df)
}
