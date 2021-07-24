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
#' @import nnet



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



#' @title replace_low_count_transcripts
#' @description Function to remove transcripts with low level counts across all three regions
#' @param edf Dataframe output of get_tr_regions()
#' @param mincount Minimum average of counts across all three region types
#' @return data frame of the filtered transcript counts for each region
#' @export

replace_low_count_transcripts <- function(edf, mincount = 2, method='average'){

    edf <- edf %>% pivot_wider(names_from = psite_region, values_from = count)
    initial_num <- length(rownames(edf))

    edf <- Ribolog::min_count_filter(edf, 2, c(3:5), method="average")
    final_num <- length(rownames(edf))

    print(paste('Number of records filtered out:', initial_num - final_num))
    edf <- edf %>%  pivot_longer(c(3:5), names_to = "psite_region", values_to = "count")
    return(edf)
}



#' @title remove_low_levels
#' @description Function to remove transcripts with not enough data to run a regression
#' @param edf Dataframe output of normalised and filtered counts for transcript and psite region
#' @param model Intended model to be used in the regression
#' @return data frame of the filtered transcript counts for each region
#' @export

remove_low_levels <- function(xd, model){
    initial_num <-  length(unique(xd$transcript))
    level_counts <- c()
    for (variable in all.vars(model) ) {
        level_counts[[variable]] <- by(xd, xd[,'transcript'], function(y) length(unique(y[[variable]])))}

    level_counts <- do.call(cbind, level_counts)
    level_counts <- as.data.frame(level_counts) %>% filter_all(all_vars(. > 1))

    xd <- xd[xd$transcript %in% rownames(level_counts), ]

    print(paste0( initial_num - length(unique(xd$transcript)),
        " transcripts were removed because there is less than 2 contrasting features in the variables of the model."))
    return(xd)
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



#' @title tr_region_multi_logit
#' @description Function to evaluate the overall effect of predictors on transcript region count through multinomial regression.
#' @param data Dataset containing transcript region read counts. This dataset must have a long shape, meaning that there should be only one
#' column containing read counts (and it MUST be named "count").
#' Other sample attributes beyond sample ID may be recorded in additional variables in this dataset, or provided separately through a design matrix
#' and a key variable (e.g. sample ID) connecting the \code{data} and \code{design} matrices.
#' @param model Regression model describing the dependence of region counts on sample attribute(s).
#' @param design (optional) Design matrix. A matrix describing sample attributes which can be used as predictors in the regression model.
#' @param sample_ID (optional) A key variable connecting the counts dataset (\code{data}) and the design matrix.
#' @param output A string either 'short' or 'long' to denote where the user wants the extra stats beyond p-values and logFC.
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return
#' Multinomial test statistics for each transcript region
#' @export

tr_region_multi_logit <- function(data, model, design, sample_ID, output = "short", adj_method = 'none'){
  data <- remove_0_1_region_transcripts(data)
  xd <- merge(data, design, by = sample_ID)
  xd <- droplevels(xd)
  xd <- xd[!is.na(xd$count),]
  xd <- remove_low_levels(xd, model)
  xd <- droplevels(xd)
  xd$psite_region <- factor(xd$psite_region, levels = c('cds', '3utr', '5utr'))

  fitx <- suppressWarnings(by(xd, xd$transcript, function(y) nnet::multinom(
      as.formula(Reduce(paste, deparse(model)), env = new.env()), data = y, weight = count, trace=FALSE)))

  sfitx <- suppressWarnings(lapply(fitx, function(x) summary(x)))
  sfitxb <- lapply(sfitx, function(x) x$'coefficients')
  sfitxse <- lapply(sfitx, function(x) x$'standard.errors')
  sfitxz <- lapply(sfitx, function(x) x$'coefficients' / x$'standard.errors')
  sfitxp <- lapply(sfitxz, function(x) (1 - pnorm(abs(x), 0, 1))*2)

  sfitxb_df <- do.call(rbind.data.frame, sfitxb)
  names(sfitxb_df) <- paste0("b.", names(sfitxb_df))

  sfitxse_df <- do.call(rbind.data.frame, sfitxse)
  names(sfitxse_df) <- paste0("se.", names(sfitxse_df))

  sfitxz_df <- do.call(rbind.data.frame, sfitxz)
  names(sfitxz_df) <- paste0("z.", names(sfitxz_df))

  sfitxp_df <- do.call(rbind.data.frame, sfitxp)
  names(sfitxp_df) <- paste0("p.", names(sfitxp_df))

  if (output == "short") {
    sfitx_com <- data.frame(cbind(transcript = sfitxb_df[,1], sfitxb_df[,-1], sfitxp_df))
  } else if (output == "long") {
    sfitx_com <- data.frame(cbind(transcript = sfitxb_df[,1], sfitxb_df[,-1], sfitxse_df, sfitxz_df, sfitxp_df))
  }

  output_df <- {}

  subset <-  sfitx_com[grepl('3utr', rownames(sfitx_com), fixed=TRUE),]
  rownames(subset) <- subset$transcript
  subset$transcript <- NULL

  subset_5 <-  sfitx_com[grepl('5utr', rownames(sfitx_com), fixed=TRUE),]
  rownames(subset_5) <- subset_5$transcript
  subset_5$transcript <- NULL


  if (output == "short"){
    subset <- Ribolog::adj_TER_p(subset, pcols = c(3:4), adj_method = adj_method)
    subset_5 <- Ribolog::adj_TER_p(subset_5, pcols = c(3:4), adj_method = adj_method)
  } else if (output == "long") {
    subset <- Ribolog::adj_TER_p(subset, pcols = c(7:8), adj_method = adj_method)
    subset_5 <- Ribolog::adj_TER_p(subset_5, pcols = c(7:8), adj_method = adj_method)
  }

  output_df[['3utr']] <- subset
  output_df[['5utr']] <- subset_5

  return(output_df)
}



#' @title visualize_orf_usage
#' @description Function to visualise the data for 3'UTR and 5'UTR usage
#' @param subset list of two datasets containing transcript region read counts for 5utr and 3utr
#' @param log_fold_change data column denoting the log fold change
#' @param p_val_column data column denoting the p-value
#' @param region A key variable denoting 3utr or 5utr
#' @param pCutoff A cutoff value for the p-values on the volcano plot
#' @param FCcutoff A cutoff value for the log fold change on the volcano plot
#' @param xlim A range for the log fold change on the volcano plot
#' @param ylim A range for the p-values on the volcano plot
#' @param gene_mapper A dataframe containing transcripts and gene names
#' @return
#' Volcano plot figure
#' @export

visualize_orf_usage  <- function(subset, log_fold_change, p_val_column, region,
                                 pCutoff=0.01, FCcutoff=1.5, xlim=c(-5,5), ylim=c(0,50), gene_mapper=NULL) {

    subset <- subset[[region]]

    if (!is.null(gene_mapper)){

        if (!('transcript' %in% colnames(gene_mapper) && 'gene_name' %in% colnames(gene_mapper))) {
            stop('The gene mapper specified does not have a transcript or a gene_name column.')}

        rownames(gene_mapper) <- gene_mapper$transcript
        rownames(subset) <- gene_mapper[rownames(subset),]$gene_name
    }

    EnhancedVolcano::EnhancedVolcano(subset, lab = rownames(subset),
                                     x=log_fold_change, xlab= 'Ln Fold Change',
                                     y=p_val_column, ylab= '- Log10 P_val',
                                     title=paste0('Volcano plot for usage of ', region),
                                     pCutoff=pCutoff, FCcutoff=FCcutoff,
                                     xlim=xlim, ylim=ylim)
}
