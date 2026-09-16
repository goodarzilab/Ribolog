#' @importFrom data.table data.table setnames setkey setcolorder setorder tstrsplit CJ .N .SD .EACHI :=
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>% count rename filter_all all_vars
#' @import corrplot
#' @import rlist
#' @import EnhancedVolcano
#' @import fitdistrplus
#' @import nnet
#' @import tidyr



#' @title get_tr_regions
#' @description Function to count the number of reads mapping to each of the 3 regions per transcript per sample.
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @return A dataframe containing the number of reads mapping to each region, for each transcript, with a
#' \code{sample_name} column identifying each sample - named to match the \code{sample_name} convention
#' used elsewhere (e.g. a sample attributes/design table), so it can be used directly as a \code{sample_ID}.
#' @examples
#' tr_regions_df <- get_tr_regions(reads_psite_list)
#' @export

get_tr_regions <- function(reads_psite_list){
  rpl <- lapply(reads_psite_list, function (x) dplyr::select(x, transcript, psite_region))
  rpl.2 <- lapply( rpl , function(x) x %>% count(transcript, psite_region)  )
  rp.df <- data.table::rbindlist(rpl.2, idcol="sample_name")
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
#' @param normalization_factors A named or ordered vector of per-sample normalization factors (e.g. size factors),
#' one per sample, in the same sample order as the columns produced when \code{edf} is pivoted wide.
#' @return data frame of the normalised transcript counts for each region
#' @export

normalize_with_ratios <- function(edf, num_samples, normalization_factors){
  edf <- edf %>% pivot_wider(names_from = sample_name, values_from = count)
  edf[is.na(edf)] = 0
  data_columns <- c(3:(2+num_samples))

  id_names <- names(edf)
  names(normalization_factors) <- names(edf[,data_columns])
  normalized_edf <- t(t(edf[,data_columns])/normalization_factors)
  edf[, data_columns] <- normalized_edf
  names(edf) <- id_names

  edf <- edf %>%  pivot_longer(all_of(data_columns), names_to = "sample_name", values_to = "count")
  return(edf)
}



#' @title replace_low_count_transcripts
#' @description Function to remove transcripts with low level counts across all three regions
#' @param edf Dataframe output of get_tr_regions()
#' @param mincount Minimum average of counts across all three region types
#' @param method Filtering method passed to \code{\link{min_count_filter}}. Default: "average".
#' @return data frame of the filtered transcript counts for each region
#' @details
#' A transcript with zero observed reads in a given region across every sample has no row for that
#' (transcript, psite_region) combination in \code{edf} at all (see \code{\link{get_tr_regions}}), so
#' pivoting to wide format introduces \code{NA} there rather than 0; these are filled with 0 before
#' filtering, since a missing combination genuinely means zero counts, not an unknown value.
#' @export

replace_low_count_transcripts <- function(edf, mincount = 2, method='average'){

    edf <- edf %>% pivot_wider(names_from = psite_region, values_from = count)
    # A transcript with zero observed reads in a region across every sample has no row for that
    # (transcript, psite_region) combination at all, so pivot_wider introduces NA here rather than 0.
    # That NA would otherwise flow into min_count_filter's row means/sums and, since indexing a data
    # frame with an NA condition returns an all-NA row (not an excluded one), silently corrupt those
    # transcripts into garbage NA rows instead of correctly filtering them out.
    edf[is.na(edf)] <- 0
    initial_num <- length(rownames(edf))

    edf <- Ribolog::min_count_filter(edf, mincount, c(3:5), method=method)
    final_num <- length(rownames(edf))

    print(paste('Number of records filtered out:', initial_num - final_num))
    edf <- edf %>%  pivot_longer(c(3:5), names_to = "psite_region", values_to = "count")
    return(edf)
}



#' @title normalize_and_clean
#' @description Normalize per-region transcript read counts using the median-of-ratios method
#' and remove transcripts with low counts across regions. Combines \code{\link{normalize_median_of_ratios}}
#' and \code{\link{replace_low_count_transcripts}} into a single step, so per-sample normalization
#' factors do not need to be supplied or computed by hand.
#' @param edf Dataframe output of get_tr_regions()
#' @param mincount Minimum average of counts across all three region types. Passed to \code{\link{replace_low_count_transcripts}}.
#' @param method Filtering method passed to \code{\link{replace_low_count_transcripts}}. Default: "average".
#' @return data frame of the normalised and filtered transcript counts for each region
#' @examples
#' tr_regions_df_pivot_norm_clean <- normalize_and_clean(tr_regions_df)
#' @export

normalize_and_clean <- function(edf, mincount = 2, method = 'average'){
  edf_wide <- edf %>% pivot_wider(names_from = sample_name, values_from = count)
  edf_wide[is.na(edf_wide)] <- 0
  data_columns <- c(3:ncol(edf_wide))

  edf_norm <- Ribolog::normalize_median_of_ratios(edf_wide, data_columns)
  edf_norm <- edf_norm %>% pivot_longer(all_of(data_columns), names_to = "sample_name", values_to = "count")

  edf_clean <- Ribolog::replace_low_count_transcripts(edf_norm, mincount = mincount, method = method)
  return(edf_clean)
}



#' @title remove_low_levels
#' @description Function to remove transcripts with not enough data to run a regression
#' @param xd Dataframe output of normalised and filtered counts for transcript and psite region, merged with the design matrix
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
#' @param adj_method P-value adjustment method. Default: "none".
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @return
#' Deviance test p-values (one per transcript).
#' @examples
#' fit.o <- tr_region_logit_dev(tr_regions_df_pivot_norm,
#'                               psite_region ~ lung_metastasis,
#'                               sample_attributes,
#'                               "sample_name")
#' @export

tr_region_logit_dev <- function(data, model, design = NULL, sample_ID = NULL, adj_method = 'none') {
  data <- Ribolog::remove_0_1_region_transcripts(data)
  xd <- merge(data, design, by = sample_ID)
  xd <- droplevels(xd)
  xd <- Ribolog::remove_low_levels(xd, model)
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



#' Rename tr_region_multi_logit's per-term output columns to a readable convention.
#' Not exported; used only internally by \code{\link{tr_region_multi_logit}}.
#' Maps e.g. "b..Intercept." -> "logFC_intercept", "p.lung_metastasisY" -> "pval_lung_metastasisY",
#' and "<adj_method>_p..Intercept." -> "pval_corrected_intercept", regardless of which term or
#' adjustment method was used.
#' @param nm Character vector of raw column names to translate.
#' @param adj_method P-value adjustment method used upstream (e.g. "fdr", "none"); determines the
#' \code{"<adj_method>_p."} prefix pattern to match.
#' @noRd
rename_multilogit_terms <- function(nm, adj_method){
  clean_term <- function(term) gsub("^\\.Intercept\\.$", "intercept", term)
  adj_pattern <- paste0("^", adj_method, "_p\\.")

  vapply(nm, function(n){
    if (grepl(adj_pattern, n)) {
      paste0("pval_corrected_", clean_term(sub(adj_pattern, "", n)))
    } else if (grepl("^b\\.", n)) {
      paste0("logFC_", clean_term(sub("^b\\.", "", n)))
    } else if (grepl("^se\\.", n)) {
      paste0("se_", clean_term(sub("^se\\.", "", n)))
    } else if (grepl("^z\\.", n)) {
      paste0("zscore_", clean_term(sub("^z\\.", "", n)))
    } else if (grepl("^p\\.", n)) {
      paste0("pval_", clean_term(sub("^p\\.", "", n)))
    } else {
      n
    }
  }, character(1), USE.NAMES = FALSE)
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
#' Multinomial test statistics for each transcript region. For each model term (e.g. "intercept",
#' "lung_metastasisY"), columns are named \code{logFC_<term>} (log fold change), \code{pval_<term>}
#' (unadjusted p-value) and \code{pval_corrected_<term>} (multiple-testing-adjusted p-value). When
#' \code{output = "long"}, \code{se_<term>} and \code{zscore_<term>} are also included.
#' @examples
#' multi_logit <- tr_region_multi_logit(data = tr_regions_df_pivot_norm_clean,
#'                                       model = psite_region ~ lung_metastasis,
#'                                       design = sample_attributes,
#'                                       sample_ID = "sample_name",
#'                                       output = "short", adj_method = "fdr")
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
    sfitx_com <- data.frame(sfitxb_df, sfitxp_df)
  } else if (output == "long") {
    sfitx_com <- data.frame(sfitxb_df, sfitxse_df, sfitxz_df, sfitxp_df)
  }

  # do.call(rbind.data.frame, <named list>) encodes each row's originating list name (the transcript ID)
  # into a "<transcript>.<outcome level>" row name (e.g. "ENST00000000233.3utr"); recover the transcript
  # ID from that suffix. A previous version of this function instead took the first coefficient column's
  # *value* (the (Intercept) beta, a number) and mislabeled it "transcript" - which silently discarded the
  # intercept term from the output, and crashed with "duplicate 'row.names' are not allowed" whenever two
  # transcripts happened to share the same intercept coefficient.
  strip_region_suffix <- function(df, suffix) {
    rownames(df) <- sub(paste0("\\.", suffix, "$"), "", rownames(df))
    df
  }

  output_df <- {}

  subset <- sfitx_com[grepl('3utr', rownames(sfitx_com), fixed=TRUE), , drop = FALSE]
  subset <- strip_region_suffix(subset, "3utr")

  subset_5 <- sfitx_com[grepl('5utr', rownames(sfitx_com), fixed=TRUE), , drop = FALSE]
  subset_5 <- strip_region_suffix(subset_5, "5utr")


  if (output == "short"){
    subset <- Ribolog::adj_TER_p(subset, pcols = c(3:4), adj_method = adj_method)
    subset_5 <- Ribolog::adj_TER_p(subset_5, pcols = c(3:4), adj_method = adj_method)
  } else if (output == "long") {
    subset <- Ribolog::adj_TER_p(subset, pcols = c(7:8), adj_method = adj_method)
    subset_5 <- Ribolog::adj_TER_p(subset_5, pcols = c(7:8), adj_method = adj_method)
  }

  names(subset) <- rename_multilogit_terms(names(subset), adj_method)
  names(subset_5) <- rename_multilogit_terms(names(subset_5), adj_method)

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
#' @param gene_mapper A dataframe containing transcripts and gene names. Since multiple transcripts
#' (isoforms) commonly share the same gene, mapped labels are disambiguated with \code{\link{make.unique}}
#' (e.g. "ABCF2", "ABCF2.1") to guarantee valid, unique row names; transcripts absent from \code{gene_mapper}
#' keep their original transcript ID as a fallback label.
#' @return
#' Volcano plot figure
#' @examples
#' visualize_orf_usage(multi_logit,
#'                      log_fold_change = 'logFC_lung_metastasisY',
#'                      p_val_column = 'pval_corrected_lung_metastasisY',
#'                      region = '3utr',
#'                      ylim = c(0, 20),
#'                      gene_mapper = mapper)
#' @export

visualize_orf_usage  <- function(subset, log_fold_change, p_val_column, region,
                                 pCutoff=0.01, FCcutoff=1.5, xlim=c(-5,5), ylim=c(0,50), gene_mapper=NULL) {

    subset <- subset[[region]]

    if (!is.null(gene_mapper)){

        if (!('transcript' %in% colnames(gene_mapper) && 'gene_name' %in% colnames(gene_mapper))) {
            stop('The gene mapper specified does not have a transcript or a gene_name column.')}

        map_vec <- setNames(as.character(gene_mapper$gene_name), as.character(gene_mapper$transcript))
        mapped_labels <- unname(map_vec[rownames(subset)])
        unmapped <- is.na(mapped_labels)
        mapped_labels[unmapped] <- rownames(subset)[unmapped]
        # Multiple transcripts (isoforms) commonly share the same gene, so gene names alone are not
        # guaranteed unique here; disambiguate (e.g. "ABCF2", "ABCF2.1") so they remain valid row names.
        rownames(subset) <- make.unique(mapped_labels)
    }

    # EnhancedVolcano emits a few warnings that are expected and cosmetic (handling of exact-zero
    # p-values, and its own internal use of a deprecated ggplot2 `size` argument); suppressed here,
    # same as Ribolog::volcano_plot.
    suppressWarnings(EnhancedVolcano::EnhancedVolcano(subset, lab = rownames(subset),
                                     x=log_fold_change, xlab= 'Ln Fold Change',
                                     y=p_val_column, ylab= '- Log10 P_val',
                                     title=paste0('Volcano plot for usage of ', region),
                                     pCutoff=pCutoff, FCcutoff=FCcutoff,
                                     xlim=xlim, ylim=ylim))
}
