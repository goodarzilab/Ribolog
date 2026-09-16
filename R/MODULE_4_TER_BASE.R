#' @importFrom data.table data.table setnames setkey setcolorder setorder tstrsplit CJ .N .SD .EACHI :=
#' @import ggplot2
#' @import ggrepel
#' @importFrom dplyr %>% count rename filter_all all_vars
#' @import corrplot
#' @import rlist
#' @import EnhancedVolcano
#' @import fitdistrplus



#' @title logit_seq
#' @description Function to perform the logistic regression test for differential translational efficiency
#' @param x Input data frame where each column contains RNA or RPF counts from a single sample.
#' Rows are genes/transcripts.
#' @param design Design matrix of the experiment describing samples and their attributes.
#' i-th row in the design matrix corresponds to the i-th column in the input data frame.
#' @param model Regression equation modeling the odds ratio of RPF/RNA counts against the selected design variables (sample attributes)
#' @param adj_method P-value adjustment method.
#' Options: "qvalue", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' "qvalue" calls the \emph{qvalue} package. Other methods are from base R.
#' @param feature_list (Optional) A vector containing IDs of genes/transcripts.
#' Must have the same length as the row number of input data frame.
#' @param id_mapper (Optional) A data frame mapping \code{feature_list} IDs (e.g. transcript IDs) to another
#' ID system (e.g. gene IDs/names), such as \code{\link{human_id_mapper}}. Must contain a \code{transcript}
#' column matching \code{feature_list} and a column named \code{id_column}. When given, the output's row
#' names are replaced with the mapped IDs (falling back to the original \code{feature_list} value for any
#' unmapped feature, and disambiguated with \code{\link{make.unique}} if the mapping is not one-to-one).
#' If \code{NULL} (default), row names are just \code{feature_list}, unchanged.
#' @param id_column Name of the column in \code{id_mapper} to map to. Default: "gene_id".
#' @param long_output If \code{TRUE} (default), report five columns per model term: \code{logFC_<term>},
#' \code{std_<term>}, \code{zvalue_<term>}, \code{pvalue_<term>}, and (if \code{adj_method != "none"})
#' \code{corrected_pvalue_<term>}. If \code{FALSE}, report only \code{logFC_<term>} and
#' \code{corrected_pvalue_<term>} (or \code{pvalue_<term>} if \code{adj_method == "none"}).
#' @return A data frame containing the output of the regression, with columns grouped together per model
#' term (including the intercept, named "intercept"). If long output is requested, five columns are
#' reported for each term: \code{logFC_<term>} (regression coefficient), \code{std_<term>} (standard
#' deviation of the estimated coefficient), \code{zvalue_<term>}, \code{pvalue_<term>} (Wald-test p-value),
#' and \code{corrected_pvalue_<term>} (multiple-testing-adjusted p-value; only present when
#' \code{adj_method != "none"}). If long output is not specified, only \code{logFC_<term>} and
#' \code{corrected_pvalue_<term>} (or \code{pvalue_<term>} if \code{adj_method == "none"}) are reported.
#' Gene/transcript IDs provided by the 'feature_list' argument is added to the output matrix as row names.
#' Generalizes to any number of model terms/predictors.
#' @details The response variable is always the read type which should be a factor variable with two levels: "RNA" and "RPF".
#' Translational efficiency of each sample if formulated as the odds of RPF vs. RNA reads.
#' The change in translational efficiency between samples or based on unit values of any predictor
#' (design variable) is given by the exponentiated regression coefficient.
#' Exponentiated intercept gives the TE for the sample with all the attributes at the reference level.
#' If the predcitor variable is categorical, its levels are sorted alphabetically, the first level is set to reference
#' and all other levels are compared to it. The reference level can be manually changed using the \code{\link[stats]{relevel}} function of R.
#' @examples
#' # Test the effect of lung metastasis on translational efficiency:
#' fit1 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # Test the effects of lung metastasis and cell line origin on translational efficiency:
#' fit2 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis + cell_line_origin, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # Test the effects of lung metastasis, cell line origin and their interaction on translational efficiency:
#' fit3 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis * cell_line_origin, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # Test the effect of cell line on translational efficiency (cell line "CN34" is used as reference because it comes first alphabetically):
#' fit4 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ cell_line, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # Test the effect of cell line on translational efficiency with cell line "MDA" set as reference:
#' sample_attributes2 <- Ribolog::sample_attributes_dummy
#' sample_attributes2$cell_line <- relevel(sample_attributes2$cell_line, ref = "MDA")
#' fit5 <- Ribolog::logit_seq(rr.v2_dummy[,-1], sample_attributes2, read_type ~ cell_line, as.vector(rr.v2_dummy$transcript), adj_method = 'none')
#' # Same test, with row names mapped from transcript IDs to gene IDs using a bundled ID mapper:
#' fit1_FDR <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                                 as.vector(rr.v2_dummy$transcript), adj_method = 'fdr', id_mapper = Ribolog::human_id_mapper)
#' @export

logit_seq <- function(x, design, model, adj_method, feature_list=NULL, long_output = TRUE,
                      id_mapper = NULL, id_column = "gene_id"){
  logit_seq_gene <- function(m){
    prep <- data.frame(design,m)
    fit <- suppressWarnings(glm(model, data=prep, family="binomial"(link="logit"), weights = m))
    sfit <- summary(fit)
    return(c(t(sfit$coefficients)))
  }

  if (!is.null(feature_list)) {
    rownames(x) <- feature_list
  } else if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
  }

  total_rpf_counts <- rowSums(x[,design$read_type=='RPF'])
  total_rna_counts <- rowSums(x[,design$read_type=='RNA'])

  empty_rpf_transcripts <- rownames(x)[total_rpf_counts == 0]
  empty_rna_transcripts <- rownames(x)[total_rna_counts == 0]

  model_string <- paste(deparse(model), collapse = " ")

  if (length(empty_rpf_transcripts) > 0){

    cat(sprintf("While running the comparison for %s, %s transcript(s) were removed because they had 0 counts across all the RPF samples being compared.\n",
                model_string, length(empty_rpf_transcripts)))
    x <- Ribolog::min_count_filter(x, mincount = 5, columns = design$read_type=='RPF' , method = "all")
  }

  if (length(empty_rna_transcripts) > 0){

    cat(sprintf("While running the comparison for %s, %s transcript(s) were removed because they had 0 counts across all the RNA samples being compared.\n",
                model_string, length(empty_rna_transcripts)))
    x <- Ribolog::min_count_filter(x, mincount = 2, columns = design$read_type=='RNA', method = "average")
  }

  logit_x <- t(apply(x, 1, logit_seq_gene))
  prep1 <- cbind.data.frame(design, data.frame(counts1 = as.numeric(x[1,])))
  fit1 <- suppressWarnings(glm(model, data=prep1, family="binomial"(link="logit"), weights = counts1))
  sfit1 <- summary(fit1)

  term_names <- rownames(sfit1$coefficients)
  clean_terms <- ifelse(term_names == "(Intercept)", "intercept", term_names)
  n_terms <- length(term_names)

  colnames(logit_x) <- apply(expand.grid(c("logFC", "std", "zvalue", "pvalue"), clean_terms), 1, paste, collapse = "_")
  rownames(logit_x) <- rownames(x)

  if (!is.null(id_mapper)) {
    if (!all(c("transcript", id_column) %in% colnames(id_mapper))) {
      stop(paste0("id_mapper must contain columns 'transcript' and '", id_column, "'."))
    }
    id_map <- setNames(as.character(id_mapper[[id_column]]), as.character(id_mapper$transcript))
    surviving_features <- rownames(logit_x)
    mapped_ids <- unname(id_map[surviving_features])
    unmapped <- is.na(mapped_ids)
    mapped_ids[unmapped] <- surviving_features[unmapped]
    if (any(unmapped)) {
      warning(sprintf("%s of %s features were not found in id_mapper; kept their original feature_list ID.",
                      sum(unmapped), length(surviving_features)))
    }
    rownames(logit_x) <- make.unique(mapped_ids)
  }

  if (adj_method != 'none'){
    pcols <- seq(4, 4 * n_terms, by = 4)
    logit_x <- Ribolog::adj_TER_p(logit_x, pcols = pcols, adj_method = adj_method)

    corrected_cols <- (4 * n_terms + 1):(4 * n_terms + n_terms)
    colnames(logit_x)[corrected_cols] <- paste0("corrected_pvalue_", clean_terms)

    reordered <- unlist(lapply(seq_len(n_terms), function(i){
      block <- ((i - 1) * 4 + 1):(i * 4)
      c(block, 4 * n_terms + i)
    }))
    logit_x <- logit_x[, reordered, drop = FALSE]
  }

  if (long_output == FALSE){
    block_size <- if (adj_method != 'none') 5 else 4
    keep <- unlist(lapply(seq_len(n_terms), function(i) c((i - 1) * block_size + 1, i * block_size)))
    logit_x <- logit_x[, keep, drop = FALSE]
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
#' fit1 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                             as.vector(rr.v2_dummy$transcript), adj_method = 'none', long_output = FALSE)
#' fit1_qval <- adj_TER_p(fit1, c(2,4), "qvalue")
#' fit3 <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis * cell_line_origin,
#'                             as.vector(rr.v2_dummy$transcript), adj_method = 'none', long_output = FALSE)
#' fit3_fdr <- adj_TER_p(fit3, c(2,4,6,8), "fdr")
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

#' @title volcano_plot
#' @description Function to make the volcano plot from the output of \code{\link{logit_seq}}.
#' @param fit Output object of a \code{\link{logit_seq}} call.
#' @param covariate (Optional) Name (or unambiguous prefix) of the model term to plot, e.g. "intercept" or
#' "lung_metastasis" - matched against the \code{logFC_<term>} columns of \code{fit} (so "lung_metastasis"
#' matches a term stored as "lung_metastasisY"). When given, this sets \code{x} to \code{logFC_<term>} and
#' \code{y} to \code{corrected_pvalue_<term>} (or \code{pvalue_<term>} if \code{fit} has no corrected p-values),
#' overriding \code{x}/\code{y}. Leave \code{NULL} to specify \code{x}/\code{y} manually instead.
#' @param x Data column denoting the log fold change. Ignored if \code{covariate} is given.
#' @param y Data column denoting the p-value. Ignored if \code{covariate} is given.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Plot title.
#' @param subtitle Plot subtitle. Default: \code{NULL} (no subtitle; EnhancedVolcano's own default would
#' otherwise print the italic word "EnhancedVolcano").
#' @param selectLab (Optional) Character vector of specific gene/transcript labels (matching
#' \code{rownames(fit)}) to display. Default: \code{NULL}, which lets EnhancedVolcano label every point
#' that passes both \code{pCutoff} and \code{FCcutoff} (i.e. the most significant hits).
#' @param labFace Font face for point labels. Default: "bold" (plain text is hard to see against the
#' point cloud); pass "plain" to disable.
#' @param pCutoff P-value significance threshold.
#' @param FCcutoff Fold change significance threshold (e.g. 1.5 = 1.5x), not log fold change. \code{fit}'s
#' \code{logFC_<term>} column is on the natural-log scale, so this is converted internally
#' (\code{log(FCcutoff)}) before being passed to EnhancedVolcano.
#' @param xlim,ylim Axis limits.
#' @param titleLabSize,border Passed to \code{EnhancedVolcano::EnhancedVolcano}.
#' @return A volcano plot visualisation
#' @details A visualisation of the volcano plot resulting from the p-values and log fold change.
#' EnhancedVolcano emits a few warnings that are expected and cosmetic (handling of exact-zero p-values,
#' and its own internal use of a deprecated ggplot2 `size` argument); these are suppressed here.
#' @examples
#' fit1_FDR <- Ribolog::logit_seq(rr.v2_dummy[,-1], Ribolog::sample_attributes_dummy, read_type ~ lung_metastasis,
#'                                 as.vector(rr.v2_dummy$transcript), adj_method = 'fdr',
#'                                 id_mapper = Ribolog::human_id_mapper, id_column = 'gene_name')
#' vplot <- volcano_plot(fit1_FDR, covariate = "lung_metastasis")
#' vplot <- volcano_plot(fit1_FDR, covariate = "lung_metastasis", selectLab = c("CYP26B1"))
#' @export

volcano_plot <- function(fit, covariate = NULL, x = "logFC_lung_metastasisY",
y = "corrected_pvalue_lung_metastasisY", xlab = "Ln fold change", ylab = "-Log10 FDR",
title = "LMCN data, metastatic vs non-metastatic", subtitle = NULL, selectLab = NULL,
labFace = "bold",
titleLabSize = 12, border = "full",
pCutoff = 0.001, FCcutoff = 1.5, xlim = c(-5, 5), ylim = c(0, 10)) {

  if (!is.null(covariate)) {
    logfc_cols <- grep("^logFC_", names(fit), value = TRUE)
    terms <- sub("^logFC_", "", logfc_cols)
    match_idx <- which(startsWith(terms, covariate))

    if (length(match_idx) == 0) {
      stop(paste0("No logFC column found matching covariate '", covariate, "'. Available terms: ",
                  paste(terms, collapse = ", ")))
    }
    if (length(match_idx) > 1) {
      stop(paste0("covariate '", covariate, "' matches more than one term: ",
                  paste(terms[match_idx], collapse = ", "), ". Use a more specific value."))
    }

    term <- terms[match_idx]
    x <- paste0("logFC_", term)
    y <- if (paste0("corrected_pvalue_", term) %in% names(fit)) paste0("corrected_pvalue_", term) else paste0("pvalue_", term)
  }

  if (! x %in% names(fit)){
    stop(print(paste('The column', x, 'does not exist in the given dataframe.')))
  }

  if (! y %in% names(fit)){
    stop(print(paste('The column', y, 'does not exist in the given dataframe.')))
  }

  # logFC_<term> is on the natural-log scale (glm binomial/logit coefficients), but FCcutoff is a real
  # fold change (e.g. 1.5 = 1.5x); convert to the matching ln-scale threshold for EnhancedVolcano, which
  # otherwise compares FCcutoff directly against the (ln-scale) x column.
  lnFCcutoff <- log(FCcutoff)

  return(suppressWarnings(EnhancedVolcano::EnhancedVolcano(fit, lab = rownames(fit), x=x, xlab=xlab, y=y, ylab=ylab,
  title=title, subtitle=subtitle, selectLab=selectLab, labFace=labFace,
  titleLabSize=titleLabSize, border=border, pCutoff=pCutoff, FCcutoff=lnFCcutoff, xlim=xlim, ylim=ylim)))
}
