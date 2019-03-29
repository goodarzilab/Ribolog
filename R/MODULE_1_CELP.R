#' @title print_read_ldist
#' @description Function to print out read length distributions of all samples from reads_list object to a pdf file
#' @param reads_list A read_list object produced by riboWaltz
#' @param outfile The path and name of the output pdf file
#' @param cl Confidence level, a number between 0 and 100, the percentage of all reads for which distribution of length is provided.
#' Use this argument to avoid printing out extremely uncommon read lengths. Default:99.
#' @details This function saves the read length distribution plots in a pdf file.
#' @import riboWaltz
#' @examples
#' print_read_ldist(reads_list, "<file.path>/RPF_Read_length_distributions2.pdf")
#' @export
print_read_ldist <- function(reads_list, outfile, cl=99){
  pdf(outfile)
  for (sample_i in names(reads_list)){
    length_dist_zoom <- riboWaltz::rlength_distr(reads_list, sample=sample_i, cl=cl)
    print(length_dist_zoom[["plot"]])
  }
  dev.off()
}



#' @title print_read_end_heatmap
#' @description Function to print out read end heatmaps of all samples from a riboWlatz reads_list object to a pdf file
#' @param reads_list A read_list object produced by riboWaltz
#' @param annotation An annotation data frame listing the transcripts and lengths of their 5'UTR, CDS and 3'UTR segments
#' @param outfile The path and name of the output pdf file
#' @param cl Confidence level, a number between 0 and 100, the percentage of all reads for which read end heat map is provided.
#' Use this argument to avoid printing out uncommon read lengths. Default:85.
#' @param utr5l Length of 5'UTR segment to be included upstream of translation start site. Default: 25.
#' @param cdsl Length of coding sequence to be included downstream of translation start site and upstream of translation termination site. Default: 40.
#' @param utr3l Length of 3'UTR segment to be included downstream of translation termination site. Default: 25.
#' @details Heatmaps of reads (stratified by read length) mapping to positions around the start and stop codon provide a visual
#' sense regarding a reasonable offset for p-site assignment as well as the 3-base periodicity of RPF reads.
#' @import riboWaltz
#' @examples
#' print_read_end_heatmap(reads_list, annotation.df, "<file.path>/RPF_read_end_heatmaps.pdf")
#' @export
print_read_end_heatmap <- function(reads_list, annotation, outfile, cl=85, utr5l = 25, cdsl = 40, utr3l = 25){
  pdf(outfile, width=20, height=10)
  for (sample_i in names(reads_list)){
    ends_heatmap_i <- riboWaltz::rends_heat(reads_list, annotation, sample=sample_i,
                                 cl=cl, utr5l = utr5l, cdsl = cdsl, utr3l = utr3l)
    print(ends_heatmap_i[["plot"]])
  }
  dev.off()
}



#' @title print_rop
#' @description Function to print out ribosome occupancy profiles of all samples from a riboWlatz reads_psite_list object to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by riboWaltz
#' @param annotation An annotation data frame listing the transcripts and lengths of their 5'UTR, CDS and 3'UTR segments
#' @param outfile The path and name of the output pdf file
#' @details Ribosome occupancy profiles provide visual sense of the positional distribution of reads mapping around start and stop codons.
#' @import riboWaltz
#' @examples
#' print_rop(reads_psite_list, annotation.df, "<file.path>/RPF_ribosome_occupancy_profiles2.pdf")
#' @export
print_rop <- function(reads_psite_list, annotation, outfile){
  pdf(outfile, width=20, height=10)
  for (sample_i in names(reads_psite_list)){
    metaprofile_psite_sample_i <- riboWaltz::metaprofile_psite(reads_psite_list, annotation, sample=sample_i)
    print(metaprofile_psite_sample_i[["plot"]])
  }
  dev.off()
}



#' @title print_period_region
#' @description Function to print out periodicity by region plots of all samples from a riboWlatz reads_psite_list object to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by riboWaltz
#' @param outfile The path and name of the output pdf file
#' @details These plots illustrate 3-base (reading-frame-dependent) periodicity of reads, separately for 5'UTR, CDS and 3'UTR.
#' @import riboWaltz
#' @examples
#' print_period_region(reads_psite_list, "<file.path>/Periodicity_by_region2.pdf")
#' @export
print_period_region <- function(reads_psite_list, outfile){
  pdf(outfile)
  for (sample_i in names(reads_psite_list)){
    frames_i <- riboWaltz::frame_psite(reads_psite_list, sample=sample_i, region="all")
    print(frames_i[["plot"]])
  }
  dev.off()
}



#' @title print_period_region_length
#' @description Function to print out periodicity by length and region plots of all samples from a riboWlatz reads_psite_list object to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by riboWaltz
#' @param outfile The path and name of the output pdf file
#' @param cl Confidence level, a number between 0 and 100, the percentage of all reads for which periodicity plots are provided.
#' Use this argument to avoid printing out uncommon read lengths. Default:90.
#' @import riboWaltz
#' @examples
#' print_period_region_length(reads_psite_list, "<file.path>/Periodicity_by_length_region2.pdf")
#' @details These plots illustrate 3-base (reading-frame-dependent) periodicity of reads, separately for 5'UTR, CDS and 3'UTR and stratified by length.
#' @export
print_period_region_length <- function(reads_psite_list, outfile, cl=90){
  pdf(outfile, width=20, height=15)
  for (sample_i in names(reads_psite_list)){
    frames_stratified_i <- riboWaltz::frame_psite_length(reads_psite_list, sample=sample_i, region="all", cl=cl)
    print(frames_stratified_i[["plot"]])
  }
  dev.off()
}



#' @title psite_to_codon_count
#' @description Function to count reads per codon for each transcript from a reads_psite_list object.
#' @param reads_psite_list A reads_psite_list object produced by riboWaltz
#' @param length_range Range of read lengths to be included in calculation of per codon RPF counts.
#' Longer and shorter reads may be excluded because of high uncertainty in offset assignment, lack of expected periodicity, etc.
#' @param annotation An annotation data frame listing the transcripts and lengths of their 5'UTR, CDS and 3'UTR segments
#' @param fasta.file The reference fasta file to which reads were mapped.
#' @import Biostrings
#' @import seqinr
#' @import plyr
#' @examples
#' tr_codon_read_count <- psite_to_codon_count(reads_psite_list, c(24:32), annotation.df, <file.path>/all_cdna_header_transcript_id_only_cds_se_longest.fa")
#' @return A list with the following structure: list$<sample.name>$<transcript.ID> data.frame: [1] codon_number [2] codon_type [3] aa_type [4] observed_count
psite_to_codon_count <- function(reads_psite_list, length_range, annotation, fasta.file){
  # Filter reads_psite_list to keep reads within the specified length range and mapping to CDS
  reads_psite_list.m <- lapply(reads_psite_list, function(x) subset(x, psite_region == "cds" & length %in% length_range))

  # Remove transcripts with CDS length non-divisible by 3.
  annotation.df.m <- annotation
  annotation.df.m$l_cds.by.3 <- annotation.df.m$l_cds %% 3
  annotation.df.m <- subset(annotation.df.m, l_cds.by.3 == 0)
  annotation.df.m <<- annotation.df.m[,-6]
  transcript_l3k <- as.character(annotation.df.m$transcript)

  # Function to create a list of data farmes, one for each transcript:
  # list$<transcript> [1] codon_number [2] codon_type [3] aa_type
  # Input: annotation data frame with transcript IDs and length of CDS and UTRs, reference fasta file and the genetic code
  # GENETIC_CODE comes from the Biostrings package.
  list_transcript_codons <- function(annotation, fasta.file, genetic.code){
    tr_codons_list <- list()
    #parse the fasta file into a list, each element named after a transcript, content is a vector of seq letters.
    fasta.parsed.t <- seqinr::read.fasta(file=fasta.file, seqtype="DNA", forceDNAtolower = FALSE)
    # Filter the list for those transcripts in the annotation
    accepted.transcripts <- as.character(annotation$transcript)
    fasta.parsed.t <- fasta.parsed.t[accepted.transcripts]
    # Create codon and aa list for each transcript
    no.tr <- length(accepted.transcripts)
    for (i in 1:no.tr){
      cds.start.i <- annotation.df.m[i,3] + 1
      cds.end.i <- annotation.df.m[i,3] + annotation.df.m[i,4]
      cds.bases.i <- fasta.parsed.t[i][[1]][cds.start.i:cds.end.i]
      codon_type <- paste0(cds.bases.i[c(TRUE, FALSE, FALSE)], cds.bases.i[c(FALSE, TRUE, FALSE)], cds.bases.i[c(FALSE, FALSE, TRUE)])
      aa_type <- genetic.code[codon_type]
      aa_type <- replace(aa_type, aa_type=='*', 'X')
      names(aa_type) <- NULL
      codon_number <- c(1:length(codon_type))
      df.i <- data.frame(codon_number, codon_type, aa_type)
      tr_codons_list[[i]] <- df.i
    }
    names(tr_codons_list) <- as.character(annotation$transcript)
    return(tr_codons_list)
  }

  tr_codons_aas <- list_transcript_codons(annotation = annotation.df.m,
                                          fasta.file = fasta.file,
                                          genetic.code = GENETIC_CODE)
  # Add a codon number column.
  calculate_codon_number <- function(sample_df){
    sample_df$codon_number <- (sample_df$psite_from_start %/% 3) + 1
    return(sample_df)
  }
  reads_psite_list.m <- lapply(reads_psite_list.m, calculate_codon_number)

  # Count the number of reads per codon for each transcript.

  # Function to modify an intermediate data frame
  trim_tr_read_df <- function(tr_read_df){
    df.temp <- tr_read_df[,-1]
    names(df.temp) <- c("codon_number", "read_count")
    row.names(df.temp) <- NULL
    return(df.temp)
  }
  tr_codon_read_count <- lapply(reads_psite_list.m, function(x) plyr::count(x, c("transcript", "codon_number")))
  tr_codon_read_count.s <- lapply(tr_codon_read_count, function(x) split(x, x$transcript))
  tr_codon_read_count.s2 <- lapply(tr_codon_read_count.s, function(x) lapply(x, trim_tr_read_df))

  # Initiate an all zero codon count base; then replace the non-zero counts from data.
  tr_codons_aas.z <- lapply(tr_codons_aas, function(x) data.frame(x, observed_count=rep(0, dim(x)[1])))
  sample_names <- names(reads_psite_list.m)
  tr_codon_read_count_zero_base <- list()
  for (s in sample_names){
    tr_codon_read_count_zero_base[[s]] <- tr_codons_aas.z
  }

  tr_codon_read_count_z_amended <- tr_codon_read_count_zero_base

  for (s in sample_names){
    for (t in transcript_l3k){
      if (t %in% names(tr_codon_read_count.s2[[s]])){
        tr_codon_read_count_z_amended[[s]][[t]]$observed_count <- tr_codon_read_count.s2[[s]][[t]]$read_count[match(tr_codon_read_count_z_amended[[s]][[t]]$codon_number, tr_codon_read_count.s2[[s]][[t]]$codon_number)]
        tr_codon_read_count_z_amended[[s]][[t]][is.na(tr_codon_read_count_z_amended[[s]][[t]])] <- 0
      }
    }
  }
  return(tr_codon_read_count_z_amended)
}



#' @title generate_rpf_loess
#' @description Function to compute per codon loess predicted read counts, bias coefficients and bias-corrected read counts
#' @param tr_codon_read_count_list A list generated by the psite_to_codon_count function
#' containing observed read counts per codon per trancript per sample
#' @param codon_radius Number of codons on either side of each codon influencing the loess prediction for the middle codon
#' @details This function is the heart of CELP method for stalling bias correction.
#' It starts with running a loess curve on per-codon read counts along the transcript to borrow information
#' from neighboring codons mitigating the uncertainty of p-site offset assignment and experimental stochasticity.
#' Then, a bias coefficient is calculated for each codon by integrating infromation on the excess of loess-corrected
#' read counts at that codons comapred to the transcript's background count across multiple samples. Finally, observed
#' read counts are divided by bias coefficients to calculate bias-corrected counts.
#' @return A list composed of two lists: 1. bias coefficients 2. bias-corrected read counts
#' The bias coefficient list has the following structure: list$<transcript.ID> data.frame: [1] codon_number [2] codon_type [3] aa_type [4] bias_coefficient.
#' The bias-corrected read count list has the following structure: list$<sample.name>$<transcript.ID> data.frame:
#' [1] codon_number [2] codon_type [3] aa_type [4] observed_count [5] loess_pred [6] loess_pred_by_nz_median [7] bias_coefficient [8] corrected_loess.
#' @examples
#' tr_codon_bias_coeff_corrected_count <- generate_rpf_loess (tr_codon_read_count, 5)
#' @export
generate_rpf_loess <- function(tr_codon_read_count_list, codon_raduis){
  tr_codon_read_count_loess_corrected <- tr_codon_read_count_list
  bias_coefficients_list <- list()
  sample_names_i <- names(tr_codon_read_count_loess_corrected)
  tr_names_i <- names(tr_codon_read_count_loess_corrected[[1]])
  for (t in tr_names_i){
    l_cds <- length(tr_codon_read_count_loess_corrected[[1]][[t]]$codon_number)
    span_tr <- (2*codon_raduis+1)/l_cds
    for (s in sample_names_i){
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred <-
        predict(loess(tr_codon_read_count_loess_corrected[[s]][[t]]$observed_count~tr_codon_read_count_loess_corrected[[s]][[t]]$codon_number, span=span_tr))
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred [tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred < 0 ] <- 0
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred_by_nz_median <-
        tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred / median(tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred[tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred>0])
    }

    bias_coefficients_list[[t]] <- data.frame(codon_number=tr_codon_read_count_loess_corrected[[1]][[t]]$codon_number)
    for (s in sample_names_i){
      bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]], tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred_by_nz_median)
    }
    names(bias_coefficients_list[[t]]) <- c("codon_number", sample_names_i)
    bias_coefficient <- apply(bias_coefficients_list[[t]][,-1], 1, function(x) gm_mean(x))
    bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]], bias_coefficient)
    for (s in sample_names_i){
      tr_codon_read_count_loess_corrected[[s]][[t]]$bias_coefficient <- bias_coefficients_list[[t]]$bias_coefficient
      tr_codon_read_count_loess_corrected[[s]][[t]]$corrected_loess <-
        tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred / tr_codon_read_count_loess_corrected[[s]][[t]]$bias_coefficient
    }
  }
  output <- list (bias_coefficients_list=bias_coefficients_list, tr_codon_read_count_loess_corrected=tr_codon_read_count_loess_corrected)

  return(output)
}



#' @title codon2transcript
#' @description Function to sum up codon counts per transcript
#' @param tr_codon_read_count_loess_corrected_list A list of per codon read counts for all samples and transcripts
#' produced by the generate_rpf_loess function. It has the following structure:
#' list$<sample.name>$<transcript.ID> data.frame:
#' [1] codon_number [2] codon_type [3] aa_type [4] observed_count [5] loess_pred [6] loess_pred_by_nz_median [7] bias_coefficient [8] corrected_loess.
#' @param count.type Options: "observed_count", "loess_pred", "corrected_loess"
#' @details Stalling bias correction is performed at the codon level but differential translational efficiency analysis
#' is ususally performed at transcript level. This function sums up the per codon counts (observed or corrected) on each transcript and
#' returns per transcript counts.
#' @examples
#' rpf.observed.sum <- codon2transcript(tr_codon_bias_coeff_loess_corrected_count$tr_codon_read_count_loess_corrected, "observed_count")
#' rpf.corrected.sum <- codon2transcript(tr_codon_bias_coeff_loess_corrected_count$tr_codon_read_count_loess_corrected, "corrected_loess")
#' @return A data frame where the first column is transcript IDs and the remaining columns contain per transcript read counts for all samples.
#' @export
codon2transcript <- function(tr_codon_read_count_loess_corrected_list, count.type){
  count.sum <- as.data.frame(sapply(tr_codon_read_count_loess_corrected_list, function(x) sapply(x, function(y) sum(y[count.type]))))

  count.sum$transcript <- rownames(count.sum)
  w <- dim(count.sum)[2]
  count.sum <- count.sum[,c(w,1:(w-1))]
  count.sum <- count.sum[order(count.sum$transcript),]
  rownames(count.sum) <- NULL
  return(count.sum)
}



