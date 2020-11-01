#' @import data.table
#' @import ggfortify
#' @import Biostrings
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import plyr
#' @import cowplot
#' @import robustbase
#' @import qvalue
#' @import nortest
#' @import fitdistrplus
#' @import matrixStats
#' @import sm
#' @import epiR
#' @import corrplot
#' @import mvmeta
#' @import DescTools
#' @import GenomicAlignments
#' @import corrplot
#' @import rlist
#' @import gdata
#' @import nlme



#' @title read_annotation
#' @description Function to create an annotation data table from a txt file.
#' @param annotation_file An annotation txt file listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' The annotation table has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' @details
#' Several \pkg{Ribolog} functions use the transcript annotation table. \code{\link{read_annotation}} reads the annotation
#' from a source file and stores it in a data table object which can be subsequently used by other functions.
#' Transcript names and segment lengths in the annotation must correspond to the reference sequences to which the reads were mapped.
#' We recommend creating the annotation txt file directly from the fasta file of cDNA sequences using the python script
#' \code{Biomart_cDNA_fasta_to_rW_annotation_and_reheadered_longest_CDS_cDNA_fasta.py} provided with the \pkg{Ribolog} package.
#' @note
#' Fasta files of cDNA sequences (only one transcript per gene with the longest CDS) and annotation tables for human and a number of popular model organisms
#' (\emph{Saccharomyces cerevisiae} yeast, \emph{Zea mays} maize, \emph{Arabidopsis thaliana} thale cress,
#' \emph{Drosophila melanogaster} fruit fly, \emph{Caenorhabditis elegans} worm, \emph{Danio rerio} zebra fish,
#' \emph{Mus musculus} mouse and \emph{Rattus norvegicus} rat) are included with the \pkg{Ribolog} package.
#' @export
#' @examples
#' annotation_human_cDNA <- read_annotation("<file.path>/Human.GRC38.96_annotation.txt")
read_annotation <- function(annotation_file){
  annotation <- as.data.table(read.table(annotation_file, header = TRUE))
  return(annotation)
}



#' @title bamtolist_rW
#' @description Function to convert bam files and an annotation file to a reads_list object.
#' @param bamfolder Path to the folder containing ribo-seq bam files (one bam file per sample is expected).
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param transcript_align A logical argument indicating whether the reads were mapped to transcripts (TRUE) or geneome + gtf (FALSE). Default: TRUE.
#' @param indel_threshold Maximum number of indels allowed per read. Default: 5.
#' @param ... ...
#' @details
#' This function takes a folder of ribo-seq bam files (one per sample) and an annotation table as input, and
#' produces a reads_list object.
#' @return
#' The output is a reads_list object each element of which is a data frame representing one sample. Each row in the data frame contains mapping
#' information (transcript name and read start and end coordinates) for one read.
#' @examples
#' reads_list_LMCN <- bamtolist_rW("<folder.path>/RPF_sorted_indexed", annotation_human_cDNA)
#' @export
bamtolist_rW <- function(bamfolder, annotation, transcript_align = TRUE, name_samples = NULL,
                         indel_threshold = 5, refseq_sep = NULL, granges = FALSE)
{
  names <- list.files(path = bamfolder, pattern = ".bam$")
  if (length(name_samples) == 0) {
    name_samples <- unlist(strsplit(names, ".bam"))
    names(name_samples) <- unlist(strsplit(names, ".bam"))
  }
  else {
    if (length(name_samples) > length(names)) {
      cat("\n")
      stop("length of name_samples greater than number of files\n\n")
    }
    if (length(name_samples) < length(names)) {
      cat("\n")
      stop("length of name_samples smaller than number of files\n\n")
    }
  }
  sample_reads_list <- list()
  for (n in names) {
    cat(sprintf("reading %s\n", n))
    sampname <- unname(name_samples[unlist(strsplit(n, ".bam"))])
    filename <- paste(bamfolder, n, sep = "/")
    dt <- as.data.table(GenomicAlignments::readGAlignments(filename))
    nreads <- nrow(dt)
    dt <- dt[, `:=`(diff_width, qwidth - width)][abs(diff_width) <=
                                                   indel_threshold]
    if (nreads != nrow(dt)) {
      cat(sprintf("%s M  (%s %%) reads removed: exceeding indel_threshold.\n",
                  format(round((nreads - nrow(dt))/1e+06, 2), nsmall = 2),
                  format(round(((nreads - nrow(dt))/nreads) * 100,
                               2), nsmall = 2)))
    }
    else {
      cat("Good! Number of indel below indel_threshold for all reads. No reads removed.\n")
    }
    dt <- dt[, .(seqnames, start, end, width, strand)]
    setnames(dt, c("transcript", "end5", "end3", "length",
                   "strand"))
    if (length(refseq_sep) != 0) {
      dt <- dt[, `:=`(transcript, tstrsplit(transcript,
                                            refseq_sep, fixed = TRUE, keep = 1))]
    }
    nreads <- nrow(dt)
    cat(sprintf("reads: %s M\n", format(round((nreads/1e+06),
                                              2), nsmall = 2)))
    dt <- dt[as.character(transcript) %in% as.character(annotation$transcript)]
    if (nreads != nrow(dt)) {
      if (nrow(dt) == 0) {
        stop(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n\n",
                     format(round((nreads - nrow(dt))/1e+06, 2),
                            nsmall = 2), format(round(((nreads - nrow(dt))/nreads) *
                                                        100, 2), nsmall = 2)))
      }
      else {
        cat(sprintf("%s M  (%s %%) reads removed: reference transcript IDs not found in annotation table.\n",
                    format(round((nreads - nrow(dt))/1e+06, 2),
                           nsmall = 2), format(round(((nreads - nrow(dt))/nreads) *
                                                       100, 2), nsmall = 2)))
      }
    }
    else {
      cat("Great! All reads' reference transcript IDs were found in annotation table. No reads removed.\n")
    }
    if (transcript_align == TRUE | transcript_align == T) {
      nreads <- nrow(dt)
      dt <- dt[strand == "+"]
      if (nreads != nrow(dt)) {
        cat(sprintf("%s M  (%s %%) reads removed: mapping on negative strand.\n",
                    format(round((nreads - nrow(dt))/1e+06, 2),
                           nsmall = 2), format(round(((nreads - nrow(dt))/nreads) *
                                                       100, 2), nsmall = 2)))
      }
      else {
        cat("Cool! All reads mapping on positive strand. No reads removed.\n")
      }
    }
    dt[annotation, on = "transcript", `:=`(c("cds_start",
                                             "cds_stop"), list(i.l_utr5 + 1, i.l_utr5 + i.l_cds))]
    dt[cds_start == 1 & cds_stop == 0, `:=`(cds_start, 0)]
    dt[, `:=`(strand, NULL)]
    if (granges == T || granges == TRUE) {
      dt <- GenomicRanges::makeGRangesFromDataFrame(dt,
                                                    keep.extra.columns = TRUE, ignore.strand = TRUE,
                                                    seqnames.field = c("transcript"), start.field = "end5",
                                                    end.field = "end3", strand.field = "strand",
                                                    starts.in.df.are.0based = FALSE)
      GenomicRanges::strand(dt) <- "+"
    }
    sample_reads_list[[sampname]] <- dt
  }
  if (granges == T || granges == TRUE) {
    sample_reads_list <- GenomicRanges::GRangesList(sample_reads_list)
  }
  return(sample_reads_list)
}



#' @title rlength_distr_rW
#' @description Function to plot ribo-seq read length distribution of a sample.
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}.
#' @param sample Sample name.
#' @param transcripts A vector containing the names of transcripts to be included. Default: NULL (all transcripts are included).
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths. Default: 99.
#' Use this argument to avoid printing out extremely uncommon read lengths. Default:99.
#' @details This function produces a read length histogram for the specified sample.
#' @examples
#' rlength_distr_rW(reads_list_LMCN, "CN34_r1_rpf")
#' @export
rlength_distr_rW <- function(reads_list, sample, transcripts = NULL, cl = 99){
  data <- reads_list
  if (length(transcripts) == 0) {
    dt <- data[[sample]]
  }
  else {
    dt <- data[[sample]][transcript %in% transcripts]
  }
  xmin <- quantile(dt$length, (1 - cl/100)/2)
  xmax <- quantile(dt$length, 1 - (1 - cl/100)/2)
  setkey(dt, length)
  dist <- dt[CJ(min(dt$length):max(dt$length)), list(count = .N),
             by = length][, `:=`(count, (count/sum(count)) * 100)]
  p <- ggplot(dist, aes(as.numeric(length), count)) + geom_bar(stat = "identity",
                                                               fill = "gray80") + labs(title = sample, x = "Read length",
                                                                                       y = "Count (%)") + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(xmin - 0.5, xmax + 0.5),
                       breaks = seq(xmin + ((xmin)%%2), xmax, by = max(c(1,
                                                                         floor((xmax - xmin)/7)))))
  output <- list()
  output[["plot"]] <- p
  output[["dt"]] <- dist
  return(output)
}



#' @title print_read_ldist
#' @description Function to print out read length distributions of all samples in a reads_list object to a pdf file.
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}.
#' @param outfile The path and name of the output pdf file.
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths.
#' Use this argument to avoid printing out extremely uncommon read lengths. Default: 99.
#' @details This function saves the read length histograms of all samples in a reads_list object to a pdf file.
#' @examples
#' print_read_ldist(reads_list_LMCN, "<file.path>/LMCN_RPF_Read_length_distributions.pdf")
#' @export
print_read_ldist <- function(reads_list, outfile, cl=99){
  pdf(outfile)
  for (sample_i in names(reads_list)){
    length_dist_zoom <- rlength_distr_rW(reads_list, sample=sample_i, cl=cl)
    print(length_dist_zoom[["plot"]])
  }
  dev.off()
}



#' @title rends_heat_rW
#' @description Function to print out read end heatmaps of a sample.
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}.
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param outfile The path and name of the output pdf file.
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths.
#' Use this argument to avoid printing out uncommon read lengths. Default:85.
#' @param utr5l Length of 5'UTR segment to be included upstream of translation start site. Default: 50.
#' @param cdsl Length of coding sequence to be included downstream of translation start site and upstream of translation termination site. Default: 50.
#' @param utr3l Length of 3'UTR segment to be included downstream of translation termination site. Default: 50.
#' @details Heatmaps of reads (stratified by read length) mapping to positions around the start and stop codon provide a visual
#' sense regarding a reasonable offset for p-site assignment as well as the 3-base periodicity of RPF reads.
#' @examples
#' print_read_end_heatmap(reads_list_LMCN, annotation_human_cDNA, "<file.path>/LMCN_RPF_read_end_heatmaps.pdf")
#' @export
rends_heat_rW <- function(reads_list, annotation, sample, transcripts = NULL, cl = 95,
                          utr5l = 50, cdsl = 50, utr3l = 50, log_colour = F, colour = "black"){

  data <- reads_list
  temp_dt <- data[[sample]]
  temp_dt[, `:=`(start_dist_end5, end5 - cds_start)][, `:=`(stop_dist_end5,
                                                            end5 - cds_stop)][, `:=`(start_dist_end3, end3 - cds_start)][,
                                                                                                                         `:=`(stop_dist_end3, end3 - cds_stop)]
  minlen <- ceiling(quantile(temp_dt$length, (1 - cl/100)/2))
  maxlen <- ceiling(quantile(temp_dt$length, 1 - (1 - cl/100)/2))
  l_transcripts <- as.character(annotation[l_utr5 >= utr5l &
                                             l_cds > 2 * (cdsl + 1) & l_utr3 >= utr3l, transcript])
  if (length(transcripts) == 0) {
    c_transcripts <- l_transcripts
  }
  else {
    c_transcripts <- intersect(l_transcripts, transcripts)
  }
  dt <- temp_dt[transcript %in% c_transcripts]
  temp_dt[, `:=`(c("start_dist_end5", "stop_dist_end5", "start_dist_end3",
                   "stop_dist_end3"), NULL)]
  start_sub <- dt[start_dist_end5 %in% seq(-utr5l, cdsl)]
  start_tab <- setkey(start_sub, length, start_dist_end5)[CJ(unique(length),
                                                             unique(start_dist_end5)), .N, by = .EACHI]
  setnames(start_tab, c("length", "dist", "count"))
  start_tab[, `:=`(region, "start")]
  stop_sub <- dt[stop_dist_end5 %in% seq(-cdsl, utr3l)]
  stop_tab <- setkey(stop_sub, length, stop_dist_end5)[CJ(unique(length),
                                                          unique(stop_dist_end5)), .N, by = .EACHI]
  setnames(stop_tab, c("length", "dist", "count"))
  stop_tab[, `:=`(region, "stop")]
  final_tab5 <- rbind(start_tab, stop_tab)
  final_tab5[, `:=`(end, "5end")]
  start_sub <- dt[start_dist_end3 %in% seq(-utr5l, cdsl)]
  start_tab <- setkey(start_sub, length, start_dist_end3)[CJ(unique(length),
                                                             unique(start_dist_end3)), .N, by = .EACHI]
  setnames(start_tab, c("length", "dist", "count"))
  start_tab[, `:=`(region, "start")]
  stop_sub <- dt[stop_dist_end3 %in% seq(-cdsl, utr3l)]
  stop_tab <- setkey(stop_sub, length, stop_dist_end3)[CJ(unique(length),
                                                          unique(stop_dist_end3)), .N, by = .EACHI]
  setnames(stop_tab, c("length", "dist", "count"))
  stop_tab[, `:=`(region, "stop")]
  final_tab3 <- rbind(start_tab, stop_tab)
  final_tab3[, `:=`(end, "3end")]
  final_tab <- rbind(final_tab5, final_tab3)
  final_tab[, `:=`(region, factor(region, levels = c("start",
                                                     "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)")))]
  final_tab[, `:=`(end, factor(end, levels = c("5end", "3end"),
                               labels = c("5' end", "3' end")))]
  max <- max(final_tab$count)
  p <- ggplot(final_tab, aes(dist, length)) + geom_tile(aes(fill = count)) +
    labs(title = paste(sample, "5' / 3' read end metaheatmaps",
                       sep = " - "), y = "Read length") + theme_bw(base_size = 20) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
          axis.title.x = element_blank()) + facet_grid(end ~
                                                         region, scales = "free", switch = "x") + theme(strip.background = element_blank(),
                                                                                                        strip.placement = "outside") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(minlen - 0.5, maxlen +
                                    0.5), breaks = seq(minlen + ((minlen)%%2), maxlen,
                                                       by = max(2, floor((maxlen - minlen)/7)))) + geom_vline(xintercept = 0,
                                                                                                              linetype = 2, color = "red")
  if (log_colour == F) {
    p <- p + scale_fill_gradient("Number\nof read\nextremities\n",
                                 low = "white", high = colour, limits = c(0.1, max),
                                 breaks = c(0.1, max/2, max), labels = c("0", floor(max/2),
                                                                         floor(max)), na.value = "white")
  }
  else {
    p <- p + scale_fill_gradient("Number\nof read\nextremities\n",
                                 low = "white", high = colour, limits = c(0.1, max),
                                 breaks = c(0.1, 10^(log10(max)/2 - 0.5), floor(max)),
                                 labels = c("0", floor(10^(log10(max)/2 - 0.5)), floor(max)),
                                 trans = "log", na.value = "transparent")
  }
  output <- list()
  output[["plot"]] <- p
  output[["dt"]] <- final_tab
  return(output)
}


#' @title print_read_end_heatmap
#' @description Function to print out metaheatmaps of reads around start and stop codons for all samples to a pdf file.
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param outfile The path and name of the output pdf file
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths.
#' Use this argument to avoid printing out uncommon read lengths. Default:85.
#' @param utr5l Length of 5'UTR segment to be included upstream of translation start site. Default: 50.
#' @param cdsl Length of coding sequence to be included downstream of translation start site and upstream of translation termination site. Default: 50.
#' @param utr3l Length of 3'UTR segment to be included downstream of translation termination site. Default: 50.
#' @details Heatmaps of reads (stratified by read length) mapping to positions around the start and stop codon provide a visual
#' sense regarding a reasonable offset for p-site assignment as well as the 3-base periodicity of RPF reads.
#' @examples
#' print_read_end_heatmap(reads_list_LMCN, annotation_human_cDNA, "<file.path>/LMCN_RPF_read_end_heatmaps.pdf")
#' @export
print_read_end_heatmap <- function(reads_list, annotation, outfile, cl=85, utr5l = 50, cdsl = 50, utr3l = 50){
  pdf(outfile, width=20, height=10)
  for (sample_i in names(reads_list)){
    ends_heatmap_i <- rends_heat_rW(reads_list, annotation, sample=sample_i,
                                    cl=cl, utr5l = utr5l, cdsl = cdsl, utr3l = utr3l)
    print(ends_heatmap_i[["plot"]])
  }
  dev.off()
}


#' @title psite_rW
#' @description Function to calculate the most likely position of p-site for each read length group.
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}
#' @param flanking Minimum number of nucleotides that have to flank the target codon in both directions. Default: 6.
#' @param start Logical argument indicating whether to use the start codon as reference. Default: TRUE.
#' If set to FALSE, the second to last codon is used.
#' @param extremity "5end", "3end" or "auto": the read end on which the correction step should be based. Default: "auto".
#' @param plot Whether read-length-specific ribosome occupancy plots are produced. Default: FALSE.
#' @param plot_dir The directory where the read-length-specific ribosome occupancy plots are saved. This argument is only
#' considered if \code{plot = TRUE}.
#' @param plot_format "png" or "pdf". It is only considered if \code{plot = TRUE}.
#' @return
#' The output is a data table, containing for all samples and read lengths the percentage of reads in the dataset,
#' the percentage of reads aligning to the start codon, and the distance of p-site from the two read extremities.
#' @examples
#' psite_offset_LMCN <- psite_rW(reads_list_LMCN)
#' @export
psite_rW <- function(reads_list, flanking = 6, start = TRUE, extremity = "auto",
                     plot = FALSE, plot_dir = NULL, plot_format = "png", cl = 99){
  data <- reads_list
  names <- names(data)
  offset <- NULL
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    lev <- sort(unique(dt$length))
    if (start == T | start == TRUE) {
      base <- 0
      dt[, `:=`(site_dist_end5, end5 - cds_start)]
      dt[, `:=`(site_dist_end3, end3 - cds_start)]
    }
    else {
      base <- -5
      dt[, `:=`(site_dist_end5, end5 - cds_stop - base)]
      dt[, `:=`(site_dist_end3, end3 - cds_stop - base)]
    }
    site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >=
                     flanking - 1]
    minlen <- min(site_sub$length)
    maxlen <- max(site_sub$length)
    t <- table(factor(site_sub$length, levels = lev))
    offset_temp <- data.table(length = as.numeric(as.character(names(t))),
                              percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    offset_temp[, `:=`(around_site, "T")][percentage == 0,
                                          `:=`(around_site, "F")]
    offset_temp5 <- site_sub[, list(offset_from_5 = as.numeric(names(which.max(table(site_dist_end5))))),
                             by = length]
    offset_temp3 <- site_sub[, list(offset_from_3 = as.numeric(names(which.max(table(site_dist_end3))))),
                             by = length]
    merge_allx <- function(x, y) merge(x, y, all.x = TRUE,
                                       by = "length")
    offset_temp <- Reduce(merge_allx, list(offset_temp, offset_temp5,
                                           offset_temp3))
    adj_off <- function(dt_site, dist_site, add, bestoff) {
      temp_v <- dt_site[[dist_site]]
      t <- table(factor(temp_v, levels = seq(min(temp_v) -
                                               2, max(temp_v) + add)))
      t[1:2] <- t[3] + 1
      locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) ==
                                                        -2)]))) + 1
      adjoff <- locmax[which.min(abs(locmax - bestoff))]
      ifelse(length(adjoff) != 0, adjoff, bestoff)
    }
    best_from5_tab <- offset_temp[, list(perc = sum(percentage)),
                                  offset_from_5][perc == max(perc)]
    best_from3_tab <- offset_temp[, list(perc = sum(percentage)),
                                  offset_from_3][perc == max(perc)]
    if (extremity == "auto" & ((best_from3_tab[, perc] >
                                best_from5_tab[, perc] & as.numeric(best_from3_tab[,
                                                                                   offset_from_3]) <= minlen - 2) | (best_from3_tab[,
                                                                                                                                    perc] <= best_from5_tab[, perc] & as.numeric(best_from5_tab[,
                                                                                                                                                                                                offset_from_5]) <= minlen - 1)) | extremity == "3end") {
      best_offset <- as.numeric(best_from3_tab[, offset_from_3])
      line_plot <- "from3"
      cat(sprintf("best offset: %i nts from the 3' end\n",
                  best_offset))
      adj_tab <- site_sub[, list(corrected_offset_from_3 = adj_off(.SD,
                                                                   "site_dist_end3", 0, best_offset)), by = length]
      offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE,
                           by = "length")
      offset_temp[is.na(corrected_offset_from_3), `:=`(corrected_offset_from_3,
                                                       best_offset)][, `:=`(corrected_offset_from_5,
                                                                            -corrected_offset_from_3 + length - 1)]
    }
    else {
      if (extremity == "auto" & ((best_from3_tab[, perc] <=
                                  best_from5_tab[, perc] & as.numeric(best_from5_tab[,
                                                                                     offset_from_5]) <= minlen - 1) | (best_from3_tab[,
                                                                                                                                      perc] > best_from5_tab[, perc] & as.numeric(best_from3_tab[,
                                                                                                                                                                                                 offset_from_3]) > minlen - 2)) | extremity ==
          "5end") {
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "from5"
        cat(sprintf("best offset: %i nts from the 5' end\n",
                    -best_offset))
        adj_tab <- site_sub[, list(corrected_offset_from_5 = adj_off(.SD,
                                                                     "site_dist_end5", 1, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE,
                             by = "length")
        offset_temp[is.na(corrected_offset_from_5), `:=`(corrected_offset_from_5,
                                                         best_offset)][, `:=`(corrected_offset_from_5,
                                                                              abs(best_offset))][, `:=`(corrected_offset_from_3,
                                                                                                        abs(corrected_offset_from_5 - length + 1))]
      }
    }
    t <- table(factor(dt$length, levels = lev))
    offset_temp[!is.na(offset_from_5), `:=`(offset_from_5,
                                            abs(offset_from_5))][, `:=`(total_percentage, as.numeric(format(round((as.vector(t)/sum(as.vector(t))) *
                                                                                                                    100, 3), nsmall = 4)))][, `:=`(percentage, as.numeric(format(round(percentage,
                                                                                                                                                                                       3), nsmall = 4)))][, `:=`(sample, n)]
    setcolorder(offset_temp, c("length", "total_percentage",
                               "percentage", "around_site", "offset_from_5", "offset_from_3",
                               "corrected_offset_from_5", "corrected_offset_from_3",
                               "sample"))
    if (start == TRUE | start == T) {
      setnames(offset_temp, c("length", "total_percentage",
                              "start_percentage", "around_start", "offset_from_5",
                              "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3",
                              "sample"))
    }
    else {
      setnames(offset_temp, c("length", "total_percentage",
                              "stop_percentage", "around_stop", "offset_from_5",
                              "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3",
                              "sample"))
    }
    if (plot == T | plot == TRUE) {
      options(warn = -1)
      if (length(plot_dir) == 0) {
        dir <- getwd()
        plot_dir <- paste(dir, "/offset_plot", sep = "")
      }
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
      }
      minlen <- ceiling(quantile(site_sub$length, (1 -
                                                     cl/100)/2))
      maxlen <- ceiling(quantile(site_sub$length, 1 - (1 -
                                                         cl/100)/2))
      for (len in minlen:maxlen) {
        progress <- ceiling(((len + 1 - minlen)/(maxlen -
                                                   minlen + 1)) * 25)
        cat(sprintf("\rplotting   %s\r", paste(paste(rep(c(" ",
                                                           "<<", "-"), c(25 - progress, 1, progress)),
                                                     collapse = ""), " ", as.character(progress *
                                                                                         4), "% ", paste(rep(c("-", ">>", " "), c(progress,
                                                                                                                                  1, 25 - progress)), collapse = ""), sep = "")))
        site_temp <- dt[site_dist_end5 %in% seq(-len +
                                                  1, 0) & length == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end5,
                                             levels = (-len + 1):(len))))
        site_temp <- dt[site_dist_end3 %in% seq(0, len -
                                                  2) & length == len]
        site_tab3 <- data.table(table(factor(site_temp$site_dist_end3,
                                             levels = (-len):(len - 2))))
        setnames(site_tab5, c("distance", "reads"))
        setnames(site_tab3, c("distance", "reads"))
        site_tab5[, `:=`(distance, as.numeric(as.character(site_tab5$distance)))][,
                                                                                  `:=`(extremity, "5' end")]
        site_tab3[, `:=`(distance, as.numeric(as.character(site_tab3$distance)))][,
                                                                                  `:=`(extremity, "3' end")]
        final_tab <- rbind(site_tab5[distance <= 0],
                           site_tab3[distance >= 0])
        final_tab[, `:=`(extremity, factor(extremity,
                                           levels = c("5' end", "3' end")))]
        p <- ggplot(final_tab, aes(distance, reads, color = extremity)) +
          geom_line() + geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) *
                                                      3, floor(max(final_tab$distance)/3) * 3, 3),
                                   linetype = 2, color = "gray90") + geom_vline(xintercept = 0,
                                                                                color = "gray50") + geom_vline(xintercept = -offset_temp[length ==
                                                                                                                                           len, offset_from_5], color = "#D55E00", linetype = 2,
                                                                                                               size = 1.1) + geom_vline(xintercept = offset_temp[length ==
                                                                                                                                                                   len, offset_from_3], color = "#56B4E9", linetype = 2,
                                                                                                                                        size = 1.1) + geom_vline(xintercept = -offset_temp[length ==
                                                                                                                                                                                             len, corrected_offset_from_5], color = "#D55E00",
                                                                                                                                                                 size = 1.1) + geom_vline(xintercept = offset_temp[length ==
                                                                                                                                                                                                                     len, corrected_offset_from_3], color = "#56B4E9",
                                                                                                                                                                                          size = 1.1) + annotate("rect", ymin = -Inf,
                                                                                                                                                                                                                 ymax = Inf, xmin = flanking - len, xmax = -flanking,
                                                                                                                                                                                                                 fill = "#D55E00", alpha = 0.1) + annotate("rect",
                                                                                                                                                                                                                                                           ymin = -Inf, ymax = Inf, xmin = flanking -
                                                                                                                                                                                                                                                             1, xmax = len - flanking - 1, fill = "#56B4E9",
                                                                                                                                                                                                                                                           alpha = 0.1) + labs(x = "Distance from start (nt)",
                                                                                                                                                                                                                                                                               y = "Number of read extremities", title = paste(n,
                                                                                                                                                                                                                                                                                                                               " - length=", len, " nts", sep = ""), color = "Extremity") +
          theme_bw(base_size = 20) + scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        if (line_plot == "from3") {
          p <- p + geom_vline(xintercept = best_offset,
                              color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset - len +
                         1, color = "black", linetype = 3, size = 1.1)
        }
        else {
          p <- p + geom_vline(xintercept = best_offset,
                              color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset + len -
                         1, color = "black", linetype = 3, size = 1.1)
        }
        p <- p + scale_x_continuous(limits = c(min(final_tab$distance),
                                               max(final_tab$distance)), breaks = seq(floor(min(final_tab$distance)/5) *
                                                                                        5, floor(max(final_tab$distance)/5) * 5, 5),
                                    labels = as.character(seq(floor(min(final_tab$distance)/5) *
                                                                5, floor(max(final_tab$distance)/5) * 5,
                                                              5) + base))
        subplot_dir <- paste(plot_dir, n, sep = "/")
        dir.create(subplot_dir)
        ggsave(paste(subplot_dir, "/", len, ".", plot_format,
                     sep = ""), plot = p, width = 15, height = 5,
               units = "in")
      }
      cat(sprintf("\rplotting   %s\n", paste(paste(rep(c(" ",
                                                         "<<", "-"), c(25 - progress, 1, progress)), collapse = ""),
                                             " ", as.character(progress * 4), "% ", paste(rep(c("-",
                                                                                                ">>", " "), c(progress, 1, 25 - progress)),
                                                                                          collapse = ""), sep = "")))
      options(warn = 0)
    }
    dt[, `:=`(c("site_dist_end5", "site_dist_end3"), NULL)]
    offset <- rbind(offset, offset_temp)
  }
  return(offset)
}



#' @title psite_info_rW
#' @description Function to add p-site offset information to a reads_list object
#' @param reads_list A reads_list object produced by \code{\link{bamtolist_rW}}
#' @param offset Offset data table produced by \code{\link{psite_rW}}
#' @param ... ...
#' @details This functions adds to each read four pieces of information regarding its p-site: distance from start
#' of the transcript, distance from start and end of the coding sequence (CDS), and region (5' UTR, CDS or 3' UTR).
#' @examples
#' reads_psite_list_LMCN <- psite_info_rW(reads_list_LMCN, psite_offset_LMCN)
#' @export
psite_info_rW <- function(reads_list, offset, site = NULL, fastapath = NULL, fasta_genome = TRUE,
                          refseq_sep = NULL, bsgenome = NULL, gtfpath = NULL, txdb = NULL,
                          dataSource = NA, organism = NA, granges = FALSE){
  data <- reads_list
  if (!(all(site %in% c("psite", "asite", "esite"))) & length(site) !=
      0) {
    cat("\n")
    stop("parameter site must be either NULL, \"psite\", \"asite\", \"esite\" or a combination of the three strings \n\n")
  }
  else {
    if (length(site) != 0 & length(fastapath) == 0 & length(bsgenome) ==
        0) {
      cat("\n")
      stop("parameter site is specified but both fastapath and bsgenome are missing \n\n")
    }
  }
  if (length(site) != 0) {
    if (((length(fastapath) != 0 & (fasta_genome == TRUE |
                                    fasta_genome == T)) | length(bsgenome) != 0) & length(gtfpath) ==
        0 & length(txdb) == 0) {
      cat("\n")
      stop("genome annotation file not specified (both GTF path and TxDb object are missing)\n\n")
    }
    if (length(fastapath) != 0 & length(bsgenome) != 0) {
      cat("\n")
      warning("both fastapath and bsgenome are specified. Only fastapath will be considered\n")
      bsgenome = NULL
    }
    if (length(gtfpath) != 0 & length(txdb) != 0) {
      cat("\n")
      warning("both gtfpath and txdb are specified. Only gtfpath will be considered\n")
      txdb = NULL
    }
    if ((length(gtfpath) != 0 | length(txdb) != 0) & ((length(fastapath) ==
                                                       0 & length(bsgenome) == 0) | (length(fastapath) !=
                                                                                     0 & (fasta_genome == FALSE | fasta_genome == F)))) {
      cat("\n")
      warning("a genome annotation file is specified but no sequences from genome assembly are provided\n")
    }
    if (length(gtfpath) != 0 | length(txdb) != 0) {
      if (length(gtfpath) != 0) {
        path_to_gtf <- gtfpath
        txdbanno <- GenomicFeatures::makeTxDbFromGFF(file = path_to_gtf,
                                                     format = "gtf", dataSource = dataSource, organism = organism)
      }
      else {
        if (txdb %in% rownames(installed.packages())) {
          library(txdb, character.only = TRUE)
        }
        else {
          source("https://bioconductor.org/biocLite.R")
          biocLite(txdb, suppressUpdates = TRUE)
          library(txdb, character.only = TRUE)
        }
        txdbanno <- get(txdb)
      }
    }
    if (length(fastapath) != 0 | length(bsgenome) != 0) {
      if (length(fastapath) != 0) {
        if (fasta_genome == TRUE | fasta_genome == T) {
          temp_sequences <- Biostrings::readDNAStringSet(fastapath,
                                                         format = "fasta", use.names = TRUE)
          if (length(refseq_sep) != 0) {
            names(temp_sequences) <- tstrsplit(names(temp_sequences),
                                               refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
          exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno,
                                                            by = "tx", use.names = TRUE))
          exon <- as.data.table(exon[unique(names(exon))])
          sub_exon_plus <- exon[as.character(seqnames) %in%
                                  names(temp_sequences) & strand == "+"]
          sub_exon_minus <- exon[as.character(seqnames) %in%
                                   names(temp_sequences) & strand == "-"][,
                                                                          `:=`(new_end, Biostrings::width(temp_sequences[as.character(seqnames)]) -
                                                                                 start + 1)][, `:=`(new_start, Biostrings::width(temp_sequences[as.character(seqnames)]) -
                                                                                                      end + 1)]
          seq_dt_plus <- sub_exon_plus[, `:=`(nt_seq,
                                              "emp")][, `:=`(nt_seq, as.character(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                                                                     start = start, end = end)))][, list(seq = paste(nt_seq,
                                                                                                                                                     collapse = "")), by = group_name]
          revcompl_temp_sequences <- reverseComplement(temp_sequences)
          seq_dt_minus <- sub_exon_minus[, `:=`(nt_seq,
                                                "emp")][, `:=`(nt_seq, as.character(Biostrings::subseq(revcompl_temp_sequences[as.character(seqnames)],
                                                                                                       start = new_start, end = new_end)))][, list(seq = paste(nt_seq,
                                                                                                                                                               collapse = "")), by = group_name]
          sequences <- Biostrings::DNAStringSet(c(seq_dt_plus$seq,
                                                  seq_dt_minus$seq))
          names(sequences) <- c(unique(sub_exon_plus$group_name),
                                unique(sub_exon_minus$group_name))
        }
        else {
          sequences <- Biostrings::readDNAStringSet(fastapath,
                                                    format = "fasta", use.names = TRUE)
          if (length(refseq_sep) != 0) {
            names(sequences) <- tstrsplit(names(sequences),
                                          refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
        }
      }
      else {
        if (bsgenome %in% installed.genomes()) {
          library(bsgenome, character.only = TRUE)
        }
        else {
          source("http://www.bioconductor.org/biocLite.R")
          biocLite(bsgenome, suppressUpdates = TRUE)
          library(bsgenome, character.only = TRUE)
        }
        sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome),
                                                            txdbanno, use.names = T)
      }
    }
  }
  names <- names(data)
  for (n in names) {
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    suboff <- offset[sample == n, .(length, corrected_offset_from_3)]
    cat("1. adding p-site position\n")
    dt[suboff, on = "length", `:=`(psite, i.corrected_offset_from_3)]
    dt[, `:=`(psite, end3 - psite)]
    setcolorder(dt, c("transcript", "end5", "psite", "end3",
                      "length", "cds_start", "cds_stop"))
    dt[, `:=`(psite_from_start, psite - cds_start)][cds_stop ==
                                                      0, `:=`(psite_from_start, 0)]
    dt[, `:=`(psite_from_stop, psite - cds_stop)][cds_stop ==
                                                    0, `:=`(psite_from_stop, 0)]
    cat("2. adding transcript region\n")
    dt[, `:=`(psite_region, "5utr")][psite_from_start >=
                                       0 & psite_from_stop <= 0, `:=`(psite_region, "cds")][psite_from_stop >
                                                                                              0, `:=`(psite_region, "3utr")][cds_stop == 0, `:=`(psite_region,
                                                                                                                                                 NA)]
    if (length(site) != 0) {
      cat("3. adding nucleotide sequence(s)\n")
      if ("psite" %in% site) {
        dt[, `:=`(p_site_codon, as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                                start = dt$psite, end = dt$psite + 2)))]
      }
      if ("asite" %in% site) {
        dt[, `:=`(a_site_codon, as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                                start = dt$psite + 3, end = dt$psite + 5)))]
      }
      if ("esite" %in% site) {
        dt[, `:=`(e_site_codon, as.character(Biostrings::subseq(sequences[as.character(dt$transcript)],
                                                                start = dt$psite - 3, end = dt$psite - 1)))]
      }
    }
    setorder(dt, transcript, end5)
    if (granges == T | granges == TRUE) {
      dt <- GenomicRanges::makeGRangesFromDataFrame(dt,
                                                    keep.extra.columns = TRUE, ignore.strand = TRUE,
                                                    seqnames.field = c("transcript"), start.field = "end5",
                                                    end.field = "end3", strand.field = "strand",
                                                    starts.in.df.are.0based = FALSE)
      GenomicRanges::strand(dt) <- "+"
    }
    data[[n]] <- dt
  }
  if (granges == T | granges == TRUE) {
    data <- GenomicRanges::GRangesList(data)
  }
  return(data)
}


#' @title metaprofile_psite_rW
#' @description Function to plot read metaprofiles around the start and stop codons for a sample
#' @param reads_psite_list A list of reads and their psite coordinates produced by \code{\link{psite_info_rW}}
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param sample Sample to be plotted.
#' @param length_range Range of read lengths to be included. Default: "all".
#' @param transcripts Vector of transcript names to be included. Default: NULL (includes all transcripts).
#' @param utr5l Length of the 5'UTR upstream of start codon to be included in the plot. Default: 25.
#' @param cdsl Length of the cds downstream of start codon and upstream od stop codon to be included in the plot. Default: 50.
#' @param utr3l Length of the 3'UTR downstream of stop codon to be included in the plot. Default: 25.
#' @param plot_title Title of the plot. Default: NULL.
#' @details Read metaprofiles (ribosome occupancy plots) visualize the positional distribution of reads mapping around start and stop codons.
#' @examples
#' print_rop(reads_psite_list_LMCN, annotation_human_cDNA, "<file.path>/LMCN_RPF_ribosome_occupancy_profiles_from_annotation.pdf")
#' @export
metaprofile_psite_rW <- function(reads_psite_list, annotation, sample, scale_factors = NULL, length_range = "all",
                                 transcripts = NULL, utr5l = 25, cdsl = 50, utr3l = 25, plot_title = NULL){
  data <- reads_psite_list
  if (!identical(length_range, "all") & !inherits(length_range,
                                                  "numeric") & !inherits(length_range, "integer")) {
    cat("\n")
    warning("class of length_range is neither numeric nor integer. Set to default \"all\"\n")
    length_range = "all"
  }
  if (!identical(length_range, "all")) {
    for (samp in sample) {
      len_check <- unique(data[[samp]]$length)
      if (sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the specified lengths: sample removed\n",
                        samp))
        sample <- sample[sample != samp]
      }
    }
  }
  if (length(sample) == 0) {
    cat("\n")
    stop("none of the data tables in sample contains any reads of the specified lengths\n\n")
  }
  if (length(scale_factors) != 0) {
    if (!all(sample %in% names(scale_factors))) {
      cat("\n")
      stop("scale factor for one or more replicates is missing\n\n")
    }
  }
  l_transcripts <- as.character(annotation[l_utr5 >= utr5l &
                                             l_cds >= 2 * (cdsl + 1) & l_utr3 >= utr3l, transcript])
  if (length(transcripts) == 0) {
    c_transcripts <- l_transcripts
    ntr <- length(c_transcripts)
  }
  else {
    c_transcripts <- intersect(l_transcripts, transcripts)
    ntr <- length(transcripts)
  }
  length_temp <- vector()
  for (samp in sample) {
    dt <- data[[samp]][as.character(transcript) %in% c_transcripts,
                       ]
    if (identical(length_range, "all")) {
      start_sub <- dt[psite_from_start %in% seq(-utr5l,
                                                cdsl)]
      stop_sub <- dt[psite_from_stop %in% seq(-cdsl, utr3l)]
    }
    else {
      start_sub <- dt[psite_from_start %in% seq(-utr5l,
                                                cdsl) & length %in% length_range]
      stop_sub <- dt[psite_from_stop %in% seq(-cdsl, utr3l) &
                       length %in% length_range]
    }
    setkey(start_sub, psite_from_start)
    start_tab <- start_sub[CJ(-utr5l:cdsl), list(reads = .N),
                           by = list(distance = psite_from_start)][, `:=`(reg,
                                                                          "start")]
    setkey(stop_sub, psite_from_stop)
    stop_tab <- stop_sub[CJ(-cdsl:utr3l), list(reads = .N),
                         by = list(distance = psite_from_stop)][, `:=`(reg,
                                                                       "stop")]
    samp_tab <- rbind(start_tab, stop_tab)
    if (length(scale_factors) != 0) {
      samp_tab[, `:=`(reads, reads * scale_factors[samp])]
    }
    if (exists("final_tab_psm")) {
      final_tab_psm[, `:=`(reads, reads + samp_tab$reads)]
    }
    else {
      final_tab_psm <- samp_tab
    }
    length_temp <- unique(c(length_temp, data[[samp]]$length))
  }
  if (!identical(length_range, "all")) {
    length_range <- sort(intersect(length_temp, length_range))
  }
  else {
    length_range <- sort(length_temp)
  }
  final_tab_psm[, `:=`(reg, factor(reg, levels = c("start",
                                                   "stop"), labels = c("Distance from start (nt)", "Distance from stop (nt)")))]
  linestart <- data.table(reg = rep(c("Distance from start (nt)",
                                      "Distance from stop (nt)"), times = c(length(c(rev(seq(-3,
                                                                                             -utr5l, -3)), seq(3, cdsl, 3))), length(c(rev(seq(-2,
                                                                                                                                               -cdsl, -3)), seq(1, utr3l, 3))))), line = c(rev(seq(-3,
                                                                                                                                                                                                   -utr5l, -3)), seq(3, cdsl, 3), rev(seq(-2, -cdsl, -3)),
                                                                                                                                                                                           seq(1, utr3l, 3)))
  linered <- data.table(reg = c("Distance from start (nt)",
                                "Distance from stop (nt)"), line = c(0, 1))
  plot <- ggplot(final_tab_psm, aes(distance, reads)) + geom_line(size = 1.05,
                                                                  color = "gray40") + geom_vline(data = linered, aes(xintercept = line),
                                                                                                 linetype = 1, color = "red") + labs(x = "", y = "P-site") +
    theme_bw(base_size = 20) + theme(panel.grid.major.x = element_blank(),
                                     panel.grid.minor.x = element_blank()) + facet_grid(. ~
                                                                                          reg, scales = "free", switch = "x") + theme(strip.background = element_blank(),
                                                                                                                                      strip.placement = "outside") + geom_vline(data = linestart,
                                                                                                                                                                                aes(xintercept = line), linetype = 3, color = "gray60")
  if (identical(plot_title, "auto")) {
    title1 <- paste0(paste(sample, collapse = "+"), " (",
                     ntr, " tr). Read length: ")
    minlr <- min(length_range)
    maxlr <- max(length_range)
    if (minlr == maxlr) {
      plottitle <- paste0(title1, min(length_range), " nts")
    }
    else {
      if (identical(length_range, minlr:maxlr) | identical(length_range,
                                                           seq(minlr, maxlr, 1))) {
        plottitle <- paste0(title1, minlr, "-", maxlr,
                            " nts")
      }
      else {
        nextl <- sort(length_range[c(which(diff(length_range) !=
                                             1), which(diff(length_range) != 1) + 1)])
        sep <- ifelse(nextl %in% length_range[which(diff(length_range) !=
                                                      1)], ", ", "-")[-length(nextl)]
        if (1 %in% which(diff(length_range) == 1)) {
          nextl <- c(length_range[1], nextl)
          sep <- c("-", sep)
        }
        if ((length(length_range) - 1) %in% which(diff(length_range) ==
                                                  1)) {
          nextl <- c(nextl, length_range[length(length_range)])
          sep <- c(sep, "-")
        }
        sep <- c(sep, "")
        plottitle <- paste0(title1, paste0(nextl, sep,
                                           collapse = ""), " nts")
      }
    }
    plot <- plot + labs(title = plottitle) + theme(plot.title = element_text(hjust = 0.5))
  }
  else {
    if (length(plot_title) != 0) {
      plot <- plot + labs(title = plot_title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  output <- list()
  output[["plot"]] <- plot
  output[["dt"]] <- final_tab_psm
  return(output)
}



#' @title print_rop
#' @description Function to print out ribosome occupancy profiles of all samples to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param outfile The path and name of the output pdf file
#' @details Ribosome occupancy profiles visualize the positional distribution of reads mapping around start and stop codons.
#' @examples
#' print_rop(LMCN_reads_psite_list, annotation_human_cDNA, "<file.path>/LMCN_RPF_ribosome_occupancy_profiles.pdf")
#' @export
print_rop <- function(reads_psite_list, annotation, outfile){
  pdf(outfile, width=20, height=10)
  for (sample_i in names(reads_psite_list)){
    metaprofile_psite_sample_i <- metaprofile_psite_rW(reads_psite_list, annotation, sample = sample_i, plot_title = sample_i)
    print(metaprofile_psite_sample_i[["plot"]])
  }
  dev.off()
}



#' @title frame_psite_rW
#' @description Function to plot the percentage of psites falling into each of the three reading frames (periodicity)
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param sample Sample to be plotted.
#' @param length_range Range of read lengths to be included. Default: "all".
#' @param transcripts Vector of transcript names to be included. Default: NULL (includes all transcripts).
#' @param region "5utr", "cds", "3utr" or "all". Default: "all".
#' @param plot_title Title of the plot. Default: NULL.
#' @details This function produces bar plots of the percentage of psites falling into each of the three possible reading frames,
#' separately for 5' UTR, CDS and 3' UTR regions. In a good dataset, psites in CDS should be enriched for frame 1.
#' @examples
#' frame_psite_rW(reads_psite_list_LMCN, "CN34_r2_rpf")
#' @export
frame_psite_rW <- function (reads_psite_list, sample = NULL, transcripts = NULL, region = "all",
                            length_range = "all", plot_title = NULL){
  data <- reads_psite_list
  if (length(sample) == 0) {
    sample <- names(data)
  }
  if (!identical(length_range, "all") & !inherits(length_range,
                                                  "numeric") & !inherits(length_range, "integer")) {
    cat("\n")
    warning("class of length_range is neither numeric nor integer. Set to default \"all\"\n")
    length_range = "all"
  }
  if (!identical(length_range, "all")) {
    for (samp in sample) {
      if (length(transcripts) == 0) {
        dt <- data[[samp]]
      }
      else {
        dt <- data[[samp]][transcript %in% transcripts]
      }
      len_check <- unique(dt$length)
      if (sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the specified lengths: sample removed\n",
                        samp))
        sample <- sample[sample != samp]
      }
    }
  }
  if (length(sample) == 0) {
    cat("\n")
    stop("none of the data tables in sample contains any reads of the specified lengths\n\n")
  }
  if (!region %in% c("all", "cds", "5utr", "3utr")) {
    cat("\n")
    warning("region is invalid. Set to default \"all\"\n")
    region = "all"
  }
  length_temp <- vector()
  for (samp in sample) {
    if (identical(length_range, "all")) {
      dt <- data[[samp]]
    }
    else {
      dt <- data[[samp]][length %in% length_range]
    }
    if (length(transcripts) != 0) {
      dt <- dt[transcript %in% transcripts]
    }
    if (region == "all") {
      frame_dt <- dt[cds_start != 0 & cds_stop != 0][,
                                                     `:=`(frame, psite_from_start%%3)][, list(count = .N),
                                                                                       by = list(region = psite_region, frame)][, `:=`(percentage,
                                                                                                                                       (count/sum(count)) * 100), by = region][is.na(percentage),
                                                                                                                                                                               `:=`(percentage, 0)][, `:=`(region, factor(region,
                                                                                                                                                                                                                          levels = c("5utr", "cds", "3utr"), labels = c("5' UTR",
                                                                                                                                                                                                                                                                        "CDS", "3' UTR")))]
    }
    else {
      frame_dt <- dt[psite_region == region][cds_start !=
                                               0 & cds_stop != 0][, `:=`(frame, psite_from_start%%3)][,
                                                                                                      list(count = .N), by = frame][, `:=`(percentage,
                                                                                                                                           (count/sum(count)) * 100)][is.na(percentage),
                                                                                                                                                                      `:=`(percentage, 0)]
    }
    frame_dt[, `:=`(sample, samp)]
    if (exists("final_frame_dt")) {
      final_frame_dt <- rbind(final_frame_dt, frame_dt)
    }
    else {
      final_frame_dt <- frame_dt
    }
    length_temp <- unique(c(length_temp, data[[samp]]$length))
  }
  if (!identical(length_range, "all")) {
    length_range <- sort(intersect(length_range, length_temp))
  }
  else {
    length_range <- sort(length_temp)
  }
  plot <- ggplot(final_frame_dt, aes(x = frame, y = percentage)) +
    geom_bar(stat = "identity") + theme_bw(base_size = 20) +
    labs(x = "Frame", y = "P-site signal (%)")
  if (region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(sample, region,
                                           frame)]
  }
  else {
    plot <- plot + facet_wrap(~sample, ncol = 3)
    final_frame_dt <- final_frame_dt[order(sample, frame)]
  }
  if (identical(plot_title, "auto")) {
    if (region == "all") {
      plottitle_region <- NULL
    }
    else {
      if (region == "5utr") {
        plottitle_region <- "Region: 5' UTR. "
      }
      if (region == "cds") {
        plottitle_region <- "Region: CDS. "
      }
      if (region == "3utr") {
        plottitle_region <- "Region: 3' UTR. "
      }
    }
    minlr <- min(length_range)
    maxlr <- max(length_range)
    if (minlr == maxlr) {
      plottitle_range <- paste0("Read length: ", minlr,
                                " nts")
    }
    else {
      if (identical(length_range, minlr:maxlr) | identical(length_range,
                                                           seq(minlr, maxlr, 1))) {
        plottitle_range <- paste0("Read lengths: ", minlr,
                                  "-", maxlr, " nts")
      }
      else {
        nextl <- sort(length_range[c(which(diff(length_range) !=
                                             1), which(diff(length_range) != 1) + 1)])
        sep <- ifelse(nextl %in% length_range[which(diff(length_range) !=
                                                      1)], ", ", "-")[-length(nextl)]
        if (1 %in% which(diff(length_range) == 1)) {
          nextl <- c(length_range[1], nextl)
          sep <- c("-", sep)
        }
        if ((length(length_range) - 1) %in% which(diff(length_range) ==
                                                  1)) {
          nextl <- c(nextl, length_range[length(length_range)])
          sep <- c(sep, "-")
        }
        sep <- c(sep, "")
        plottitle_range <- paste0("Read lengths: ", paste0(nextl,
                                                           sep, collapse = ""), " nts")
      }
    }
    plottitle <- paste0(plottitle_region, plottitle_range)
    plot <- plot + labs(title = plottitle) + theme(plot.title = element_text(hjust = 0.5))
  }
  else {
    if (length(plot_title) != 0) {
      plot <- plot + labs(title = plot_title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  ret_list <- list()
  ret_list[["dt"]] <- final_frame_dt
  ret_list[["plot"]] <- plot
  return(ret_list)
}



#' @title print_period_region
#' @description Function to print out psite periodicity by region plots of all samples to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param outfile The path and name of the output pdf file
#' @details These plots illustrate the 3-base (reading frame-dependent) periodicity of reads, separately for 5'UTR, CDS and 3'UTR regions.
#' In a good dataset, psites in the CDS should be enriched for frame 1.
#' @examples
#' print_period_region(reads_psite_list_LMCN, "<file.path>/LMCN_Periodicity_by_region.pdf")
#' @export
print_period_region <- function(reads_psite_list, outfile){
  pdf(outfile)
  for (sample_i in names(reads_psite_list)){
    frames_i <- frame_psite_rW(reads_psite_list, sample=sample_i, region="all")
    print(frames_i[["plot"]])
  }
  dev.off()
}



#' @title frame_psite_length_rW
#' @description Function to plot the percentage of psites falling into each of the three reading frames (periodicity) stratified by length
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param sample Sample to be plotted.
#' @param length_range Range of read lengths to be included. Default: "all".
#' @param transcripts Vector of transcript names to be included. Default: NULL (includes all transcripts).
#' @param region "5utr", "cds", "3utr" or "all". Default: "all".
#' @param plot_title Title of the plot. Default: NULL.
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths.
#' Use this argument to avoid printing out uncommon read lengths. Default:95.
#' @examples
#' print_period_region_length(reads_psite_list_LMCN, "<file.path>/LMCN_Periodicity_by_length_region.pdf")
#' @export

frame_psite_length_rW <- function (reads_psite_list, sample = NULL, transcripts = NULL, region = "all",
                                   cl = 95, length_range = "all", plot_title = NULL){
  data <- reads_psite_list
  if (length(sample) == 0) {
    sample <- names(data)
  }
  if (!identical(length_range, "all") & !inherits(length_range,
                                                  "numeric") & !inherits(length_range, "integer")) {
    cat("\n")
    warning("class of length_range is neither numeric nor integer. Set to default \"all\"\n")
    length_range = "all"
  }
  if (!identical(length_range, "all")) {
    for (samp in sample) {
      if (length(transcripts) == 0) {
        dt <- data[[samp]]
      }
      else {
        dt <- data[[samp]][transcript %in% transcripts]
      }
      len_check <- unique(dt$length)
      if (sum(length_range %in% len_check) == 0) {
        cat("\n")
        warning(sprintf("\"%s\" doesn't contain any reads of the specified lengths: sample removed\n",
                        samp))
        sample <- sample[sample != samp]
      }
    }
  }
  if (length(sample) == 0) {
    cat("\n")
    stop("none of the data tables in sample contains any reads of the specified lengths\n\n")
  }
  if (!identical(length_range, "all")) {
    minl <- min(length_range)
    maxl <- max(length_range)
  }
  if (!region %in% c("all", "cds", "5utr", "3utr")) {
    cat("\n")
    warning("region is invalid. Set to default \"all\"\n")
    region = "all"
  }
  for (samp in sample) {
    if (length(transcripts) == 0) {
      dt <- data[[samp]]
    }
    else {
      dt <- data[[samp]][transcript %in% transcripts]
    }
    if (region == "all") {
      dt <- dt[cds_start != 0 & cds_stop != 0][, `:=`(frame,
                                                      psite_from_start%%3)]
      if (identical(length_range, "all")) {
        minl <- quantile(dt$length, (1 - cl/100)/2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100)/2)
        length_range <- minl:maxl
      }
      dt[, `:=`(psite_region, factor(psite_region, levels = c("5utr",
                                                              "cds", "3utr")))][, `:=`(frame, factor(frame,
                                                                                                     levels = c(0, 1, 2)))][, `:=`(length, factor(length,
                                                                                                                                                  levels = unique(length)))]
      setkey(dt, length, psite_region, frame)
      frame_dt <- dt[CJ(levels(length), levels(psite_region),
                        levels(frame)), list(count = .N), by = .EACHI][as.numeric(as.character(length)) %in%
                                                                         length_range][, `:=`(percentage, (count/sum(count)) *
                                                                                                100), by = list(length, psite_region)][is.na(percentage),
                                                                                                                                       `:=`(percentage, 0)][, `:=`(psite_region, factor(psite_region,
                                                                                                                                                                                        levels = c("5utr", "cds", "3utr"), labels = c("5' UTR",
                                                                                                                                                                                                                                      "CDS", "3' UTR")))]
      setnames(frame_dt, "psite_region", "region")
    }
    else {
      dt <- dt[psite_region == region][, `:=`(frame, psite_from_start%%3)]
      if (identical(length_range, "all")) {
        minl <- quantile(dt$length, (1 - cl/100)/2)
        maxl <- quantile(dt$length, 1 - (1 - cl/100)/2)
        length_range <- minl:maxl
      }
      dt[, `:=`(frame, factor(frame, levels = c(0, 1, 2)))][,
                                                            `:=`(length, factor(length, levels = unique(length)))]
      setkey(dt, length, frame)
      frame_dt <- dt[CJ(levels(length), levels(frame)),
                     list(count = .N), by = .EACHI][as.numeric(as.character(length)) %in%
                                                      length_range][, `:=`(percentage, (count/sum(count)) *
                                                                             100), by = length][is.na(percentage), `:=`(percentage,
                                                                                                                        0)]
    }
    frame_dt$sample <- samp
    if (exists("final_frame_dt")) {
      final_frame_dt <- rbind(final_frame_dt, frame_dt)
    }
    else {
      final_frame_dt <- frame_dt
    }
  }
  mins <- min(final_frame_dt$percentage)
  maxs <- max(final_frame_dt$percentage)
  plot <- ggplot(final_frame_dt, aes(frame, as.numeric(as.character(length)))) +
    geom_tile(aes(fill = percentage)) + scale_fill_gradient("P-site signal (%)  ",
                                                            low = "white", high = "#104ec1", breaks = c(mins, mins/2 +
                                                                                                          maxs/2, maxs), labels = c(round(mins), round(mins/2 +
                                                                                                                                                         maxs/2), round(maxs))) + labs(x = "Frame", y = "Read length") +
    theme_bw(base_size = 20) + theme(panel.grid.major.x = element_blank(),
                                     panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
                                     panel.grid.minor.y = element_blank()) + theme(legend.position = "top",
                                                                                   legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(5,
                                                                                                                                                  0, -5, 0)) + scale_y_continuous(limits = c(minl -
                                                                                                                                                                                               0.5, maxl + 0.5), breaks = seq(minl + ((minl)%%2), maxl,
                                                                                                                                                                                                                              by = max(2, floor((maxl - minl)/7))))
  if (region == "all") {
    plot <- plot + facet_grid(sample ~ region)
    final_frame_dt <- final_frame_dt[order(sample, region,
                                           length, frame)]
  }
  else {
    plot <- plot + facet_wrap(~sample, ncol = 3)
    final_frame_dt <- final_frame_dt[order(sample, length,
                                           frame)]
  }
  if (identical(plot_title, "auto") & !identical(region, "all")) {
    if (region == "5utr") {
      plottitle <- "Region: 5' UTR"
    }
    if (region == "cds") {
      plottitle <- "Region: CDS"
    }
    if (region == "3utr") {
      plottitle <- "Region: 3' UTR"
    }
    plot <- plot + labs(title = plottitle) + theme(plot.title = element_text(hjust = 0.5))
  }
  else {
    if (length(plot_title) != 0 & !identical(plot_title,
                                             "auto")) {
      plot <- plot + labs(title = plot_title) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  ret_list <- list()
  ret_list[["dt"]] <- final_frame_dt
  ret_list[["plot"]] <- plot
  return(ret_list)
}


#' @title print_period_region_length
#' @description Function to print out periodicity by length and region plots of all samples to a pdf file
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param outfile The path and name of the output pdf file
#' @param cl An integer in [1,100] specifying a confidence level to restrict the plot to a sub-range of read lengths.
#' Use this argument to avoid printing out uncommon read lengths. Default:95.
#' @details These plots illustrate the 3-base (reading-frame-dependent) periodicity of reads, separately for 5'UTR, CDS and 3'UTR and stratified by read length.
#' They can inform the choice of read length range for further analysis of ribo-seq data. Good-length reads are abundant and show enrichment for frame 1 in the CDS region.
#' @examples
#' print_period_region_length(reads_psite_list, "<file.path>/Periodicity_by_length_region2.pdf")
#' @export
print_period_region_length <- function(reads_psite_list, outfile, cl=95){
  pdf(outfile, width=20, height=15)
  for (sample_i in names(reads_psite_list)){
    frames_stratified_i <- frame_psite_length_rW(reads_psite_list, sample=sample_i, region="all", cl=cl)
    print(frames_stratified_i[["plot"]])
  }
  dev.off()
}



#' @title psite_to_codon_count
#' @description Function to count reads per codon for each transcript.
#' @param reads_psite_list A reads_psite_list object produced by \code{\link{psite_info_rW}}
#' @param length_range Range of read lengths to be included in calculation of per codon RPF counts.
#' Longer and shorter reads may be excluded because of high uncertainty in offset assignment, lack of expected periodicity, etc.
#' @param annotation Annotation data table produced by \code{\link{read_annotation}} listing transcript names and lengths of their 5'UTR, CDS and 3'UTR segments.
#' It has five columns: transcript, l_tr, l_utr5, l_cds and l_utr3.
#' Transcript names and segment lengths must correspond to the reference sequences to which the reads were mapped.
#' @param fasta.file The reference fasta file to which reads were mapped.
#' @examples
#' tr_codon_read_count_LMCN <- psite_to_codon_count(reads_psite_list_LMCN, c(24:32), annotation_human_cDNA, "<file.path>/Human.GRC38.96_cDNA_longest_CDS.txt")
#' @return A list of lists with the following structure: list$<sample.name>$<transcript.ID> data.frame: [1] codon_number [2] codon_type [3] aa_type [4] observed_count
#' @export
psite_to_codon_count <- function(reads_psite_list, length_range, annotation, fasta.file){
  # Filter reads_psite_list to keep reads within the specified length range and mapping to CDS
  reads_psite_list.m <- lapply(reads_psite_list, function(x) subset(x, psite_region == "cds" & length %in% length_range))

  # Remove transcripts with CDS length non-divisible by 3.
  annotation.df.m <- as.data.frame(annotation)
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



#' @title CELP_detect_bias
#' @description Function to compute codon-level stalling bias coefficients using the CELP (Consistent Excess of Loess Preds) method
#' @param tr_codon_read_count_list A list generated by \code{\link{psite_to_codon_count}}
#' containing observed read counts per codon per trancript per sample.
#' @param codon_radius Number of codons on either side of each codon influencing the loess prediction for the middle codon. Default: 5.
#' @param loess_method Determines whether the fitted surface should be computed exactly ("direct") or via interpolation from a kd tree ("interpolate").
#' A third option is "none" which means than raw observed counts to calculate the bias coefficient (no loess). Default: "interpolate".
#' @details This function starts with running a loess curve on per-codon read counts along the transcript to borrow information
#' from neighboring codons mitigating the uncertainty of p-site offset assignment and experimental stochasticity.
#' Loess span parameter is calculated from the user-defined codon radius and CDS length.
#' Then, a bias coefficient is calculated for each codon by integrating information on the excess of loess-predicted
#' read counts at that codon comapred to the transcript's background across all samples. Finally, loess predicted counts
#' read counts are divided by bias coefficients to calculate bias-corrected counts.
#' The "direct" fitting method takes longer to complete but does not run into kd-tree-related memory issues.
#' The loess_method "none" option is included mainly to enable comparison and show effectiveness of using loess in
#' finding bona fide bias positions. It is NOT recommended for actual bias detection.
#' @return A list of data frames (one per transcript) containing codon-level bias coefficients.
#' It has the following structure: list$<transcript.ID> data.frame: [1] codon_number [2] codon_type [3] aa_type [4] bias_coefficient.
#' @examples
#' bias_coeff_list_LMCN_noloess <- CELP_detect_bias(tr_codon_read_count_LMCN, loess_method = "none")
#' bias_coeff_list_LMCN_withloess <- CELP_detect_bias(tr_codon_read_count_LMCN, loess_method = "interpolate")
#' @export
CELP_detect_bias <- function(tr_codon_read_count_list, loess_method = "interpolate", codon_raduis = 5){

  x <- tr_codon_read_count_list
  bias_coefficients_list <- list()
  sample_names_i <- names(x)
  tr_names_i <- names(x[[1]])

  # Run loess and compute loess predicted values
  for (t in tr_names_i){
    if (loess_method == "none"){
      for (s in sample_names_i){
        x[[s]][[t]]$excess_ratio <- x[[s]][[t]]$observed_count / median(x[[s]][[t]]$observed_count[x[[s]][[t]]$observed_count>0])
      }
    } else {
      l_cds <- dim(x[[1]][[t]])[1]
      span_tr <- (2 * codon_raduis + 1) / l_cds
      for (s in sample_names_i){
        x[[s]][[t]]$loess_pred <- suppressWarnings(predict(loess(x[[s]][[t]]$observed_count ~ x[[s]][[t]]$codon_number, span = span_tr, se = FALSE, control = loess.control(surface = loess_method))))
        x[[s]][[t]]$loess_pred [x[[s]][[t]]$loess_pred < 0 ] <- 0
        x[[s]][[t]]$excess_ratio <- x[[s]][[t]]$loess_pred / median(x[[s]][[t]]$loess_pred[x[[s]][[t]]$loess_pred>0])
      }
    }
  }

  # Calculate position-specific bias coefficients
  for (t in tr_names_i){
    bias_coefficients_list[[t]] <- data.frame(codon_number = x[[1]][[t]]$codon_number, codon_type = x[[1]][[t]]$codon_type, aa_type = x[[1]][[t]]$aa_type)
    for (s in sample_names_i){
      bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]], x[[s]][[t]]$excess_ratio)
    }
    names(bias_coefficients_list[[t]]) <- c("codon_number", "codon_type", "aa_type", sample_names_i)
    bias_coefficient <- apply(bias_coefficients_list[[t]][,-c(1:3)], 1, function(y) gm_mean(y))
    bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]][,c(1:3)], bias_coefficient)
  }
  return(bias_coefficients_list)
}




#' @title CELP_bias
#' @description Function to compute codon-level stalling bias coefficients and bias-corrected read counts using the CELP (Consistent Excess of Loess Preds) method
#' @param tr_codon_read_count_list A list generated by \code{\link{psite_to_codon_count}}
#' containing observed read counts per codon per trancript per sample.
#' @param codon_radius Number of codons on either side of each codon influencing the loess prediction for the middle codon. Default: 5.
#' @param loess_method Determines whether the fitted surface should be computed exactly ("direct") or via interpolation from a kd tree ("interpolate"). Default: "interpolate".
#' @param gini_moderation Logical argument. If set to TRUE, (bias_coefficient)^(gini_index) is used as correction factor.
#' If set to FALSE, bias_coefficient is used as correction factor. Default: FALSE.
#' @details This function is the heart of CELP method for stalling bias detection and correction.
#' It starts with running a loess curve on per-codon read counts along the transcript to borrow information
#' from neighboring codons mitigating the uncertainty of p-site offset assignment and experimental stochasticity.
#' Loess span parameter is calculated from the user-defined codon radius and CDS length.
#' Then, a bias coefficient is calculated for each codon by integrating information on the excess of loess-predicted
#' read counts at that codon comapred to the transcript's background across all samples. Finally, loess predicted counts
#' read counts are divided by bias coefficients to calculate bias-corrected counts.
#' The "direct" fitting method takes longer to complete but does not run into kd-tree-related memory issues.
#' Gini index for each transcript is calculated from the bias coefficients of all of its codons.
#' Gini moderation ensures that the strength of bias correction is proportional to the original level of heterogenity in read distribution along the transcript.
#' @return A list composed of two lists: 1. bias coefficients 2. bias-corrected read counts
#' The bias coefficient list has the following structure: list$<transcript.ID> data.frame: [1] codon_number [2] codon_type [3] aa_type [4] bias_coefficient.
#' The bias-corrected read count list has the following structure: list$<sample.name>$<transcript.ID> data.frame:
#' [1] codon_number [2] codon_type [3] aa_type [4] observed_count [5] bias_coefficient [6] corrected_count.
#' Gini moderation ensures that the strength of correction is proportional to the original level of heterogenity in read distribution along the transcript.
#' @examples
#' tr_codon_bias_coeff_corrected_count_LMCN <- CELP_bias(tr_codon_read_count_LMCN)
#' tr_codon_bias_coeff_corrected_count_LMCN_gini_moderated <- CELP_bias(tr_codon_read_count_LMCN, gini_moderation = TRUE)
#' tr_codon_bias_coeff_corrected_count_LMCN_direct_fit <- CELP_bias(tr_codon_read_count_LMCN, loess_method = "direct")
#' @export
CELP_bias <- function(tr_codon_read_count_list, codon_raduis = 5, loess_method = "interpolate", gini_moderation = FALSE){

  tr_codon_read_count_loess_corrected <- tr_codon_read_count_list
  bias_coefficients_list <- list()
  sample_names_i <- names(tr_codon_read_count_loess_corrected)
  tr_names_i <- names(tr_codon_read_count_loess_corrected[[1]])

  # Run loess and compute loess predicted values
  for (t in tr_names_i){
    l_cds <- dim(tr_codon_read_count_loess_corrected[[1]][[t]])[1]
    span_tr <- (2*codon_raduis+1)/l_cds
    for (s in sample_names_i){
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred <-
        suppressWarnings(predict(loess(tr_codon_read_count_loess_corrected[[s]][[t]]$observed_count ~ tr_codon_read_count_loess_corrected[[s]][[t]]$codon_number, span = span_tr, se = FALSE, control = loess.control(surface = loess_method))))
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred [tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred < 0 ] <- 0
      tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred_by_nz_median <-
        tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred / median(tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred[tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred>0])
    }

    # Calculate position-specific bias coefficients
    bias_coefficients_list[[t]] <- data.frame(codon_number = tr_codon_read_count_loess_corrected[[1]][[t]]$codon_number,
                                              codon_type = tr_codon_read_count_loess_corrected[[1]][[t]]$codon_type,
                                              aa_type = tr_codon_read_count_loess_corrected[[1]][[t]]$aa_type)
    for (s in sample_names_i){
      bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]], tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred_by_nz_median)
    }
    names(bias_coefficients_list[[t]]) <- c("codon_number", "codon_type", "aa_type", sample_names_i)
    bias_coefficient <- apply(bias_coefficients_list[[t]][,-c(1:3)], 1, function(y) gm_mean(y))
    bias_coefficient_gini <- DescTools::Gini(bias_coefficient)
    bias_coefficients_list[[t]] <- data.frame(bias_coefficients_list[[t]][,c(1:3)], bias_coefficient)

    # calculate bias-corrected read counts
    for (s in sample_names_i){
      tr_codon_read_count_loess_corrected[[s]][[t]]$bias_coefficient <- bias_coefficients_list[[t]]$bias_coefficient
      if (gini_moderation == TRUE){
        correction_power <- bias_coefficient_gini
      } else{
        correction_power <- 1
      }
      tr_codon_read_count_loess_corrected[[s]][[t]]$corrected_count <-
        tr_codon_read_count_loess_corrected[[s]][[t]]$loess_pred / (tr_codon_read_count_loess_corrected[[s]][[t]]$bias_coefficient)^correction_power
      tr_codon_read_count_loess_corrected[[s]][[t]] <- subset(tr_codon_read_count_loess_corrected[[s]][[t]], select = -c(loess_pred, loess_pred_by_nz_median))
    }
  }

  output <- list(bias_coefficients_list = bias_coefficients_list, tr_codon_read_count_loess_corrected = tr_codon_read_count_loess_corrected)
  return(output)
}



#' @title codon2transcript
#' @description Function to sum up codon counts per transcript.
#' @param tr_codon_read_count_loess_corrected_list A list of codon level read counts for all samples and transcripts.
#' It is the second element of a tr_codon_bias_coeff_corrected_count object produced by \code{\link{CELP_bias}}.
#' It has the following structure:
#' list$<sample.name>$<transcript.ID> data.frame:
#' [1] codon_number [2] codon_type [3] aa_type [4] observed_count [5] bias_coefficient [6] corrected_count.
#' @param count_type Options: "observed_count", "corrected_count".
#' @details Stalling bias correction is performed at the codon level but differential translational efficiency analysis
#' is ususally performed at transcript level. This function sums up codon level counts (observed or corrected) for each transcript.
#' @examples
#' rpf_observed_sum_LMCN <- codon2transcript(tr_codon_bias_coeff_loess_corrected_count_LMCN$tr_codon_read_count_loess_corrected, "observed_count")
#' rpf_corrected_sum_LMCN <- codon2transcript(tr_codon_bias_coeff_loess_corrected_count_LMCN$tr_codon_read_count_loess_corrected, "corrected_count")
#' @return A data frame where the first column is transcript IDs and the remaining columns contain per transcript read counts for all samples.
#' @export
codon2transcript <- function(tr_codon_read_count_loess_corrected_list, count_type){
  count_sum <- as.data.frame(sapply(tr_codon_read_count_loess_corrected_list, function(x) sapply(x, function(y) sum(y[count_type]))))

  count_sum$transcript <- rownames(count_sum)
  w <- dim(count_sum)[2]
  count_sum <- count_sum[,c(w,1:(w-1))]
  count_sum <- count_sum[order(count_sum$transcript),]
  rownames(count_sum) <- NULL
  return(count_sum)
}



#' @title plot_mirrors_CELP
#' @description Functions to visualize bias coefficient for the specified transcript in one sample.
#' @param x Data frame containing CELP output for a trancript in one sample
#' @param ylim_low_i Lower bound on y axis
#' @param ylim_up_i Upper bound on y axis
#' @param xlim_low_i Lower bound on x axis
#' @param xlim_up_i Upper bound on x axis
#' @param sample Name of the sample to be plotted
plot_mirrors_CELP <- function(x, ylim_low_i, ylim_up_i, xlim_low_i = NULL, xlim_up_i = NULL, sample){
  plot(x$codon_number, y = x$observed_count, type="h",
       xlab = "Codon number", ylab = "Read counts", main = sample,
       ylim = c(ylim_low_i, ylim_up_i), xlim = c(xlim_low_i, xlim_up_i))
  lines(x$bias_coefficient, col = "red")
  lines((-1)*(x$corrected_count), col = "darkorchid1", type = "h")
}



#' @title visualize_CELP
#' @description Function to visualize translational stalling (CELP bias) overlaid on observed and corrected read counts.
#' @param tr_codon_read_count_loess_corrected_list A list of codon level read counts for all samples and transcripts.
#' It is the second element of a tr_codon_bias_coeff_corrected_count object produced by \code{\link{CELP_bias}}.
#' It has the following structure:
#' list$<sample.name>$<transcript.ID> data.frame:
#' [1] codon_number [2] codon_type [3] aa_type [4] observed_count [5] bias_coefficient [8] corrected_count.
#' @param transcript Name of the transcript to be plotted
#' @param panel_rows Number of rows in the case of plotting multiple samples on one page. Default: 1.
#' @param panel_cols Number of columns in the case of plotting multiple samples on one page. Default: 1.
#' @param from_codon Plot starts from this codon. Use this and to_codon to zoom in on particular regions of the transcript. Default: NULL (plot starts at the start codon).
#' @param to_codon Plot ends with this codon. Use this and from_codon to zoom in on particular regions of the transcript. Default: NULL (plot stops at the stop codon).
#' @param outfile Path and name of the output pdf file. If it is not provided, plots will be printed to standard output. Default: NULL.
#' @details This function plots the CELP bias coefficient as a curve overlaid on barplots of observed read counts upward
#' (positive y) and corrected read counts downward (in the nominal negative y range). This allows visual inspection of the
#' prominent bias positions and a comparison of read count heterogeneity along the transcript before and after
#' CELP bias correction.
#' @examples
#' visualize_CELP(tr_codon_bias_coeff_loess_corrected_count_LMCN$tr_codon_read_count_loess_corrected, "ENST00000000233", "<file.path>/ENST00000000233.CELP.bias.plots.pdf")
#' visualize_CELP(tr_codon_bias_coeff_loess_corrected_count_LMCN$tr_codon_read_count_loess_corrected, "ENST00000000233", panel_rows = 2, panel_cols = 4, "<file.path>/ENST00000000233.CELP.bias.all.in.one.page.plots.pdf")
#' @export
visualize_CELP <- function(tr_codon_read_count_loess_corrected_list, transcript, panel_rows = 1, panel_cols = 1, from_codon = NULL, to_codon = NULL, outfile=NULL){
  x_tr <- lapply(tr_codon_read_count_loess_corrected_list, function(y) y[[transcript]])
  ylim_up <- max(unlist(lapply(x_tr, function(y) max(y$observed_count))))
  ylim_low <- (-1) * max(unlist(lapply(x_tr, function(y) max(y$corrected_count))))
  if (is.null(outfile)){
    par(mfrow = c(panel_rows, panel_cols))
    for (s in names(x_tr)){
      plot_mirrors_CELP(x_tr[[s]], ylim_low_i = ylim_low, ylim_up_i = ylim_up, xlim_low_i = from_codon, xlim_up_i = to_codon, s)
    }
  } else {
    pdf(outfile, height = 10, width = 20)
    par(mfrow = c(panel_rows, panel_cols))
    for (s in names(x_tr)){
      plot_mirrors_CELP(x_tr[[s]], ylim_low_i = ylim_low, ylim_up_i = ylim_up, xlim_low_i = from_codon, xlim_up_i = to_codon, s)
    }
    dev.off()
  }
  return()
}
