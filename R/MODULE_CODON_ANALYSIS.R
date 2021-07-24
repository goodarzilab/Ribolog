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

codon_usage_psite_rW <- function(data, annotation, sample, site = "psite",
                              fastapath = NULL, fasta_genome = TRUE,
                              refseq_sep = NULL, bsgenome = NULL, gtfpath = NULL,
                              txdb = NULL, dataSource = NA, organism = NA,
                              transcripts = NULL, frequency_normalization = TRUE,
                              codon_values = NULL, label_scatter = FALSE,
                              label_number = 64, label_aminoacid = FALSE) {

  if(site != "psite" & site != "asite" & site != "esite"){
    cat("\n")
    stop("parameter site must be either \"psite\", \"asite\" or \"esite\" \n\n")
  }

  if(length(sample) > 2){
    cat("\n")
    stop("max two samples can be specified \n\n")
  }

  if(length(codon_values) != 0){
    rownames(codon_values) <- as.character(codon_values$codon)
    if (nrow(codon_values) < 64 | length(unique(rownames(codon_values))) < 64){
      cat("\n")
      stop("number of different triplets in codon_values < 64. At least one codon is missing\n\n")
    } else {
      if (length(unique(rownames(codon_values))) > 64){
        cat("\n")
        stop("number of different triplets in codon_values > 64. Too many codons\n\n")
      }
    }
  }

  if(length(sample) == 2 & length(codon_values) != 0){
    cat("\n")
    warning("codon_values will be ingnored: two samples are specified\n")
    codon_values = NULL
  }

  if((sum(lapply(data[sample], function(x) (gsub("site", "_site_codon", site) %in% colnames(x))) == TRUE) < length(sample)) |
     ((sum(lapply(data[sample], function(x) (gsub("site", "_site_codon", site) %in% colnames(x))) == TRUE) == length(sample)) &
      (frequency_normalization == TRUE | frequency_normalization == T))){

    if(length(fastapath) == 0 & length(bsgenome) == 0){
      cat("\n")
      if(sum(lapply(data[sample], function(x) (gsub("site", "_site_codon", site) %in% colnames(x))) == TRUE) < length(sample)){
        stop(paste0("unable to add ", gsub("site", "_site_codon", site), " column. Either fastpath or bsgenome must be specified\n\n"))
      } else {
        warning("no nucleotide sequences provided. Either fastpath or bsgenome must be specified. Non-normalized codon usage indexes will be returned\n")
        frequency_normalization = FALSE
      }
    }

    if(((length(fastapath) != 0 & (fasta_genome == TRUE | fasta_genome == T)) |
        length(bsgenome) != 0) &
       length(gtfpath) == 0 & length(txdb) == 0){
      cat("\n")
      stop("genome annotation file not specified. Either gtfpath or txdb object must be specified\n\n")
    }

    if(length(fastapath) != 0 & length(bsgenome) != 0){
      cat("\n")
      warning("both fastapath and bsgenome are specified. Only fastapath will be considered\n")
      bsgenome = NULL
    }

    if(length(gtfpath) != 0 & length(txdb) != 0){
      cat("\n")
      warning("both gtfpath and txdb are specified. Only gtfpath will be considered\n")
      txdb = NULL
    }

    if((length(gtfpath) != 0 | length(txdb) != 0) &
       ((length(fastapath) == 0 & length(bsgenome) == 0) |
        (length(fastapath) != 0 & (fasta_genome == FALSE | fasta_genome == F)))){
      cat("\n")
      stop("genome annotation file specified but no sequences from genome assembly provided\n\n")
    }

    if(length(gtfpath) != 0 | length(txdb) != 0){
      if(length(gtfpath) != 0){
        path_to_gtf <- gtfpath
        suppressWarnings(txdbanno <- GenomicFeatures::makeTxDbFromGFF(file = path_to_gtf,
                                                                      format = "gtf",
                                                                      dataSource = dataSource,
                                                                      organism = organism))
      } else {
        if(txdb %in% rownames(installed.packages())){
          library(txdb, character.only = TRUE)
        } else {
          source("https://bioconductor.org/biocLite.R")
          biocLite(txdb, suppressUpdates = TRUE)
          library(txdb, character.only = TRUE)
        }
        txdbanno <- get(txdb)
      }
    }

    if(length(fastapath) != 0 | length(bsgenome) != 0){
      if(length(fastapath) != 0) {
        if(fasta_genome == TRUE | fasta_genome == T){
          temp_sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(temp_sequences) <- tstrsplit(names(temp_sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
          exon <- suppressWarnings(GenomicFeatures::exonsBy(txdbanno, by = "tx", use.names = TRUE))
          exon <- data.table::as.data.table(exon[unique(names(exon))])
          sub_exon_plus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "+"]
          sub_exon_minus <- exon[as.character(seqnames) %in% names(temp_sequences) & strand == "-"
                                 ][, new_end := Biostrings::width(temp_sequences[as.character(seqnames)]) - start + 1
                                   ][, new_start := Biostrings::width(temp_sequences[as.character(seqnames)]) - end + 1]

          seq_dt_plus <- sub_exon_plus[, nt_seq := "emp"
                                       ][, nt_seq := as.character(Biostrings::subseq(temp_sequences[as.character(seqnames)],
                                                                                     start = start,
                                                                                     end = end))
                                         ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]

          revcompl_temp_sequences <- Biostrings::reverseComplement(temp_sequences)
          seq_dt_minus <- sub_exon_minus[, nt_seq := "emp"
                                         ][, nt_seq := as.character(Biostrings::subseq(revcompl_temp_sequences[as.character(seqnames)],
                                                                                       start = new_start,
                                                                                       end = new_end))
                                           ][, list(seq = paste(nt_seq, collapse = "")), by = group_name]

          sequences <- Biostrings::DNAStringSet(c(seq_dt_plus$seq, seq_dt_minus$seq))
          names(sequences) <- c(unique(sub_exon_plus$group_name), unique(sub_exon_minus$group_name))
        } else {
          sequences <- Biostrings::readDNAStringSet(fastapath, format = "fasta", use.names = TRUE)
          if(length(refseq_sep) != 0){
            names(sequences) <- tstrsplit(names(sequences), refseq_sep, fixed = TRUE, keep = 1)[[1]]
          }
        }
      } else {
        if(bsgenome %in% installed.genomes()){
          library(bsgenome, character.only = TRUE)
        } else {
          source("http://www.bioconductor.org/biocLite.R")
          biocLite(bsgenome, suppressUpdates = TRUE)
          library(bsgenome, character.only = TRUE)
        }
        sequences <- GenomicFeatures::extractTranscriptSeqs(get(bsgenome), txdbanno, use.names = T)
      }
    }

    for(samp in sample){
      if(!(gsub("site", "_site_codon", site) %in% colnames(data[[samp]]))){
        if(site == "psite" & !("p_site_codon" %in% colnames(data[[samp]]))){
          data[[samp]] <- data[[samp]][, p_site_codon := as.character(Biostrings::subseq(sequences[as.character(data[[samp]]$transcript)],
                                                                                         start = data[[samp]]$psite,
                                                                                         end = data[[samp]]$psite + 2))]
        }
        if(site == "asite" & !("a_site_codon" %in% colnames(data[[samp]]))){
          data[[samp]] <- data[[samp]][, a_site_codon := as.character(Biostrings::subseq(sequences[as.character(data[[samp]]$transcript)],
                                                                                         start = data[[samp]]$psite + 3,
                                                                                         end = data[[samp]]$psite + 5))]
        }
        if(site == "esite" & !("e_site_codon" %in% colnames(data[[samp]]))){
          data[[samp]] <- data[[samp]][, e_site_codon := as.character(Biostrings::subseq(sequences[as.character(data[[samp]]$transcript)],
                                                                                         start = data[[samp]]$psite - 3,
                                                                                         end = data[[samp]]$psite - 1))]
        }
      }
    }
  }

  l_transcripts <- as.character(annotation[l_cds > 0 & l_cds %% 3 == 0, transcript])

  cod_aa <- data.table(codon=c("GCC", "GCG", "GCU", "GCA", "AGA", "CGG", "AGG", "CGA", "CGC", "CGU", "AAC", "AAU", "GAC", "GAU", "UGC", "UGU", "CAA", "CAG", "GAG", "GAA", "GGC", "GGU", "GGA", "GGG", "CAC", "CAU", "AUA", "AUC", "AUU", "CUG", "CUA", "UUA", "CUU", "UUG", "CUC", "AAA", "AAG", "AUG", "UUC", "UUU", "CCG", "CCC", "CCU", "CCA", "AGC", "UCG", "UCU", "UCA", "UCC", "AGU", "UAG", "UAA", "UGA", "ACA", "ACC", "ACG", "ACU", "UGG", "UAU", "UAC", "GUA", "GUG", "GUU", "GUC"),
                       aa=c("A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D", "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S", "S", "S", "*", "*", "*", "T", "T", "T", "T", "W", "Y", "Y", "V", "V", "V", "V"))
  cod_lev <- gsub("U", "T", cod_aa$codon)

  if (length(transcripts) == 0) {
    c_transcript <- l_transcripts
  } else {
    c_transcript <- intersect(l_transcripts, transcripts)
  }

  if(frequency_normalization == TRUE | frequency_normalization == T){
    if(setequal(names(sequences), as.character(annotation$transcript)) == FALSE){
      exc_seq <- length(setdiff(names(sequences), as.character(annotation$transcript)))
      exc_anno <- length(setdiff(as.character(annotation$transcript), names(sequences)))

      if(exc_seq > 0){
        cat("\n")
        warning(sprintf("more reference sequences (from FASTA or BSgenome) than transcript IDs (from annotation data table)\n  %s reference sequences discarded", exc_seq))
      }

      if(exc_anno > 0){
        cat("\n")
        warning(sprintf("more transcript IDs (from annotation data table) than sequences (from FASTA or BSgenome)\n  %s reference sequences discarded", exc_seq))
      }
    }

    c_transcript <- intersect(c_transcript, names(sequences))
    sub_sequences <- sequences[c_transcript]
    cds_biost <- Biostrings::subseq(sub_sequences,
                                    start = annotation[transcript %in% names(sub_sequences), l_utr5] + 1,
                                    end = annotation[transcript %in% names(sub_sequences), l_utr5] +
                                      annotation[transcript %in% names(sub_sequences), l_cds])

    seq_freq <- (Biostrings::trinucleotideFrequency(cds_biost, step = 3,
                                                    as.prob = TRUE,
                                                    with.labels = TRUE,
                                                    simplify.as = "collapsed")) * 1000
    names(seq_freq) <- gsub("T", "U", names(seq_freq))
  }

  if(exists("norm_table")){
    rm(norm_table)
  }
  for(samp in sample) {
    dt <- data[[samp]][as.character(transcript) %in% c_transcript & psite_from_start %% 3 == 0]

    if(site == "psite"){
      dt <- dt[psite_region == "cds"]
      temp_table <- data.table(table(factor(dt$p_site_codon, levels = cod_lev)))
    } else {
      if(site == "asite"){
        dt <- dt[psite_from_start >= -3 & psite_from_stop <= -3]
        temp_table <- data.table(table(factor(dt$a_site_codon, levels = cod_lev)))
      } else {
        dt <- dt[psite_from_start >= 3 & psite_from_stop <= 3]
        temp_table <- data.table(table(factor(dt$e_site_codon, levels = cod_lev)))
      }
    }
    colnames(temp_table) <- c("codon", "raw_value")

    temp_table <- temp_table[, codon := gsub("T", "U", codon)]

    if(frequency_normalization == TRUE | frequency_normalization == T){
      temp_table <- temp_table[, normalized_value := raw_value / seq_freq[codon]]
      colplot <- "normalized_value"
    } else {
      colplot <- "raw_value"
    }

    temp_table <- temp_table[, plot_value := get(colplot) - min(get(colplot))
                             ][, plot_value := plot_value / max(plot_value)]

    if(!exists("norm_table")){
      norm_table <- temp_table
    } else {
      if(frequency_normalization == TRUE | frequency_normalization == T){
        norm_table[, "comp_normalized_value" := temp_table$normalized_value]
      }
      norm_table[, c("comp_raw_value", "comp_plot_value") := list(temp_table$raw_value, temp_table$plot_value)]
    }
  }

  norm_table <- norm_table[, class := "cds"
                           ][codon == "AUG", class := "start"
                             ][codon %in% c("UAA", "UGA", "UAG"), class := "stop"
                               ][, class := factor(class, levels = c("start", "cds", "stop"),
                                                   labels=c("Start codon", "cds", "Stop codon"))
                                 ][cod_aa, on = "codon"]

  if(length(sample) == 1){
    output <- list()
    if(frequency_normalization == TRUE | frequency_normalization == T){
      output[["dt"]] <- norm_table[, c("codon",  "aa", "class", "raw_value", "normalized_value", "plot_value")]
    } else {
      output[["dt"]] <- norm_table[, c("codon",  "aa", "class", "raw_value", "plot_value")]
    }

    if(length(codon_values) != 0){
      codon_values <- codon_values[, list(codon, value)]
      setnames(codon_values, old = c("codon", "value"), new =  c("codon", "comp_values"))
      codon_values[, codon := gsub("T", "U", codon)]
      norm_table <- norm_table[codon_values, on = "codon"
                               ][, comp_plot_value := comp_values - min(comp_values)
                                 ][, comp_plot_value := comp_plot_value / max(comp_plot_value)]

      output[["dt"]] <- cbind(output[["dt"]], norm_table[, c("comp_values", "comp_plot_value")])
      setnames(output[["dt"]], old = c("comp_values", "comp_plot_value"), new = c("user_value", "plot_user_value"))
      nc <- ncol(output[["dt"]])
      output[["dt"]] <- output[["dt"]][, c(1:(nc - 3), nc -1 , nc - 2, nc), with = FALSE]

      plotnamex <- "Codon usage index"
      plotnamey <- "User value"
    }
  } else {
    output <- list()
    if(frequency_normalization == TRUE | frequency_normalization == T){
      output[["dt"]] <- norm_table[, c("codon", "aa", "class", "raw_value", "comp_raw_value",
                                       "normalized_value", "comp_normalized_value", "plot_value", "comp_plot_value")]
      setnames(output[["dt"]],
               old = c("normalized_value", "comp_normalized_value"),
               new = paste0("normalized_value_", c(sample[1], sample[2])))
    } else {
      output[["dt"]] <- norm_table[, c("codon", "aa", "class", "raw_value", "comp_raw_value",
                                       "plot_value", "comp_plot_value")]
    }

    setnames(output[["dt"]],
             old = c("raw_value", "comp_raw_value", "plot_value", "comp_plot_value"),
             new = paste(rep(c("raw_value", "plot_value"),each=2), c(sample[1], sample[2]), sep = "_"))

    plotnamex <- sample[1]
    plotnamey <- sample[2]
  }

  oldw <- getOption("warn")
  options(warn=-1)

  bs <- 30
  plot_list <- list()
  i = 1
  for(samp in sample){
    plot_col <- ifelse(i == 1, "plot_value", "comp_plot_value")
    i <- i + 1

    norm_table <- norm_table[order(get(plot_col))
                             ][, codon := factor(codon, levels = codon)]

    colour_codon <- ifelse(norm_table$codon == "AUG", "#104ec1",
                           ifelse(norm_table$codon %in% c("UAA", "UGA", "UAG"),
                                  "darkred", "gray40"))

    bp <- ggplot(norm_table, aes_string(x = "codon", y = plot_col, fill = "class")) +
      geom_bar(stat = "identity", alpha = 0.9) +
      scale_fill_manual(name = "", breaks=c("Start codon","Stop codon"), values = c("#104ec1", "gray60", "darkred")) +
      theme_bw(base_size = bs) +
      theme(legend.position = "top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(5,0,-15,0)) +
      theme(legend.text = element_text(margin = margin(l = -12, unit = "pt"))) +
      scale_x_discrete("Codon") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = bs * 0.5)) +
      geom_text(aes(label = aa), vjust = -0.5, size = bs/6) +
      theme(axis.text.x = element_text(colour = colour_codon), panel.grid.major.x = element_blank(),
                                       panel.grid.minor.y = element_blank()) +
      scale_y_continuous("Usage index", limits = c(0,1.075), breaks = c(0,0.25,0.5,0.75,1))
    plot_list[[paste0("plot_", samp)]] <- bp
  }

  if(length(sample) == 1){
    output[["plot"]] <- plot_list[[paste0("plot_", sample)]]
  }

  if(length(codon_values) != 0 | length(sample) == 2){

    norm_table <- norm_table[order(plot_value)
                             ][, codon := factor(codon, levels = codon)]

    correlation <- round(cor(norm_table$plot_value, norm_table$comp_plot_value), 3)

    bs <- 25
    pcomp <- ggplot(norm_table,aes(x = plot_value, y = comp_plot_value, colour = class)) +
      geom_smooth(method = "lm", se=T, color="gray80", fill="gray80", linetype = 1, formula = y ~ x, level = 0.99,  fullrange = TRUE) +
      geom_point(alpha = 0.9, size = bs * 0.14) +
      scale_colour_manual(name = "", breaks = c("Start codon","Stop codon"), values = c("#104ec1", "gray40", "darkred")) +
      theme_bw(base_size = bs) +
      theme(legend.position = "top", legend.margin=margin(0,0,0,0), legend.box.margin=margin(5,0,-15,0)) +
      theme(legend.text = element_text(margin = margin(l = -12, unit = "pt"))) +
      theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) +
      scale_x_continuous(name = plotnamex, limits = c(-0.3,1.3), breaks = c(0,0.25,0.5,0.75,1), expand = c(0,0)) +
      scale_y_continuous(name = plotnamey, limits = c(-0.3,1.3), breaks = c(0,0.25,0.5,0.75,1), expand = c(0,0)) +
      coord_cartesian(xlim = c(-0.05,1.05), ylim = c(-0.05,1.05)) +
      annotate("text", x = 1, y = 0, label = paste0("R=",correlation), vjust = -0.2, size = bs * 0.2, hjust = 1, color = "black")

    if (label_scatter == T || label_scatter == TRUE) {
      if (label_number != 64){
        fit <- lm(norm_table$comp_plot_value ~ norm_table$plot_value)
        outlier_tab <- data.table(predict.lm(fit, interval = "confidence", level = 0.95)
        )[, codon := norm_table$codon
          ][, comp_plot_value := norm_table$comp_plot_value
            ][comp_plot_value < lwr | comp_plot_value > upr]
        outlier_cod <- as.character(outlier_tab[, max_dist := min(abs(comp_plot_value - lwr), abs(comp_plot_value - upr)),
                                                by = 1:nrow(outlier_tab)
                                                ][order(-max_dist)
                                                  ][1:label_number, codon])
        lab_table <- norm_table[, codon_out := ""
                                 ][as.character(codon) %in% outlier_cod, codon_out := codon
                                   ][, aa_out := ""
                                     ][codon_out != "", aa_out := aa]

        if (label_aminoacid == T || label_aminoacid == TRUE) {
          label_col <- "aa_out"
        } else {
          label_col <- "codon_out"
        }
      } else {
        lab_table <- norm_table
        if (label_aminoacid == T || label_aminoacid == TRUE) {
          label_col <- "aa"
        } else {
          label_col <- "codon"
        }
      }

      pcomp <- pcomp +
        ggrepel::geom_text_repel(data = lab_table,
                                 aes_string("plot_value", "comp_plot_value", label = as.character(label_col)), show.legend = F)
    }

    if(length(sample) == 2){
      output[[paste0("plot_", sample[1])]] <- plot_list[[paste0("plot_", sample[1])]]
      output[[paste0("plot_", sample[2])]] <- plot_list[[paste0("plot_", sample[2])]]
    }

    output[["plot_comparison"]] <- pcomp
  }

  options(warn = oldw)
  return(output)
}



#' @title codon_barplot
#' @description Function to visualise the codon usage in every sample
#' @param reads_psite_list Output object of a psite_info_rW function
#' @param annotation A dataframe denoting the annotations
#' @param fasta_path A path to the fasta file
#' @param outfile file name to output visualisations
#' @return Barplots for each sample
#' @export

codon_barplot <- function(reads_psite_list, annotation, fasta_path, outfile=NULL){

  if (!is.null(outfile)) { pdf(outfile)}

  for (sample_i in names(reads_psite_list)){
    output <- codon_usage_psite_rW(reads_psite_list, annotation, sample = c(sample_i),
                                fastapath = fasta_path, fasta_genome = FALSE, frequency_normalization = FALSE)
    print(output[["plot"]]  + ggtitle(sample_i) + theme(plot.title = element_text(hjust = 0.5)))
  }

  if (!is.null(outfile)) {
    dev.off()
    sprintf("PDF (%s) created and saved", outfile)
  }
}



#' @title codon_comparative_scatterplot
#' @description Function to visualise the comparative codon usage for every sample
#' @param reads_psite_list Output object of a psite_info_rW function
#' @param annotation A dataframe denoting the annotations
#' @param fasta_path A path to the fasta file
#' @param outfile file name to output visualisations
#' @return Barplots for each sample
#' @export

codon_comparative_scatterplot <- function(reads_psite_list, annotation, fasta_path, outfile=NULL){

  if (!is.null(outfile)) { pdf(outfile)}

  for (sample_i in names(reads_psite_list)){

     for (sample_2 in names(reads_psite_list)[names(reads_psite_list) != sample]){

    output <- codon_usage_psite_rW(reads_psite_list, annotation, sample = c(sample_i, sample_2),
                                fastapath = fasta_path, fasta_genome = FALSE, frequency_normalization = FALSE)

    options(repr.plot.width = 10, repr.plot.height = 10)
    print(output[["plot_comparison"]])

    }
  }

  if (!is.null(outfile)) {
    dev.off()
    sprintf("PDF (%s) created and saved", outfile)
  }
}
