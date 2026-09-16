#' @title Bundled transcript annotation tables
#' @description Transcript annotation tables listing transcript names and lengths of their 5'UTR, CDS and
#' 3'UTR segments, one per model organism, bundled with the package so \code{\link{load_annotation_and_cdna}}
#' can retrieve them without a separate download (or by loading \code{data("<organism>_annotation")}
#' directly). Produced from a Biomart-derived cDNA reference the same way \code{\link{read_annotation}}
#' parses a user-supplied annotation file; see \code{\link{read_annotation}} for column details.
#' @format A data frame/\code{data.table} with columns \code{transcript}, \code{l_tr}, \code{l_utr5},
#' \code{l_cds}, \code{l_utr3} (one row per transcript).
#' @docType data
#' @keywords datasets
#' @name annotation-datasets
NULL

#' @rdname annotation-datasets
"arabidopsis_annotation"

#' @rdname annotation-datasets
"fly_annotation"

#' @rdname annotation-datasets
"human_annotation"

#' @rdname annotation-datasets
"maize_annotation"

#' @rdname annotation-datasets
"mouse_annotation"

#' @rdname annotation-datasets
"rat_annotation"

#' @rdname annotation-datasets
"worm_annotation"

#' @rdname annotation-datasets
"yeast_annotation"

#' @rdname annotation-datasets
"zebrafish_annotation"


#' @title Bundled cDNA sequences
#' @description cDNA reference sequences, one per model organism, bundled with the package so
#' \code{\link{load_annotation_and_cdna}} can retrieve them without a separate download. Only one
#' transcript per gene (the one with the longest CDS) is included, matching the corresponding
#' \code{<organism>_annotation} table.
#' @format A named list, one element per transcript (named by transcript ID, matching
#' \code{<organism>_annotation$transcript}), each a character vector of individual nucleotides
#' (e.g. \code{c("A", "T", "G", ...)}) making up that transcript's cDNA sequence.
#' @docType data
#' @keywords datasets
#' @name cdna-datasets
NULL

#' @rdname cdna-datasets
"arabidopsis_cdna"

#' @rdname cdna-datasets
"fly_cdna"

#' @rdname cdna-datasets
"human_cdna"

#' @rdname cdna-datasets
"maize_cdna"

#' @rdname cdna-datasets
"mouse_cdna"

#' @rdname cdna-datasets
"rat_cdna"

#' @rdname cdna-datasets
"worm_cdna"

#' @rdname cdna-datasets
"yeast_cdna"

#' @rdname cdna-datasets
"zebrafish_cdna"


#' @title Bundled transcript-to-gene ID mappers
#' @description Mapping between transcript IDs, gene IDs, and gene names, one per model organism,
#' bundled for use with \code{\link{logit_seq}}'s \code{id_mapper}/\code{id_column} arguments (and any
#' other function that accepts an \code{id_mapper}, e.g. \code{\link{visualize_orf_usage}}).
#' @format A data frame with columns:
#' \describe{
#'   \item{transcript}{Transcript ID (e.g. "ENST00000000233" for human).}
#'   \item{gene_id}{Gene ID (e.g. "ENSG00000004059" for human).}
#'   \item{gene_name}{Gene symbol (e.g. "ARF5" for human).}
#' }
#' @docType data
#' @keywords datasets
#' @name id_mapper-datasets
NULL

#' @rdname id_mapper-datasets
"arabidopsis_id_mapper"

#' @rdname id_mapper-datasets
"fly_id_mapper"

#' @rdname id_mapper-datasets
"human_id_mapper"

#' @rdname id_mapper-datasets
"maize_id_mapper"

#' @rdname id_mapper-datasets
"mouse_id_mapper"

#' @rdname id_mapper-datasets
"rat_id_mapper"

#' @rdname id_mapper-datasets
"worm_id_mapper"

#' @rdname id_mapper-datasets
"yeast_id_mapper"

#' @rdname id_mapper-datasets
"zebrafish_id_mapper"


#' @title sample_attributes_dummy
#' @description Example sample design/attributes table for the bundled LMCN demonstration dataset
#' (breast cancer cell lines CN34, LM1a, LM2, and MDA, each with RNA-seq and ribosome profiling
#' replicates), used in package documentation examples as the \code{design} argument to functions such
#' as \code{\link{logit_seq}} and \code{\link{build_empirical_null}}. Named with a \code{_dummy} suffix
#' (rather than the more natural \code{sample_attributes}) specifically to avoid colliding with a user's
#' own same-named variable: since package data is lazy-loaded, a generically-named bundled object would
#' silently shadow/overwrite a user's own \code{sample_attributes} in their workspace.
#' @format A data frame with 16 rows (one per sample) and columns:
#' \describe{
#'   \item{sample_name}{Sample identifier (e.g. "CN34_r1_rna"), matching count matrix column names.}
#'   \item{read_type}{"RNA" or "RPF".}
#'   \item{lung_metastasis}{"Y" or "N": whether the cell line is a lung-metastatic derivative.}
#'   \item{cell_line}{Cell line name: "CN34", "LM1a", "LM2", or "MDA".}
#'   \item{replicate_no}{Replicate number (1 or 2).}
#'   \item{replicate_name}{Replicate identifier shared between the RNA and RPF sample of the same
#'   biological preparation (e.g. "CN34_r1"), used as \code{uniqueID} in \code{\link{partition_to_uniques}}
#'   and related functions.}
#'   \item{cell_line_origin}{Parental cell line group: "CN34" (CN34/LM1a) or "MDA" (MDA/LM2).}
#' }
#' @docType data
#' @keywords datasets
"sample_attributes_dummy"


#' @title rr.v2_dummy
#' @description Combined, normalized RNA-seq and ribosome profiling (RPF) transcript read counts for the
#' bundled LMCN demonstration dataset (cell lines CN34, LM1a, LM2, MDA; 2 replicates each; matching
#' \code{\link{sample_attributes_dummy}}), used in package documentation examples as the \code{x} argument
#' to functions such as \code{\link{logit_seq}} and \code{\link{build_empirical_null}} (typically as
#' \code{rr.v2_dummy[, -1]}, excluding the \code{transcript} ID column). Named with a \code{_dummy} suffix
#' (rather than the more natural \code{rr.v2}) specifically to avoid colliding with a user's own
#' same-named variable: since package data is lazy-loaded, a generically-named bundled object would
#' silently shadow/overwrite a user's own \code{rr.v2} in their workspace.
#' @format A data frame with 11661 rows (one per transcript) and 17 columns: \code{transcript}, followed
#' by 8 RNA-seq columns and 8 RPF columns, one per sample in \code{\link{sample_attributes_dummy}} (column
#' names match \code{sample_attributes_dummy$sample_name}).
#' @docType data
#' @keywords datasets
"rr.v2_dummy"


#' @title reads_list_dummy
#' @description A small (subsampled) RPF reads_list object - the same shape produced by
#' \code{\link{bamtolist_rW}} - for the bundled LMCN demonstration dataset's 8 RPF samples, used in
#' package documentation examples as the starting point for the p-site pipeline
#' (\code{\link{psite_rW}}, \code{\link{psite_info_rW}}, and downstream functions). Subsampled to 100
#' reads per sample so examples run quickly. Named with a \code{_dummy} suffix (rather than the more
#' natural \code{reads_list}) specifically to avoid colliding with a user's own same-named variable:
#' since package data is lazy-loaded, a generically-named bundled object would silently shadow/overwrite
#' a user's own \code{reads_list} in their workspace.
#' @format A named list of 8 data frames/\code{data.table}s, one per RPF sample (e.g. "CN34_r1_rpf"),
#' each with 100 rows (one per read) and columns \code{transcript}, \code{end5}, \code{end3},
#' \code{length}, \code{cds_start}, \code{cds_stop} - matching \code{\link{bamtolist_rW}}'s output - plus
#' the p-site columns (\code{psite}, \code{psite_from_start}, \code{psite_from_stop},
#' \code{psite_region}) that \code{\link{psite_info_rW}} would add. These extra columns are harmless
#' leftovers: \code{\link{psite_rW}} ignores them and \code{\link{psite_info_rW}} recomputes and
#' overwrites them, so \code{reads_list_dummy} can be fed into the p-site pipeline exactly like a raw
#' \code{\link{bamtolist_rW}} output.
#' @docType data
#' @keywords datasets
"reads_list_dummy"
