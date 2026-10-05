# Package imports. NAMESPACE is generated from these tags by roxygen2.
#' @rawNamespace import(Biostrings, except = pattern)
#' @import GenomicAlignments
#' @import rtracklayer
#' @import BSgenome
#' @rawNamespace import(BiocGenerics, except = c(combine, Position))
#' @import multitaper
#' @import reshape2
#' @import ggplot2
#' @import cowplot
#' @import grid
#' @import gridExtra
#' @import GenomicFeatures
#' @importClassesFrom Rsamtools FaFile
#' @importFrom AnnotationDbi saveDb
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges %over% CharacterList IRanges IntegerList LogicalList
#' @importFrom IRanges NumericList RleList countOverlaps disjoin nearest
#' @importFrom IRanges quantile reduce resize restrict shift subsetByOverlaps
#' @importFrom IRanges trim
#' @importFrom Rsamtools BamFile ScanBamParam scanBamFlag scanFa scanFaIndex
#' @importFrom S4Vectors mcols<- DataFrame Rle aggregate cbind.DataFrame
#' @importFrom S4Vectors endoapply na.omit queryHits runLength runValue
#' @importFrom S4Vectors setMethods splitAsList subjectHits
#' @importFrom Seqinfo Seqinfo isCircular seqlengths<- seqlengths seqlevels
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices colorRampPalette dev.off pdf rgb
#' @importFrom methods as is new representation setMethod setRefClass
#' @importFrom stats as.ts optim pf setNames
#' @importFrom utils read.table write.table
NULL

# Column names used in aes() and subset(), which R CMD check takes for globals.
utils::globalVariables(c("Freq", "ORFs_tx.ORF_pct_P_sites", "ORFs_tx.ORFs_pM",
  "Position_from_splice_site", "Var1", "coll", "feature", "gene_id", "group",
  "n_tx", "pval", "pval_uniq", "ref", "transcript", "value", "width.ORFs_tx.",
  "x1", "x2", "y1", "y2"))
# Globals: load_annotation() assigns GTF_annotation and genome_seq with <<-, and
# ORFquant() and others read them; other functions load() GTF_annotation and
# ORFquant_results. globalVariables() does not cover <<-, so R CMD check still
# reports load_annotation()'s to genome_seq and run_ORFquant()'s to for_ORFquant.
utils::globalVariables(c("GTF_annotation", "genome_seq", "ORFquant_results"))
