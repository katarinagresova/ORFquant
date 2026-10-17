# The annotation step is in RiboseQC. ORFquant exports it again, so scripts that call
# prepare_annotation_files(), load_annotation() or FaFile_Circ() after library(ORFquant)
# don't change. The methods of getSeq() and seqinfo() for FaFile_Circ are in RiboseQC.

#' @importFrom RiboseQC prepare_annotation_files
#' @export
RiboseQC::prepare_annotation_files

#' @importFrom RiboseQC load_annotation
#' @export
RiboseQC::load_annotation

#' @importFrom RiboseQC FaFile_Circ
#' @export
RiboseQC::FaFile_Circ

#' @importClassesFrom RiboseQC FaFile_Circ
#' @exportClass FaFile_Circ
#' @exportMethod getSeq
#' @exportMethod seqinfo
#' @name FaFile_Circ-class
#' @rdname reexports
NULL
