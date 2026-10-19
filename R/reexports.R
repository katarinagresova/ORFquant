# The annotation step is in RiboseQC. ORFquant exports it again, so scripts that call
# prepare_annotation_files(), load_annotation() or FaFile_Circ() after library(ORFquant)
# don't change. The methods of getSeq() and seqinfo() for FaFile_Circ are in RiboseQC.
# The same for the P-sites of spliced reads: get_ps_fromspliceplus() and
# get_ps_fromsplicemin(), which prepare_for_ORFquant() uses.

#' @importFrom RiboseQC prepare_annotation_files
#' @export
RiboseQC::prepare_annotation_files

#' @importFrom RiboseQC load_annotation
#' @export
RiboseQC::load_annotation

#' @importFrom RiboseQC FaFile_Circ
#' @export
RiboseQC::FaFile_Circ

#' @importFrom RiboseQC get_ps_fromspliceplus
#' @export
RiboseQC::get_ps_fromspliceplus

#' @importFrom RiboseQC get_ps_fromsplicemin
#' @export
RiboseQC::get_ps_fromsplicemin

# ORFquant has no S4 classes or methods of its own, so R ignores exportClasses() and
# exportMethods() in its NAMESPACE. It exports the generics getSeq() and seqinfo() as
# functions, so they work after library(ORFquant), as before.

#' @importFrom Biostrings getSeq
#' @export
Biostrings::getSeq

#' @importFrom Seqinfo seqinfo
#' @export
Seqinfo::seqinfo

#' @importClassesFrom RiboseQC FaFile_Circ
#' @name FaFile_Circ-class
#' @rdname reexports
NULL
