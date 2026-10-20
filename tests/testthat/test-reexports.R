# The annotation step (prepare_annotation_files(), load_annotation() and the class
# FaFile_Circ) is in RiboseQC, and its tests are there. ORFquant exports it again.
# The same for get_ps_fromspliceplus(), get_ps_fromsplicemin() and prepare_for_ORFquant().

test_that("ORFquant exports the annotation step of RiboseQC", {
  expect_identical(ORFquant::prepare_annotation_files, RiboseQC::prepare_annotation_files)
  expect_identical(ORFquant::load_annotation, RiboseQC::load_annotation)
  expect_identical(ORFquant::FaFile_Circ, RiboseQC::FaFile_Circ)
  expect_true(all(c("prepare_annotation_files", "load_annotation", "FaFile_Circ") %in%
                    getNamespaceExports("ORFquant")))
})

test_that("ORFquant exports the P-site functions of RiboseQC for spliced reads", {
  expect_identical(ORFquant::get_ps_fromspliceplus, RiboseQC::get_ps_fromspliceplus)
  expect_identical(ORFquant::get_ps_fromsplicemin, RiboseQC::get_ps_fromsplicemin)
  expect_true(all(c("get_ps_fromspliceplus", "get_ps_fromsplicemin") %in%
                    getNamespaceExports("ORFquant")))
})

test_that("ORFquant exports prepare_for_ORFquant() of RiboseQC", {
  expect_identical(ORFquant::prepare_for_ORFquant, RiboseQC::prepare_for_ORFquant)
  expect_true("prepare_for_ORFquant" %in% getNamespaceExports("ORFquant"))
})

# run_ORFquant() needs this method (R/helpers.R). It is in RiboseQC.
test_that("unlist() of a GRanges gives the GRanges back", {
  gr <- GenomicRanges::GRanges(c("chr1:1-10:+", "chr1:20-30:-"))
  expect_identical(unlist(gr), gr)
  expect_identical(unlist(GenomicRanges::GRanges()), GenomicRanges::GRanges())
})

# The tests run in the namespace of ORFquant, so they find getSeq() and seqinfo() also
# when ORFquant doesn't export them. The vignette needs the exports.
test_that("ORFquant exports getSeq() and seqinfo()", {
  expect_identical(ORFquant::getSeq, Biostrings::getSeq)
  expect_identical(ORFquant::seqinfo, Seqinfo::seqinfo)
})

test_that("a FaFile_Circ gives its sequence and its circular chromosomes", {
  genome <- FaFile_Circ(Rsamtools::FaFile(example_file("chr22_example.fa.gz")), circularRanges = "chr22")
  expect_true(isCircular(seqinfo(genome))[["chr22"]])
  region <- GenomicRanges::GRanges("chr22", IRanges::IRanges(1000001, 1000010))
  expect_identical(as.character(unlist(getSeq(genome, region))),
                   as.character(unlist(getSeq(Rsamtools::FaFile(example_file("chr22_example.fa.gz")), region))))
})

# data/ORFquant_1.03.0_chr22_example.Rannot is the annotation of the example that
# prepare_annotation_files() of ORFquant 1.03.0 made, with its file paths reduced to the
# file names. Its genome holds a FaFile_Circ of the class of ORFquant, not of RiboseQC.
test_that("an annotation that ORFquant made before RiboseQC took over loads, with a genome that works", {
  where <- new.env()
  loaded <- NULL
  suppressMessages(loaded <- load_annotation(test_path("data", "ORFquant_1.03.0_chr22_example.Rannot"),
                                             envir = where))
  expect_named(loaded, c("GTF_annotation", "genome_seq"))
  expect_identical(get("GTF_annotation", envir = where), loaded$GTF_annotation)
  expect_identical(get("genome_seq", envir = where), loaded$genome_seq)
  genome <- loaded$genome_seq
  expect_identical(attr(class(genome), "package"), "ORFquant")
  expect_true(is(genome, "FaFile_Circ"))
  # the file has the name of the genome, not its path: point it to the example
  fasta <- example_file("chr22_example.fa.gz")
  genome$field("path", fasta)
  genome$field("index", paste0(fasta, ".fai"))
  genome$field("gzindex", paste0(fasta, ".gzi"))
  expect_equal(unname(seqlengths(seqinfo(genome))), 50818468)
  cds <- unlist(loaded$GTF_annotation$cds_txs)[1:2]
  expect_identical(as.character(getSeq(genome, cds)),
                   as.character(getSeq(Rsamtools::FaFile(fasta), cds)))
})
