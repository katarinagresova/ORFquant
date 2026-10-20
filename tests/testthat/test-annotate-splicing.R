# The splice features of annotate_splicing(): the - strand must give the labels
# of the + strand. Before, new and missing exons at the 5' end of a - strand ORF
# got no _5prime, and an ORF that overlaps no reference exon got _3prime, not
# _notoverl (R/annotate_orfs.R L460 @ 6abfe9d).

gr <- function(s, e, st) GenomicRanges::GRanges("chrT", IRanges::IRanges(s, e), strand = st)
# The mirror image of + strand exons on the - strand, in ascending order as in
# run_ORFquant()
mirror <- function(x) sort(gr(1001 - BiocGenerics::end(x), 1001 - BiocGenerics::start(x), "-"))
# spl_type of each exon in 5' to 3' order
spl_5to3 <- function(r) {
  r$spl_type[order(BiocGenerics::start(r), decreasing = all(BiocGenerics::strand(r) == "-"))]
}

splicing_cases <- list(
  # ORF, reference CDS and spl_type in 5' to 3' order, on the + strand
  list(orf = gr(100, 150, "+"), ref = gr(c(300, 500), c(400, 600), "+"),
       spl = c("new_monoCDS_notoverl", "missing_CDS_notoverl", "missing_CDS_notoverl")),
  list(orf = gr(c(100, 300, 500), c(150, 400, 600), "+"), ref = gr(c(300, 500), c(400, 600), "+"),
       spl = c("new_firstCDS_5prime", "same_5ss;same_3ss", "same_5ss;same_lastCDS")),
  list(orf = gr(c(320, 500), c(400, 600), "+"), ref = gr(c(100, 300, 500), c(200, 400, 600), "+"),
       spl = c("missing_CDS_5prime", "down_firstCDS;same_3ss", "same_5ss;same_lastCDS")),
  list(orf = gr(c(100, 300), c(200, 380), "+"), ref = gr(c(100, 300, 500), c(200, 400, 600), "+"),
       spl = c("same_firstCDS;same_3ss", "same_5ss;up_lastCDS", "missing_CDS_3prime")),
  list(orf = gr(c(50, 300, 700), c(80, 400, 750), "+"), ref = gr(c(300, 500), c(400, 600), "+"),
       spl = c("new_firstCDS_5prime", "same_5ss;same_3ss", "missing_CDS_3prime", "new_lastCDS_3prime")),
  list(orf = gr(c(300, 500), c(400, 600), "+"), ref = gr(c(300, 500), c(400, 600), "+"),
       spl = c("same_firstCDS;same_3ss", "same_5ss;same_lastCDS")))

test_that("annotate_splicing() gives the same spl_type on both strands", {
  for (k in seq_along(splicing_cases)) {
    p <- splicing_cases[[k]]
    expect_identical(spl_5to3(annotate_splicing(p$orf, p$ref)), p$spl, info = paste("+ strand, case", k))
    expect_identical(spl_5to3(annotate_splicing(mirror(p$orf), mirror(p$ref))), p$spl, info = paste("- strand, case", k))
  }
})

test_that("annotate_splicing() does not depend on the order of the ORF exons", {
  # Before, exons in descending order got the overlap of another exon: wrong
  # labels, or an error when an exon got none (R/annotate_orfs.R L49 @ 4c1fd44)
  for (k in seq_along(splicing_cases)) {
    p <- splicing_cases[[k]]
    expect_identical(annotate_splicing(rev(p$orf), p$ref), annotate_splicing(p$orf, p$ref), info = paste("+ strand, case", k))
    m <- mirror(p$orf); mref <- mirror(p$ref)
    expect_identical(annotate_splicing(rev(m), mref), annotate_splicing(m, mref), info = paste("- strand, case", k))
  }
})

test_that("cds_without_stop() removes the stop codon from the CDS exons, also when an intron splits it", {
  # The CDS exons with the stop codon, and without it, on the + strand
  cases <- list(
    list(cds = gr(c(100, 300), c(200, 350), "+"), out = gr(c(100, 300), c(200, 347), "+")),
    # 1 nt of the stop codon before the intron, 2 nt after it
    list(cds = gr(c(100, 300), c(200, 301), "+"), out = gr(100, 199, "+")),
    # the stop codon is the last exon
    list(cds = gr(c(100, 300), c(200, 302), "+"), out = gr(100, 200, "+")))
  for (k in seq_along(cases)) {
    p <- cases[[k]]
    expect_identical(cds_without_stop(p$cds, 3), p$out, info = paste("+ strand, case", k))
    expect_identical(cds_without_stop(mirror(p$cds), 3), mirror(p$out), info = paste("- strand, case", k))
    # the exons keep their order
    expect_identical(cds_without_stop(rev(p$cds), 3), rev(p$out), info = paste("+ strand, descending, case", k))
  }
  # A CDS without its stop codon (stop_in_gtf is NA) does not change
  expect_identical(cds_without_stop(cases[[1]]$cds, 0), cases[[1]]$cds)
  expect_identical(cds_without_stop(GenomicRanges::GRanges(), 3), GenomicRanges::GRanges())
  # Before, the reference of ORFs_spl_feat_longest kept the stop codon, so the last
  # exon of an ORF that ends at the annotated stop codon was up_lastCDS
  orf <- gr(c(100, 300), c(200, 347), "+")
  expect_identical(spl_5to3(annotate_splicing(orf, cds_without_stop(cases[[1]]$cds, 3))),
                   c("same_firstCDS;same_3ss", "same_5ss;same_lastCDS"))
  expect_identical(spl_5to3(annotate_splicing(gr(100, 199, "+"), cds_without_stop(cases[[2]]$cds, 3))),
                   "same_5monoCDS;same_3monoCDS")
})
