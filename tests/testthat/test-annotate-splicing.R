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
