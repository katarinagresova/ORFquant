# The helper functions that make the pipeline faster. Each one must give what
# the code it replaces gave, so most tests compare the two.

gr <- GenomicRanges::GRanges
rg <- IRanges::IRanges

test_that("translate_solve() equals translate() with if.fuzzy.codon = 'solve'", {
  set.seed(1)
  dna <- function(n, letters) {
    Biostrings::DNAString(paste(sample(letters, n, replace = TRUE), collapse = ""))
  }
  sequences <- list(dna(300, c("A", "C", "G", "T")),       # no fuzzy codon: the fast path
                    dna(300, c("A", "C", "G", "T", "N")),  # fuzzy codons: "solve"
                    Biostrings::DNAString("ATGGGNCTNANNTAG"))
  for (x in sequences) {
    expect_identical(translate_solve(x),
                     Biostrings::translate(x, if.fuzzy.codon = "solve"))
  }
  # A fuzzy codon gives an amino acid if all its codons give the same one
  expect_equal(as.character(translate_solve(Biostrings::DNAString("GGNCTNANN"))), "GLX")
  # Sets of sequences, and the other arguments of translate()
  set <- Biostrings::DNAStringSet(c(a = "ATGAAATAG", b = "ATGGGNTAG"))
  expect_equal(as.character(translate_solve(set)), c(a = "MK*", b = "MG*"))
  mito <- Biostrings::getGeneticCode("SGC1")
  expect_equal(as.character(translate_solve(Biostrings::DNAString("ATATGAAGA"),
                                            genetic.code = mito)), "MW*")
})

# GRangesList_fast() must give the result of GRangesList(), or the same error
same_as_GRangesList <- function(x) {
  outcome <- function(f) tryCatch(f(x), error = conditionMessage)
  expect_identical(outcome(GRangesList_fast), outcome(GenomicRanges::GRangesList))
}

test_that("GRangesList_fast() gives what GRangesList() gives", {
  # Columns, an element without ranges and columns, several strands
  x <- list(a = gr("chr1", rg(c(1, 10), width = 5), "+", score = c(1, 2), type = c("x", "y")),
            b = gr("chr2", rg(3, width = 4), "-", score = 3, type = "z"),
            c = gr(),
            d = gr("chr1", rg(7, width = 2), "*", score = 4, type = "w"))
  same_as_GRangesList(x)
  expect_equal(lengths(GRangesList_fast(x)), c(a = 2, b = 1, c = 0, d = 1))

  # The Seqinfo are merged in the order of merging them one by one
  si <- Seqinfo::Seqinfo
  y <- list(gr("chr2", rg(1, 5), seqinfo = si(c("chr2", "chr1"), c(200, 100)), v = 1),
            gr("chr3", rg(2, 6), seqinfo = si(c("chr3", "chr2"), c(300, 200)), v = 2),
            gr("chr4", rg(3, 7), seqinfo = si(c("chr1", "chr4"), c(100, 400)), v = 3))
  same_as_GRangesList(y)
  expect_equal(Seqinfo::seqlevels(GRangesList_fast(y)), c("chr2", "chr1", "chr3", "chr4"))

  # Seqinfo that don't merge give the error of GRangesList()
  z <- list(gr("chr1", rg(1, 5), seqinfo = si("chr1", 100), v = 1),
            gr("chr1", rg(2, 6), seqinfo = si("chr1", 200), v = 2))
  same_as_GRangesList(z)
  expect_error(GRangesList_fast(z), "incompatible seqlengths")
})

test_that("bind_GRanges() gives NULL where GRangesList() must do the work", {
  factors <- list(gr("chr1", rg(1, 5), v = factor("a")), gr("chr1", rg(2, 6), v = factor("b")))
  columns <- list(gr("chr1", rg(1, 2), a = 1), gr("chr1", rg(3, 4), b = 1))
  other_class <- list(gr("chr1", rg(1, 2), a = 1), rg(3, 4))
  for (x in list(factors, columns, other_class, list())) {
    expect_null(bind_GRanges(x))
  }
  same_as_GRangesList(factors)
  expect_identical(GRangesList_fast(list()), GenomicRanges::GRangesList(list()))
})

test_that("mc_lapply() is lapply() with each function in a forked process", {
  skip_on_os("windows")
  fs <- list(function() 1, function() "a", function() NULL)
  expected <- lapply(fs, function(f) f())
  expect_identical(mc_lapply(fs, 1), expected)
  expect_identical(mc_lapply(fs, 2), expected)

  # The warning of a forked process is given again, and its error stops
  expect_warning(res <- mc_lapply(list(function() { warning("careful"); 1 }, function() 2), 2),
                 "careful")
  expect_identical(res, list(1, 2))
  expect_error(suppressWarnings(mc_lapply(list(function() 1, function() stop("boom")), 2)), "boom")
  # A process that ends without a result, as when the system kills it for memory
  expect_error(suppressWarnings(mc_lapply(list(function() 1,
                                               function() system2("kill", c("-9", Sys.getpid()))), 2)),
               "ended without a result")
})

test_that("mc_job() evaluates in a forked process with more than 1 core", {
  skip_on_os("windows")
  expect_equal(mc_value(mc_job(1 + 1, 1)), 2)
  expect_equal(mc_value(mc_job(1 + 1, 2)), 2)
  expect_warning(res <- mc_value(mc_job({ warning("jobwarn"); 5 }, 2)), "jobwarn")
  expect_equal(res, 5)
  expect_error(mc_value(mc_job(stop("jobboom"), 2)), "jobboom")
  job <- mc_job({ Sys.sleep(2); 1 }, 2)
  expect_null(mc_value(job, wait = FALSE))   # not done yet
  expect_equal(mc_value(job), 1)
})

test_that("get_orfs() finds the ORFs of a transcript in the three frames", {
  # Frame 2 has ATG AAA TTT TAG, frame 1 has ATG CCC TGA, and the last ATG has
  # no stop codon. An ORF ends before its stop codon.
  seq <- Biostrings::DNAString("CCATGAAATTTTAGGGATGCCCTGAAATGGGTAA")
  orfs <- get_orfs("tx", seq, genetic_code_table = Biostrings::GENETIC_CODE)
  expect_named(orfs, c("frame_0", "frame_1", "frame_2"))
  coordinates <- vapply(orfs, function(g) {
    paste(sprintf("%d-%d", IRanges::start(g), IRanges::end(g)), collapse = ",")
  }, "")
  expect_equal(coordinates, c(frame_0 = "", frame_1 = "17-22", frame_2 = "3-11"))
  expect_equal(as.character(GenomicRanges::seqnames(orfs$frame_2)), "tx_frame_2")
})

test_that("get_orfs() reports every start codon, or one ORF for each stop codon", {
  seq <- Biostrings::DNAString("ATGATGAAATAG")
  all_starts <- get_orfs("tx", seq, genetic_code_table = Biostrings::GENETIC_CODE)$frame_0
  expect_equal(IRanges::start(all_starts), c(1, 4))
  expect_equal(IRanges::end(all_starts), c(9, 9))
  merged <- get_orfs("tx", seq, get_all_starts = FALSE,
                     genetic_code_table = Biostrings::GENETIC_CODE)$frame_0
  expect_equal(IRanges::start(merged), 1)
  expect_equal(IRanges::end(merged), 9)
})

test_that("get_reathr_seq() gives no readthrough region for an ORF without a stop codon after it", {
  # The ORF ends at the end of the transcript, as an ORF without a start codon
  # (stn.orf_find.nostarts = TRUE) can; its frame has one stop codon, before it.
  # Before, this stopped with "replacement has length zero"
  orf <- GenomicRanges::GRanges("tx", IRanges::IRanges(4, 15), "+")
  expect_length(get_reathr_seq("tx", orf, Biostrings::DNAString("TAAGCCGCCGCCGCC"), Biostrings::GENETIC_CODE), 0)
  # After an ORF before a stop codon, the readthrough regions start at the stop codon
  orf <- GenomicRanges::GRanges("tx", IRanges::IRanges(1, 6), "+")
  regions <- get_reathr_seq("tx", orf, Biostrings::DNAString("ATGGCCTAAGCCGCCTAGGCC"), Biostrings::GENETIC_CODE)
  expect_equal(IRanges::start(regions), c(7, 16))
})

test_that("from_tx_togen() gives each ORF its own genomic exons, under its ORF_id_tr", {
  # A transcript of 2 exons (30 nt) on each strand. The ORFs are not sorted,
  # and one of them crosses the intron.
  coordinates <- function(x) vapply(x, function(g) {
    paste(IRanges::start(g), IRanges::end(g), sep = "-", collapse = ",")
  }, "")
  orfs <- gr("tx", rg(c(13, 1, 4), c(27, 9, 18)), "+")
  orfs$ORF_id_tr <- paste("tx", IRanges::start(orfs), IRanges::end(orfs), sep = "_")
  expected <- list("+" = c(tx_13_27 = "203-217", tx_1_9 = "101-109", tx_4_18 = "104-110,201-208"),
                   "-" = c(tx_13_27 = "104-110,201-208", tx_1_9 = "212-220", tx_4_18 = "203-217"))
  for (s in c("+", "-")) {
    exons <- GenomicRanges::GRangesList(tx = gr("chr1", rg(c(101, 201), c(110, 220)), s))
    introns <- gr("chr1", rg(111, 200), s)
    res <- from_tx_togen(orfs, exons, introns)
    expect_equal(coordinates(res), expected[[s]])
    expect_true(all(unlist(GenomicRanges::strand(res), use.names = FALSE) == s))

    # mapFromTranscripts() drops an ORF outside the transcript. Then
    # from_tx_togen() stops; it does not give the next ORF's exons to it.
    outside <- gr("tx", rg(c(1, 22, 13), c(9, 33, 27)), "+")
    outside$ORF_id_tr <- paste("tx", IRanges::start(outside), IRanges::end(outside), sep = "_")
    expect_error(from_tx_togen(outside, exons, introns), "exactly one transcript of 'exons': tx_22_33$")
    # With the transcript name 2 times, each ORF maps 2 times. Before, the
    # second ORF got the second range of the first ORF, without an error.
    twice <- c(exons, GenomicRanges::GRangesList(tx = gr("chr1", rg(301, 330), s)))
    expect_error(from_tx_togen(orfs, twice, introns), "tx_13_27, tx_1_9, tx_4_18$")
  }
})
