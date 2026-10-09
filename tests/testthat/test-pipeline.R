# The three steps on the example data in inst/extdata (7 genes of human chr22):
# with 1 core, with 2 cores and with the command-line script. All must give the
# same results, and the results of golden/, which are what ORFquant gave when
# the files were made.
#
# If a change that must not change results makes this test fail, the change is
# wrong. If the results must change, make the golden files again and review
# their diff, from the package directory:
#   ORFQUANT_UPDATE_GOLDEN=1 Rscript -e 'testthat::test_local()'

skip_on_os("windows")   # run_ORFquant() forks its workers

work <- tempfile("orfquant_pipeline_")
dir.create(work)
runs <- new.env()   # the prefix of the files of each run
quiet <- function(expr) invisible(utils::capture.output(suppressMessages(suppressWarnings(expr))))
# run_ORFquant() leaves its annotation in two global variables
withr::defer(rm(list = intersect(c("GTF_annotation", "genome_seq"), ls(globalenv())),
                envir = globalenv()), teardown_env())

# Prepares the P-sites of the example and runs ORFquant on them; the prefix
run_example <- function(name, n_cores, chunk_size) {
  dir.create(file.path(work, name))
  prefix <- file.path(work, name, "example")
  quiet(prepare_for_ORFquant(annotation_file = runs$annotation,
                             bam_file = example_file("chr22_example.bam"),
                             path_to_rl_cutoff_file = example_file("chr22_example_cutoffs.tsv"),
                             dest_name = prefix, n_cores = n_cores, chunk_size = chunk_size))
  quiet(run_ORFquant(for_ORFquant_file = paste0(prefix, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = n_cores,
                     prefix = prefix, interactive = FALSE))
  prefix
}

test_that("the annotation of the example is prepared", {
  quiet(prepare_annotation_files(annotation_directory = work,
                                 gtf_file = example_file("chr22_example.gtf.gz"),
                                 genome_seq = example_file("chr22_example.fa.gz"),
                                 forge_BSgenome = FALSE))
  runs$annotation <- file.path(work, "chr22_example.gtf.gz_Rannot")
  expect_true(file.exists(runs$annotation))
})

test_that("the example runs with 1 core", {
  skip_if(is.null(runs$annotation))
  runs$one_core <- run_example("one_core", n_cores = 1, chunk_size = 5e6)
  expect_true(all(file.exists(unlist(result_files(runs$one_core)))))
})

test_that("the example runs with 2 cores, and 6 chunks of alignments", {
  skip_if(is.null(runs$annotation))
  # 525 alignments in the BAM file: the P-sites of 6 chunks are calculated
  # in forked processes
  runs$two_cores <- run_example("two_cores", n_cores = 2, chunk_size = 100)
  expect_true(all(file.exists(unlist(result_files(runs$two_cores)))))
})

test_that("the P-sites do not depend on the number of cores or of chunks", {
  skip_if(is.null(runs$one_core) || is.null(runs$two_cores))
  expect_identical(get(load(paste0(runs$one_core, "_for_ORFquant"))),
                   get(load(paste0(runs$two_cores, "_for_ORFquant"))))
})

test_that("the command-line script runs, with 2 cores and the annotation", {
  skip_if(is.null(runs$annotation))
  # Its steps run in new R processes, which need ORFquant to be installed
  loads <- system2(rscript, c("-e", shQuote("quit(status = !requireNamespace('ORFquant', quietly = TRUE))")),
                   stdout = FALSE, stderr = FALSE) == 0
  skip_if_not(loads, "ORFquant is not installed in the library of Rscript")
  out <- system2(rscript,
                 shQuote(c(system.file("scripts", "run_orfquant.R", package = "ORFquant"),
                           "--annotation", runs$annotation,
                           "--bam", example_file("chr22_example.bam"),
                           "--offsets", example_file("chr22_example_cutoffs.tsv"),
                           "--outdir", file.path(work, "cli"), "--cores", "2",
                           "--sample", "example")),
                 stdout = TRUE, stderr = TRUE)
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  if (is.null(attr(out, "status"))) {
    runs$cli <- file.path(work, "cli", "example")
    expect_true(all(file.exists(unlist(result_files(runs$cli)))))
  }
})

test_that("1 core, 2 cores and the command-line script give the same results", {
  skip_if(is.null(runs$one_core) || is.null(runs$two_cores))
  expect_same_results(result_files(runs$two_cores), result_files(runs$one_core))
  skip_if(is.null(runs$cli), "the command-line script did not run")
  expect_same_results(result_files(runs$cli), result_files(runs$one_core))
})

test_that("the results are those of the golden files", {
  skip_if(is.null(runs$one_core))
  new <- result_files(runs$one_core)
  if (nzchar(Sys.getenv("ORFQUANT_UPDATE_GOLDEN"))) {
    file.copy(new$tsv, test_path("golden", "expected_ORFs.tsv"), overwrite = TRUE)
    file.copy(new$fasta, test_path("golden", "expected_proteins.fasta"), overwrite = TRUE)
    gz <- gzfile(test_path("golden", "expected_ORFs.gtf.gz"), "w")
    writeLines(gtf_lines(new$gtf), gz)
    close(gz)
  }
  expect_same_results(new, golden_files())
  # Golden files made from an empty result would pass the comparison
  expect_gt(nrow(read_tsv(new$tsv)), 0)
})

test_that("the output dictionary documents every column and category of the results", {
  skip_if(is.null(runs$one_core))
  rd <- tryCatch(tools::Rd_db("ORFquant")[["ORFquant_output.Rd"]], error = function(e) NULL)
  skip_if(is.null(rd), "ORFquant is not installed, so it has no help pages")
  text <- paste(as.character(rd), collapse = "")   # as.character() gives the tokens of the Rd file
  documented <- function(x) {
    vapply(x, function(i) grepl(paste0("\\code{", i, "}"), text, fixed = TRUE), TRUE)
  }
  orfs <- read_tsv(result_files(runs$one_core)$tsv)
  columns <- setdiff(colnames(orfs), c("seqnames", "start", "end", "width", "strand"))
  expect_true(all(documented(columns)), info = paste("not documented:",
                                                    paste(columns[!documented(columns)], collapse = ", ")))
  for (category in c("ORF_category_Tx", "ORF_category_Tx_compatible", "ORF_category_Gen")) {
    values <- unique(na.omit(orfs[[category]]))
    expect_true(all(documented(values)), info = paste(category, "not documented:",
                                                      paste(values[!documented(values)], collapse = ", ")))
  }
})

test_that("plot_orfquant_locus() plots a gene of the example without warnings", {
  skip_if(is.null(runs$one_core))
  for (pkg in c("Gviz", "lemon", "dplyr", "GenomeInfoDb")) skip_if_not_installed(pkg)
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  plotfile <- file.path(work, "locus.pdf")
  # before, ggplot2 warned that qplot() is deprecated, and the pdf device that
  # the font width of the tabs in the title of the P-site track is unknown
  expect_no_warning(utils::capture.output(suppressMessages(
    plot_orfquant_locus(locus = res$ORFs_tx$gene_id[1], orfquant_results = res, plotfile = plotfile))))
  expect_true(file.exists(plotfile))
})

# Input that cannot work stops the run before the BAM file is read, or before
# the regions are quantified; these use the annotation and P-sites made above

# The error message of a call, and the warnings that came before it
refusal <- function(expr) {
  warned <- character()
  msg <- tryCatch(withCallingHandlers(suppressMessages(expr),
                                      warning = function(w) {
                                        warned <<- c(warned, conditionMessage(w))
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) conditionMessage(e))
  list(message = if (is.character(msg)) msg else NA_character_, warnings = warned)
}
# The example annotation with changed seqinfo
changed_annotation <- function(change) {
  ann <- get(load(runs$annotation))
  ann$seqinfo <- change(ann$seqinfo)
  f <- tempfile("changed_Rannot_", work)
  save(ann, file = f)
  f
}
prepare_with <- function(cutoffs, annotation = runs$annotation, dest = file.path(work, "refused")) {
  refusal(prepare_for_ORFquant(annotation_file = annotation,
                               bam_file = example_file("chr22_example.bam"),
                               path_to_rl_cutoff_file = cutoffs, dest_name = dest))
}
write_cutoffs <- function(lines) {
  f <- tempfile("cutoffs_", work, fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("prepare_for_ORFquant() refuses an output directory that does not exist", {
  skip_if(is.null(runs$annotation))
  msg <- prepare_with(example_file("chr22_example_cutoffs.tsv"), dest = file.path(work, "no_dir", "sample"))$message
  expect_match(msg, "Cannot write the output files")
})

test_that("prepare_for_ORFquant() refuses a wrong cutoff table", {
  skip_if(is.null(runs$annotation))
  head <- "read_length\tcutoff\tcompartment"
  expect_match(prepare_with(write_cutoffs(head))$message, "has no rows")
  expect_match(prepare_with(write_cutoffs(c(head, "twenty\t12\tnucl")))$message, "must be numbers")
  out <- prepare_with(write_cutoffs(c(head, "28\t12\tmito")))
  expect_match(out$message, "No compartment of the rl_cutoff file is 'nucl'")
  expect_match(out$message, "mito")
})

test_that("prepare_for_ORFquant() refuses a BAM file with other chromosomes than the annotation", {
  skip_if(is.null(runs$annotation))
  other <- changed_annotation(function(si) {
    GenomeInfoDb::seqlevels(si) <- paste0("other", seq_along(GenomeInfoDb::seqlevels(si)))
    si
  })
  out <- prepare_with(write_cutoffs(c("read_length\tcutoff\tcompartment", "28\t12\tnucl", "29\t12\tmito")),
                      annotation = other)
  expect_match(out$message, "No chromosome of the BAM file is in the annotation")
  expect_match(out$message, "chr22")
  # the compartment that is not used is a warning, not an error
  expect_match(out$warnings, "mito")
})

test_that("prepare_for_ORFquant() refuses a cutoff table that gives no P-sites", {
  skip_if(is.null(runs$annotation))
  msg <- prepare_with(write_cutoffs(c("read_length\tcutoff\tcompartment", "99\t12\tnucl")))$message
  expect_match(msg, "No P-sites were found")
  expect_false(file.exists(file.path(work, "refused_for_ORFquant")))
})

test_that("run_ORFquant() refuses P-sites that are not on the genome of the annotation", {
  skip_if(is.null(runs$one_core))
  longer <- changed_annotation(function(si) {
    GenomeInfoDb::seqlengths(si) <- GenomeInfoDb::seqlengths(si) + 1L
    si
  })
  msg <- refusal(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                              annotation_file = longer, n_cores = 1,
                              prefix = file.path(work, "refused_run"), interactive = FALSE))$message
  expect_match(msg, "are not on the genome of the annotation")
})

test_that("run_ORFquant() refuses an output directory that does not exist", {
  skip_if(is.null(runs$one_core))
  msg <- refusal(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                              annotation_file = runs$annotation, n_cores = 1,
                              prefix = file.path(work, "no_dir", "sample"), interactive = FALSE))$message
  expect_match(msg, "Cannot write the output files")
})

test_that("prepare_for_ORFquant() sorts the P-sites of bigWig files by the chromosomes of the annotation", {
  skip_if(is.null(runs$one_core))
  # The annotation has chr22 and then chr1; the bigWig files have chr1 first
  annotation <- changed_annotation(function(si) {
    Seqinfo::seqlevels(si) <- c(Seqinfo::seqlevels(si), "chr1")
    Seqinfo::seqlengths(si)[["chr1"]] <- Seqinfo::seqlengths(si)[["chr22"]]
    si
  })
  ps <- get(load(paste0(runs$one_core, "_for_ORFquant")))$P_sites_all
  in_bw <- GenomicRanges::GRanges(rep(c("chr1", "chr22"), each = length(ps)), rep(IRanges::ranges(ps), 2),
                                  rep(BiocGenerics::strand(ps), 2), score = rep(ps$score, 2),
                                  seqinfo = get(load(annotation))$seqinfo[c("chr1", "chr22")])
  bw <- file.path(work, c("plus.bw", "minus.bw"))
  rtracklayer::export(in_bw[BiocGenerics::strand(in_bw) == "+"], bw[1], format = "bigWig")
  rtracklayer::export(in_bw[BiocGenerics::strand(in_bw) == "-"], bw[2], format = "bigWig")
  dest <- file.path(work, "bigwig")
  quiet(prepare_for_ORFquant(annotation_file = annotation, bam_file = example_file("chr22_example.bam"),
                             path_to_P_sites_plus_bw = bw[1], path_to_P_sites_minus_bw = bw[2], dest_name = dest))
  got <- get(load(paste0(dest, "_for_ORFquant")))$P_sites_all
  expect_false(BiocGenerics::is.unsorted(got))
  expected <- in_bw
  Seqinfo::seqlevels(expected) <- c("chr22", "chr1")
  expected <- BiocGenerics::sort(expected)
  expected$score <- as.numeric(expected$score)   # the scores of bigWig files are double
  expect_identical(as.data.frame(got), as.data.frame(expected))
})

test_that("run_ORFquant() gives no ORFs for a gene when stn.orf_quant.cutoff_P_sites removes all its ORFs", {
  skip_if(is.null(runs$one_core))
  # Before, the run stopped with "1 elements in value to replace 0 elements"
  # in the first region where the cutoff removed all ORFs
  dir.create(file.path(work, "cutoff"))
  prefix <- file.path(work, "cutoff", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = prefix,
                     stn.orf_quant.cutoff_P_sites = 30, interactive = FALSE))
  orfs <- read_tsv(result_files(prefix)$tsv)
  expect_true(all(orfs$P_sites >= 30))
  # The genes with an ORF of 30 P-sites or more without the cutoff keep ORFs;
  # in the example, some genes have no such ORF
  all_orfs <- read_tsv(result_files(runs$one_core)$tsv)
  max_P_sites <- tapply(all_orfs$P_sites, all_orfs$gene_id, max)
  expect_true(any(max_P_sites < 30))
  expect_setequal(unique(orfs$gene_id), names(max_P_sites)[max_P_sites >= 30])
})

test_that("run_ORFquant() finds ORFs without a start codon with stn.orf_find.nostarts = TRUE", {
  skip_if(is.null(runs$one_core))
  # Before, stn.orf_find.nostarts had no effect
  dir.create(file.path(work, "nostarts"))
  prefix <- file.path(work, "nostarts", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = prefix,
                     stn.orf_find.nostarts = TRUE, interactive = FALSE))
  orfs <- read_tsv(result_files(prefix)$tsv)
  no_start <- orfs[!startsWith(orfs$Protein, "M"), ]
  expect_gt(nrow(no_start), 0)
  # In the example, each is the part of the frame before the start codon of an ORF
  expect_true(all(paste(no_start$seqnames, no_start$end + 1) %in% paste(orfs$seqnames, orfs$start)))
})
