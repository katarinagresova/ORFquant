# The command-line script that the RiboVerse skill runs: what it checks before
# it starts the pipeline. The pipeline itself is in test-pipeline.R.

cli_script <- system.file("scripts", "run_orfquant.R", package = "ORFquant")

# The exit status and the output (standard output and error) of the script
run_cli <- function(...) {
  out <- suppressWarnings(system2(rscript, shQuote(c(cli_script, ...)),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = paste(out, collapse = "\n"))
}

# The options of a run that is correct, as a vector of arguments. An argument
# of this function replaces an option, and NULL removes it.
cli_args <- function(...) {
  opts <- list(gtf = example_file("chr22_example.gtf.gz"),
               fasta = example_file("chr22_example.fa.gz"),
               bam = example_file("chr22_example.bam"),
               offsets = example_file("chr22_example_cutoffs.tsv"),
               outdir = tempfile("orfquant_cli_"))
  opts <- utils::modifyList(opts, list(...))
  as.vector(rbind(paste0("--", names(opts)), unlist(opts)))
}

test_that("--help prints the usage and succeeds", {
  r <- run_cli("--help")
  expect_equal(r$status, 0L)
  expect_match(r$output, "^Usage: Rscript run_orfquant.R")
  expect_match(r$output, "--annotation FILE", fixed = TRUE)
})

test_that("without arguments it prints the usage and fails", {
  r <- run_cli()
  expect_equal(r$status, 1L)
  expect_match(r$output, "^Usage: Rscript run_orfquant.R")
})

test_that("it stops with a message on wrong options, before it runs a step", {
  expect_cli_error <- function(args, message) {
    r <- run_cli(args)
    expect_equal(r$status, 1L, info = r$output)
    expect_match(r$output, message, fixed = TRUE)
    expect_false(grepl("failed", r$output, fixed = TRUE))   # no step ran
  }
  gtf <- example_file("chr22_example.gtf.gz")
  expect_cli_error(c(cli_args(), "--nonsense", "1"), "unknown option --nonsense")
  expect_cli_error(cli_args(fasta = NULL), "missing --fasta")
  expect_cli_error(cli_args(bam = NULL), "missing --bam")
  expect_cli_error(cli_args(annotation = gtf), "give either --gtf and --fasta, or --annotation")
  expect_cli_error(cli_args(gtf = NULL, fasta = NULL),
                   "give either --gtf and --fasta, or --annotation")
  expect_cli_error(cli_args(bam = "no_such_file.bam"), "--bam: file not found: no_such_file.bam")
  expect_cli_error(c(cli_args(), "--cores", "0"), "--cores must be a whole number above 0, not 0")
  expect_cli_error(c(cli_args(), "--cores", "two"),
                   "--cores must be a whole number above 0, not two")
})
