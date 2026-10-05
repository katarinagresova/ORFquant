# A file of the example data in inst/extdata
example_file <- function(f) system.file("extdata", f, package = "ORFquant")

# The Rscript of the R that runs the tests, to run the command-line script
rscript <- file.path(R.home("bin"), "Rscript")

# The files ORFquant writes with prefix `prefix`, and the golden files
result_files <- function(prefix) {
  list(tsv = paste0(prefix, "_Detected_ORFs.tsv"),
       fasta = paste0(prefix, "_Protein_sequences.fasta"),
       gtf = paste0(prefix, "_Detected_ORFs.gtf"))
}
golden_files <- function() {
  list(tsv = test_path("golden", "expected_ORFs.tsv"),
       fasta = test_path("golden", "expected_proteins.fasta"),
       gtf = test_path("golden", "expected_ORFs.gtf.gz"))
}

read_tsv <- function(f) utils::read.delim(f, stringsAsFactors = FALSE)
# The lines of a GTF file without its header, which has the date
gtf_lines <- function(f) {
  lines <- readLines(f)
  lines[!startsWith(lines, "#")]
}

# Checks that the three files of `new` have the content of those of `old`:
# the numbers of the TSV file to 6 digits, the others exactly
expect_same_results <- function(new, old) {
  expect_equal(read_tsv(new$tsv), read_tsv(old$tsv), tolerance = 1e-6)
  expect_identical(readLines(new$fasta), readLines(old$fasta))
  expect_identical(gtf_lines(new$gtf), gtf_lines(old$gtf))
}
