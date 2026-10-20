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

# The tables of the results that the three files don't have
result_tables <- c("ORFs_feat", "ORFs_spl_feat_longest", "ORFs_spl_feat_maxORF",
                   "ORFs_txs_feats", "ORFs_readthroughs")
# Writes these tables of the results of `prefix` as TSV files, as run_ORFquant() writes
# ORFs_tx, with the name of each range and, for ORFs_feat, the number of its list
# element; the paths of the files
write_result_tables <- function(prefix) {
  res <- get(load(paste0(prefix, "_final_ORFquant_results")))
  files <- stats::setNames(paste0(prefix, "_", result_tables, ".tsv"), result_tables)
  for (t in result_tables) {
    x <- res[[t]]
    element <- NULL
    if (is(x, "GRangesList")) {
      element <- rep(seq_along(x), lengths(x))
      x <- unlist(x, use.names = FALSE)
    }
    df <- ORFs_tx_as_table(x)
    if (!is.null(names(x))) df <- cbind(name = names(x), df)
    if (!is.null(element)) df <- cbind(element = element, df)
    utils::write.table(df, files[[t]], sep = "\t", quote = FALSE, row.names = FALSE)
  }
  files
}
golden_tables <- function() {
  stats::setNames(test_path("golden", paste0("expected_", result_tables, ".tsv")), result_tables)
}
# Checks that the tables of `new` have the content of those of `old`, the numbers to 6 digits
expect_same_tables <- function(new, old) {
  for (t in names(old)) {
    expect_equal(read_tsv(new[[t]]), read_tsv(old[[t]]), tolerance = 1e-6, info = t)
  }
}
