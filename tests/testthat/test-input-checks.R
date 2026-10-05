# Input that cannot work must stop a run at once, with a message that says what
# is wrong. These checks need no annotation, so they run before any long step.

work <- tempfile("orfquant_input_")
dir.create(work)
gtf_example <- readLines(gzfile(example_file("chr22_example.gtf.gz")))
fasta <- example_file("chr22_example.fa.gz")
write_gtf <- function(lines, name) {
  f <- file.path(work, name)
  writeLines(lines, f)
  f
}
# Runs prepare_annotation_files() on a GTF file and gives its error message
annotation_error <- function(gtf, genome_seq = fasta, out = file.path(work, "annotation")) {
  tryCatch({
    suppressMessages(prepare_annotation_files(annotation_directory = out, gtf_file = gtf,
                                              genome_seq = genome_seq, forge_BSgenome = FALSE))
    NA_character_
  }, error = function(e) conditionMessage(e))
}

test_that("prepare_annotation_files() names the files that do not exist, and writes nothing", {
  out <- file.path(work, "annotation_none")
  expect_match(annotation_error(file.path(work, "no.gtf"), out = out), "don't exist.*no\\.gtf")
  expect_match(annotation_error(write_gtf(gtf_example, "ok.gtf"), file.path(work, "no.fa"), out = out),
               "don't exist.*no\\.fa")
  expect_false(dir.exists(out))
})

test_that("prepare_annotation_files() asks for a genome sequence", {
  expect_error(prepare_annotation_files(annotation_directory = file.path(work, "annotation_none"),
                                        gtf_file = write_gtf(gtf_example, "ok.gtf")),
               "genome_seq")
})

test_that("prepare_annotation_files() refuses a GTF file without exons, or without ids", {
  expect_match(annotation_error(write_gtf(character(), "empty.gtf")), "empty or has no exon lines")
  expect_match(annotation_error(write_gtf(gtf_example[!grepl("\texon\t", gtf_example)], "noexon.gtf")),
               "no exon lines")
  expect_match(annotation_error(write_gtf(sub('transcript_id "[^"]*"; ', "", gtf_example), "notx.gtf")),
               "no transcript_id or gene_id")
  expect_match(annotation_error(write_gtf(sub('gene_id "[^"]*"; ', "", gtf_example), "nogene.gtf")),
               "no transcript_id or gene_id")
})

test_that("prepare_annotation_files() says when the chromosomes of the GTF are not in the genome", {
  msg <- annotation_error(write_gtf(sub("^chr22", "22", gtf_example), "nochr.gtf"))
  expect_match(msg, "Reading the GTF file .* failed")
  expect_match(msg, "chromosome names must be those of the genome sequence \\(e\\.g\\. chr22\\)")
})

test_that("prepare_for_ORFquant() and run_ORFquant() refuse wrong n_cores", {
  expect_error(prepare_for_ORFquant(annotation_file = "a", bam_file = "b", path_to_rl_cutoff_file = "c",
                                    n_cores = 0), "n_cores must be a number")
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a", n_cores = 0),
               "n_cores must be a number")
  expect_error(run_ORFquant(for_ORFquant_file = "f", annotation_file = "a"),
               "n_cores must be a number")
})

test_that("prepare_for_ORFquant() names the files that do not exist", {
  msg <- tryCatch(prepare_for_ORFquant(annotation_file = file.path(work, "no_Rannot"),
                                       bam_file = file.path(work, "no.bam"),
                                       path_to_rl_cutoff_file = file.path(work, "no.tsv")),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "don't exist")
  for (f in c("no_Rannot", "no.bam", "no.tsv")) expect_match(msg, f, fixed = TRUE)
})
