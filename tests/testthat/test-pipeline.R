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
  # txdbmaker warns that it ignores the phase of the stop_codon lines, and that the
  # TxDb has no genome version
  quiet(prepare_annotation_files(annotation_directory = work,
                                 gtf_file = example_file("chr22_example.gtf.gz"),
                                 genome_seq = example_file("chr22_example.fa.gz"),
                                 forge_BSgenome = FALSE),
        expected = c('The "phase" metadata column contains non-NA values',
                     "genome version information is not available"))
  runs$annotation <- file.path(work, "chr22_example.gtf.gz_Rannot")
  expect_true(file.exists(runs$annotation))
})

test_that("the example runs with 1 core", {
  skip_if(is.null(runs$annotation))
  rm(list = intersect(c("GTF_annotation", "genome_seq"), ls(globalenv())), envir = globalenv())
  runs$one_core <- run_example("one_core", n_cores = 1, chunk_size = 5e6)
  expect_true(all(file.exists(unlist(result_files(runs$one_core)))))
  # prepare_for_ORFquant() and run_ORFquant() leave the annotation in the global variables, as before
  expect_true(all(c("GTF_annotation", "genome_seq") %in% ls(globalenv())))
})

test_that("the example runs with 2 cores, and 6 chunks of alignments", {
  skip_if(is.null(runs$annotation))
  # 525 alignments in the BAM file: the P-sites of 6 chunks are calculated
  # in forked processes
  runs$two_cores <- run_example("two_cores", n_cores = 2, chunk_size = 100)
  expect_true(all(file.exists(unlist(result_files(runs$two_cores)))))
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
  one_core <- write_result_tables(runs$one_core)
  expect_same_results(result_files(runs$two_cores), result_files(runs$one_core))
  expect_same_tables(write_result_tables(runs$two_cores), one_core)
  skip_if(is.null(runs$cli), "the command-line script did not run")
  expect_same_results(result_files(runs$cli), result_files(runs$one_core))
  expect_same_tables(write_result_tables(runs$cli), one_core)
})

test_that("the results are those of the golden files", {
  skip_if(is.null(runs$one_core))
  new <- result_files(runs$one_core)
  tables <- write_result_tables(runs$one_core)
  if (nzchar(Sys.getenv("ORFQUANT_UPDATE_GOLDEN"))) {
    file.copy(new$tsv, test_path("golden", "expected_ORFs.tsv"), overwrite = TRUE)
    file.copy(new$fasta, test_path("golden", "expected_proteins.fasta"), overwrite = TRUE)
    gz <- gzfile(test_path("golden", "expected_ORFs.gtf.gz"), "w")
    writeLines(gtf_lines(new$gtf), gz)
    close(gz)
    file.copy(tables, golden_tables(), overwrite = TRUE)
  }
  expect_same_results(new, golden_files())
  expect_same_tables(tables, golden_tables())
  # Golden files made from an empty result would pass the comparison
  expect_gt(nrow(read_tsv(new$tsv)), 0)
})

# Before, an ORF could be alone in its group, or in two groups, and keep the P-sites that it
# shares with other ORFs: the TCN2 uORF ENST00000450638.5_250_336 had 22 P-sites, and the
# group of TCN2 had 171 P-sites on bins with 149
test_that("the P-sites of a group of ORFs that share features add up to the P-sites on their exonic bins", {
  skip_if(is.null(runs$one_core))
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  p_sites <- stats::setNames(res$ORFs_tx$P_sites, res$ORFs_tx$ORF_id_tr)
  for (feats in as.list(res$ORFs_feat)) {
    orfs <- lapply(as.list(feats$ORF_id_tr_selected), intersect, names(p_sites))
    groups <- list()
    for (x in orfs[lengths(orfs) > 0]) {
      joined <- vapply(groups, function(g) any(g %in% x), logical(1))
      groups <- c(groups[!joined], list(union(unlist(groups[joined]), x)))
    }
    for (g in groups) {
      on_bins <- feats$type == "E" & vapply(orfs, function(x) any(x %in% g), logical(1))
      expect_equal(sum(p_sites[g]), sum(feats$reads[on_bins]))
    }
  }
  expect_false("ENST00000450638.5_250_336" %in% names(p_sites))
})

test_that("with average_coverage, the coverage of the ORFs of a gene adds up to the coverage of its features", {
  skip_if(is.null(runs$one_core))
  dir.create(file.path(work, "average_coverage"))
  prefix <- file.path(work, "average_coverage", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = prefix,
                     stn.orf_quant.scaling = "average_coverage", interactive = FALSE))
  res <- get(load(paste0(prefix, "_final_ORFquant_results")))
  coverage <- stats::setNames(res$ORFs_tx$P_sites / width(res$ORFs_tx), res$ORFs_tx$ORF_id_tr)
  n <- 0
  for (feats in as.list(res$ORFs_feat)) {
    # an ORF below stn.orf_quant.cutoff_pct is not in the results
    if (!all(unlist(feats$ORF_id_tr_selected) %in% names(coverage))) next
    n <- n + 1
    # One factor for the ORFs of a region makes their coverage that of the features,
    # at most 1 (in the example, each such gene has a factor below 1). Before, the
    # reads of a junction were divided by its width and then by 60
    on_feats <- vapply(as.list(feats$ORF_id_tr_selected), function(x) sum(coverage[x]), numeric(1))
    expect_equal(sum(on_feats), sum(ifelse(feats$type == "J", feats$reads / 60, feats$reads / width(feats))))
  }
  expect_gt(n, 0)
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
  # Before, it left the option ucscChromosomeNames of Gviz at FALSE. A test of
  # test-input-checks.R sets it to FALSE with the code before the fix, so set it here
  withr::local_options(ucscChromosomeNames = TRUE)
  # before, ggplot2 warned that qplot() is deprecated, and the pdf device that
  # the font width of the tabs in the title of the P-site track is unknown
  expect_no_warning(utils::capture.output(suppressMessages(
    plot_orfquant_locus(locus = res$ORFs_tx$gene_id[1], orfquant_results = res, plotfile = plotfile))))
  expect_true(file.exists(plotfile))
  expect_true(getOption("ucscChromosomeNames"))
})

test_that("plot_ORFquant_results() plots the splice sites of the example with coverage files", {
  skip_if(is.null(runs$one_core))
  ps <- get(load(paste0(runs$one_core, "_for_ORFquant")))$P_sites_all
  cov <- file.path(work, c("coverage_plus.bw", "coverage_minus.bw"))
  rtracklayer::export(ps[strand(ps) == "+"], cov[1], format = "bigWig")
  rtracklayer::export(ps[strand(ps) == "-"], cov[2], format = "bigWig")
  # before, it stopped when a type of splice site had no regions, or regions on
  # one strand only, as a type of the example has (S24). Some groups of ORFs of the
  # example are too small for the distributions of the plots, and ggplot2 warns
  quiet(plot_ORFquant_results(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                              ORFquant_output_file = paste0(runs$one_core, "_final_ORFquant_results"),
                              annotation_file = runs$annotation,
                              coverage_file_plus = cov[1], coverage_file_minus = cov[2],
                              output_plots_path = file.path(work, "plots"), prefix = "example"),
        expected = c("Groups with fewer than two datapoints have been dropped",
                     "Computation failed in `stat_ydensity()`",
                     "no non-missing arguments to min", "no non-missing arguments to max"))
  expect_true(file.exists(file.path(work, "plots", "example_Meta_splicing_coverage.pdf")))
  plots <- get(load(file.path(work, "plots", "example_ORFquant_plots_RData")))
  # Before, all genes of ORFs_maxiso_tpm were in the bin NA
  expect_false(anyNA(plots$ORFs_maxiso_tpm$data$maxiso))
  expect_identical(levels(plots$ORFs_maxiso_tpm$data$maxiso), levels(plots$ORFs_maxiso$data$maxiso))
  # Before, a gene with 4 selected transcripts was in the bin "3-4", and the bins
  # were in the order of their numbers of genes
  sel <- plots$Sel_txs_genes_tpm$data
  bins <- c("1-3", "4-6", "7-9", ">9")
  expect_identical(as.character(sel$qnt), bins[findInterval(sel$elementNROWS.ch_txs_sel., c(1, 4, 7, 10))])
  expect_identical(levels(plots$Sel_txs_genes$data$qnt), intersect(bins, as.character(sel$qnt)))
  # Before, the file also held the annotation, the results and the P-sites, in the
  # environment of the plots. The plots do not need them
  expect_false(any(c("GTF_annotation", "ORFquant_results", "for_ORFquant") %in% ls(plots$ORFs_found$plot_env)))
  for (n in setdiff(names(plots), "Meta_splicing_coverage")) {
    expect_no_error(suppressWarnings(ggplot_build(plots[[n]])))
  }
  # Before, the log axes had no labels above 1000, and the y label of ORFs_pM was P-sites_pNpM
  expect_identical(plots$ORFs_found_ORFs_pM$labels$y, "ORFs_pM")
  y_columns <- c(ORFs_found_ORFs_pM = "ORFs_tx.ORFs_pM", ORFs_genes_tpm = "tpms",
                 ORFs_maxiso_tpm = "tpms_mult", Sel_txs_genes_tpm = "tpm")
  for (n in names(y_columns)) {
    y <- suppressWarnings(ggplot_build(plots[[n]]))$layout$panel_params[[1]]$y
    top <- max(as.numeric(y$get_labels()[!is.na(y$get_breaks())]))
    expect_gt(top * 10, max(plots[[n]]$data[[y_columns[[n]]]]), label = n)
  }
})

test_that("plot_ORFquant_results() plots the splice sites with bedGraph coverage files that end before a splice site", {
  skip_if(is.null(runs$one_core))
  ps <- get(load(paste0(runs$one_core, "_for_ORFquant")))$P_sites_all
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  spl <- unlist(res$ORFs_spl_feat_maxORF)
  cov <- file.path(work, c("coverage_plus.bedGraph", "coverage_minus.bedGraph"))
  # The P-sites before the middle of the splice features of each strand. Before, the
  # coverage of a bedGraph (it has no chromosome lengths) ended at its last P-site,
  # and a window around a splice site after it stopped the plots with "subscript
  # contains out-of-bounds ranges"
  for (i in 1:2) {
    p <- ps[strand(ps) == c("+", "-")[i]]
    cut <- stats::median(start(spl[strand(spl) == c("+", "-")[i]]))
    rtracklayer::export(p[end(p) < cut], cov[i], format = "bedGraph")
  }
  quiet(plot_ORFquant_results(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                              ORFquant_output_file = paste0(runs$one_core, "_final_ORFquant_results"),
                              annotation_file = runs$annotation,
                              coverage_file_plus = cov[1], coverage_file_minus = cov[2],
                              output_plots_path = file.path(work, "plots_bedgraph"), prefix = "example"),
        expected = c("Groups with fewer than two datapoints have been dropped",
                     "Computation failed in `stat_ydensity()`",
                     "no non-missing arguments to min", "no non-missing arguments to max"))
  expect_true(file.exists(file.path(work, "plots_bedgraph", "example_Meta_splicing_coverage.pdf")))
})

test_that("plot_ORFquant_results() plots the splice sites of an ORF exon that spans an intron of the reference", {
  skip_if(is.null(runs$one_core))
  load(paste0(runs$one_core, "_final_ORFquant_results"))
  ps <- get(load(paste0(runs$one_core, "_for_ORFquant")))$P_sites_all
  # The example has no event CDS_spanning (an ORF exon over an intron of the CDS of
  # its reference). Here an ORF gets one, from its second to its third exon
  spl <- ORFquant_results$ORFs_spl_feat_maxORF
  id <- "ENST00000402752.5_227_826"
  x <- spl[names(spl) == id]
  expect_length(x, 6)
  span <- GRanges(seqnames(x)[2], IRanges(start(x[2]), end(x[3])), strand(x)[2], seqinfo = seqinfo(spl))
  names(span) <- id
  span$ref <- GRangesList(granges(x[2:3], use.mcols = FALSE))
  span$spl_type <- "CDS_spanning"
  ORFquant_results$ORFs_spl_feat_maxORF <- c(spl, span)
  dir.create(file.path(work, "spanning"))
  res_file <- file.path(work, "spanning", "example_final_ORFquant_results")
  save(ORFquant_results, file = res_file)
  cov <- file.path(work, "spanning", c("coverage_plus.bw", "coverage_minus.bw"))
  rtracklayer::export(ps[strand(ps) == "+"], cov[1], format = "bigWig")
  rtracklayer::export(ps[strand(ps) == "-"], cov[2], format = "bigWig")
  # Before, the plot of these events gave facet_grid() the argument facets, which
  # ggplot2 4.0 no longer accepts, and the function stopped before it wrote a plot
  quiet(plot_ORFquant_results(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                              ORFquant_output_file = res_file,
                              annotation_file = runs$annotation,
                              coverage_file_plus = cov[1], coverage_file_minus = cov[2],
                              output_plots_path = file.path(work, "spanning", "plots"), prefix = "example"),
        expected = c("Groups with fewer than two datapoints have been dropped",
                     "Computation failed in `stat_ydensity()`",
                     "no non-missing arguments to min", "no non-missing arguments to max"))
  expect_true(file.exists(file.path(work, "spanning", "plots", "example_Meta_splicing_coverage.pdf")))
})

test_that("create_ORFquant_html_report() makes the HTML report of the plots of the example", {
  plots <- file.path(work, "plots", "example_ORFquant_plots_RData")
  skip_if_not(file.exists(plots))
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not installed")
  template <- system.file("rmd", package = "ORFquant")
  template_files <- list.files(template, all.files = TRUE, recursive = TRUE)
  sinks <- sink.number()
  withr::defer(while (sink.number() > sinks) sink())
  # A sample without coverage files has no plot of the splice sites
  no_coverage <- file.path(work, "no_coverage_ORFquant_plots_RData")
  list_ORFquant_plots <- get(load(plots))
  list_ORFquant_plots$Meta_splicing_coverage <- NULL
  save(list_ORFquant_plots, file = no_coverage)
  report <- file.path(work, "example_report.html")
  quiet(create_ORFquant_html_report(input_files = c(plots, no_coverage),
                                    input_sample_names = c("example", "no_coverage"),
                                    output_file = report))
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  # The input file of the sample, the sections of the template and their plots
  expect_match(html, normalizePath(plots), fixed = TRUE)
  for (section in c("Transcript selection", "Summary of detected ORFs", "ORFs per gene")) {
    expect_match(html, section, fixed = TRUE)
  }
  expect_match(html, "data:image/png;base64,", fixed = TRUE)
  # Before, section 3.5 was empty for a sample without coverage files. pandoc wraps
  # the lines of the HTML text
  expect_match(gsub("\\s+", " ", html), "made the plots of this sample without coverage files",
               fixed = TRUE)
  # Before, it rendered the template in the installed package, and it did not end
  # its sink() when rendering failed
  expect_identical(list.files(template, all.files = TRUE, recursive = TRUE), template_files)
  not_rdata <- file.path(work, "not_an_RData_file")
  writeLines("text", not_rdata)
  expect_error(suppressMessages(create_ORFquant_html_report(input_files = not_rdata,
                                                            input_sample_names = "none",
                                                            output_file = file.path(work, "no_report.html"))))
  expect_identical(sink.number(), sinks)
  # Before, a missing file stopped inside knitr, with an error about get()
  expect_error(create_ORFquant_html_report(input_files = file.path(work, "no_such_file"), input_sample_names = "none",
                                           output_file = file.path(work, "no_report.html")),
               "don't exist")
  expect_error(create_ORFquant_html_report(input_files = plots, input_sample_names = c("one", "two"),
                                           output_file = file.path(work, "no_report.html")),
               "same length")
})

# Input that cannot work stops run_ORFquant() before the regions are quantified;
# these use the annotation and P-sites made above. The tests of
# prepare_for_ORFquant() are in RiboseQC, which has its code now.

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

test_that("run_ORFquant() names the gene_id and gene_name values that are not in the annotation", {
  skip_if(is.null(runs$one_core))
  run_genes <- function(...) {
    refusal(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                         annotation_file = runs$annotation, n_cores = 1,
                         prefix = file.path(work, "unknown_genes"), interactive = FALSE, ...))
  }
  # Before, an unknown gene_id stopped with an S4 subsetting error
  known_id <- read_tsv(result_files(runs$one_core)$tsv)$gene_id[1]
  expect_match(run_genes(gene_id = c(known_id, "no_such_id"))$message, "1 gene_id values are not in the annotation.*no_such_id")
  # and an unknown gene_name was not used, without a message
  out <- run_genes(gene_name = "no_such_name")
  expect_match(out$warnings, "1 gene_name values are not in the annotation.*no_such_name")
  expect_match(out$message, "Incorrect gene_id \\| gene_name \\| genomic_region")
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

test_that("run_ORFquant() writes the protein FASTA file without the GTF file, and the GTF file without the FASTA file", {
  skip_if(is.null(runs$one_core))
  # Before, write_GTF_file = FALSE also turned off the FASTA file
  dir.create(file.path(work, "no_gtf"))
  dir.create(file.path(work, "no_fasta"))
  no_gtf <- file.path(work, "no_gtf", "example")
  no_fasta <- file.path(work, "no_fasta", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = no_gtf,
                     write_GTF_file = FALSE, interactive = FALSE))
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = no_fasta,
                     write_protein_fasta = FALSE, interactive = FALSE))
  expect_false(file.exists(result_files(no_gtf)$gtf))
  expect_identical(readLines(result_files(no_gtf)$fasta), readLines(result_files(runs$one_core)$fasta))
  expect_false(file.exists(result_files(no_fasta)$fasta))
  expect_identical(gtf_lines(result_files(no_fasta)$gtf), gtf_lines(result_files(runs$one_core)$gtf))
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

test_that("each CDS line of the GTF has the ORF_id of its ORF, also on a transcript with more than one ORF", {
  # Before, all CDS lines of a transcript had the ORF_id of its first ORF
  # (lcalviell/ORFquant#13). The example has no transcript with more than one
  # ORF, but the run of the test above has: the part before the start codon is
  # an ORF of the same transcript
  prefix <- file.path(work, "nostarts", "example")
  skip_if(!file.exists(result_files(prefix)$gtf))
  expect_true(anyDuplicated(read_tsv(result_files(prefix)$tsv)$transcript_id) > 0)
  gtf <- rtracklayer::import(result_files(prefix)$gtf, format = "gtf")
  cds <- gtf[gtf$type == "CDS"]
  ORFs_gen <- get(load(paste0(prefix, "_final_ORFquant_results")))$ORFs_gen
  expect_setequal(paste(as.character(cds), cds$ORF_id), paste(as.character(ORFs_gen), names(ORFs_gen)))
})

test_that("the splice annotation of an ORF outside the CDS does not depend on the ORFs before it", {
  skip_if(is.null(runs$one_core))
  # Before, annotate_ORFs() compared such an ORF with the CDS of the last ORF that
  # overlapped one, so with no CDS for the first ORF of a region. In the example,
  # ENST00000495655.2_175_255 overlaps no CDS and is the 2nd of 3 ORFs of its gene
  orf <- "ENST00000495655.2_175_255"
  quiet(load_annotation(runs$annotation, envir = globalenv()))
  for_ORFquant <- get(load(paste0(runs$one_core, "_for_ORFquant")))
  region <- GenomicRanges::reduce(unlist(GTF_annotation$txs_gene["ENSG00000184470.21"]))
  genetic_code <- Biostrings::getGeneticCode("SGC0")
  attributes(genetic_code)$alt_init_codons <- names(which(genetic_code == "M"))
  res <- NULL
  quiet(res <- ORFquant(region = region, for_ORFquant = for_ORFquant, genetic_code_region = genetic_code))
  expect_true(orf %in% names(res$ORFs_tx_position))
  # annotate_ORFs() again, with the ORF first or last
  spl_type_of_orf <- function(first) {
    ids <- names(res$ORFs_tx_position)
    ids <- if (first) c(orf, setdiff(ids, orf)) else c(setdiff(ids, orf), orf)
    res$ORFs_tx_position <- res$ORFs_tx_position[ids]
    res$ORFs_genomic_position <- res$ORFs_genomic_position[ids]
    annotated <- NULL
    quiet(annotated <- annotate_ORFs(res, GTF_annotation, genome_seq, region, genetic_code))
    annotated$ORFs_splice_feats$annotation_wrt_longest[[orf]]$spl_type
  }
  first <- spl_type_of_orf(first = TRUE)
  expect_identical(first, spl_type_of_orf(first = FALSE))
  # The reference is the CDS of the gene: its exons that the ORF does not overlap are missing
  expect_true(any(grepl("missing_CDS", first)))
})

test_that("ORFs_spl_feat_longest compares an ORF with the CDS of ref_id without its stop codon", {
  skip_if(is.null(runs$one_core))
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  # The ORFs with the start and the stop of the CDS of ref_id: TCN2 (+ strand) and TOP3B (- strand)
  exact <- res$ORFs_tx$ORF_id_tr[res$ORFs_tx$ORF_category_Gen == "exact_start_stop"]
  expect_setequal(exact, c("ENST00000698268.1_172_1479", "ENST00000357179.10_157_2742"))
  spl <- res$ORFs_spl_feat_longest
  for (id in exact) {
    types <- spl$spl_type[names(spl) == id]
    # Before, the CDS had its stop codon, and the ORF not, so the last exon was up_lastCDS
    expect_true(any(grepl("same_lastCDS", types)), info = id)
    expect_false(any(grepl("up_lastCDS|missing_CDS_3prime", types)), info = id)
  }
})

test_that("The ORF categories, NC_protein_isoform and ORFs_spl_feat_longest keep the last codon of a CDS that ends without a stop codon", {
  skip_if(is.null(runs$one_core))
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  loaded <- NULL
  quiet(loaded <- load_annotation(runs$annotation, envir = new.env()))
  for_ORFquant <- get(load(paste0(runs$one_core, "_for_ORFquant")))
  genetic_code <- Biostrings::getGeneticCode("SGC0")
  attributes(genetic_code)$alt_init_codons <- names(which(genetic_code == "M"))
  # The ORFs of the test above. Here the CDS of the gene that end at the stop codon of ref_id
  # (also the ORF's transcript and its compatible_tx) end without it, as a CDS with the GENCODE
  # tag cds_end_NF does, so the ORF ends at the end of the CDS. All of them: TOP3B has two
  # transcripts with this CDS, and the other one would become ref_id
  orfs <- res$ORFs_tx[res$ORFs_tx$ORF_category_Gen == "exact_start_stop"]
  expect_length(orfs, 2)
  expect_identical(orfs$compatible_tx, orfs$transcript_id)
  cds_end <- function(cds) if (as.vector(strand(cds))[1] == "+") max(end(cds)) else min(start(cds))
  for (k in seq_along(orfs)) {
    annotation <- loaded$GTF_annotation
    txs <- names(annotation$cds_txs)
    txs <- txs[annotation$trann$gene_id[match(txs, annotation$trann$transcript_id)] %in% orfs$gene_id[k]]
    txs <- txs[vapply(txs, function(t) cds_end(annotation$cds_txs[[t]]), numeric(1)) == cds_end(annotation$cds_txs[[orfs$ref_id[k]]])]
    expect_true(orfs$transcript_id[k] %in% txs)
    for (t in txs) {
      annotation$cds_txs[[t]] <- cds_without_stop(annotation$cds_txs[[t]], 3)
    }
    on_tx <- as.character(seqnames(annotation$cds_txs_coords)) %in% txs
    end(annotation$cds_txs_coords)[on_tx] <- end(annotation$cds_txs_coords)[on_tx] - 3
    region <- GenomicRanges::reduce(unlist(annotation$txs_gene[orfs$gene_id[k]]))
    out <- NULL
    quiet(out <- ORFquant(region = region, for_ORFquant = for_ORFquant, genetic_code_region = genetic_code,
                          annotation = annotation, genome_sequence = loaded$genome_seq))
    types <- out$ORFs_splice_feats$annotation_wrt_longest[[orfs$ORF_id_tr[k]]]$spl_type
    # Before, the last codon of this CDS was removed too, and the last exon was down_lastCDS
    expect_true(any(grepl("same_lastCDS", types)), info = orfs$ORF_id_tr[k])
    expect_false(any(grepl("down_lastCDS|new_lastCDS_3prime", types)), info = orfs$ORF_id_tr[k])
    # Before, C_extension, Alt3_stop and C
    o <- out$ORFs_tx_position[[orfs$ORF_id_tr[k]]]
    expect_identical(o$ref_id, orfs$ref_id[k])
    expect_identical(c(o$ORF_category_Tx, o$ORF_category_Tx_compatible, o$ORF_category_Gen, o$NC_protein_isoform),
                     c("ORF_annotated", "ORF_annotated", "exact_start_stop", "same"), info = orfs$ORF_id_tr[k])
  }
})

test_that("NMD_candidate is TRUE when Distance_to_lastExEx is more than 3, also on the compatible transcripts", {
  skip_if(is.null(runs$one_core))
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  orfs <- res$ORFs_tx
  expect_identical(unname(orfs$NMD_candidate), unname(orfs$Distance_to_lastExEx > 3))
  d <- as.integer(unlist(orfs$Distance_to_lastExEx_compatible_txs, use.names = FALSE))
  # Distances of 10 or more: as text, "10" > "3" is FALSE
  expect_true(any(d >= 10))
  expect_identical(unlist(orfs$NMD_candidate_compatible_txs, use.names = FALSE), d > 3)
})

test_that("ORFquant() takes the annotation that load_annotation() returns, without the global variables", {
  skip_if(is.null(runs$one_core))
  loaded <- NULL
  quiet(loaded <- load_annotation(runs$annotation, envir = globalenv()))
  expect_named(loaded, c("GTF_annotation", "genome_seq"))
  # load_annotation() also assigns the variables in its argument envir
  expect_identical(get("GTF_annotation", envir = globalenv()), loaded$GTF_annotation)
  expect_identical(get("genome_seq", envir = globalenv()), loaded$genome_seq)
  for_ORFquant <- get(load(paste0(runs$one_core, "_for_ORFquant")))
  region <- GenomicRanges::reduce(unlist(loaded$GTF_annotation$txs_gene["ENSG00000184470.21"]))
  genetic_code <- Biostrings::getGeneticCode("SGC0")
  attributes(genetic_code)$alt_init_codons <- names(which(genetic_code == "M"))
  with_globals <- NULL
  quiet(with_globals <- ORFquant(region = region, for_ORFquant = for_ORFquant, genetic_code_region = genetic_code))
  rm(list = c("GTF_annotation", "genome_seq"), envir = globalenv())
  without_globals <- NULL
  quiet(without_globals <- ORFquant(region = region, for_ORFquant = for_ORFquant, genetic_code_region = genetic_code,
                                    annotation = loaded$GTF_annotation, genome_sequence = loaded$genome_seq))
  expect_true(length(with_globals$ORFs_tx_position) > 0)
  expect_equal(without_globals, with_globals)
})

test_that("run_ORFquant() with genomic_region analyzes the regions it overlaps whole", {
  skip_if(is.null(runs$one_core))
  res <- get(load(paste0(runs$one_core, "_final_ORFquant_results")))
  # One nucleotide of the first ORF. In the example, each region has one gene
  one_nt <- GenomicRanges::resize(res$ORFs_gen[names(res$ORFs_gen) == res$ORFs_tx$ORF_id_tr[1]][1], 1)
  dir.create(file.path(work, "genomic_region"))
  prefix <- file.path(work, "genomic_region", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = runs$annotation, n_cores = 1, prefix = prefix,
                     genomic_region = one_nt, interactive = FALSE))
  orfs <- orfs_table(prefix)
  all_orfs <- orfs_table(runs$one_core)
  all_orfs <- all_orfs[all_orfs$gene_id == res$ORFs_tx$gene_id[1], ]
  expect_gt(nrow(all_orfs), 1)
  # The ORFs of the gene, with the values of the run on all genes, except ORFs_pM,
  # which sums to a million over the ORFs of a run
  columns <- setdiff(colnames(all_orfs), "ORFs_pM")
  expect_equal(orfs[columns], all_orfs[columns], ignore_attr = TRUE)
})

test_that("run_ORFquant() with unique_reads_only = TRUE finds ORFs by the P-sites of unique reads, with 1 core and with 2 cores", {
  skip_if(is.null(runs$one_core))
  # In the example, 3 of the 5 ORFs of DDT and RANBP1 have pval_uniq of 0.05 or more
  genes <- c("DDT", "RANBP1")
  all_orfs <- read_tsv(result_files(runs$one_core)$tsv)
  expect_true(any(all_orfs$pval_uniq[all_orfs$gene_name %in% genes] >= 0.05))
  dir.create(file.path(work, "unique_reads"))
  prefixes <- file.path(work, "unique_reads", c("one_core", "two_cores"))
  for (n_cores in 1:2) {
    quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                       annotation_file = runs$annotation, n_cores = n_cores, prefix = prefixes[n_cores],
                       gene_name = genes, unique_reads_only = TRUE, interactive = FALSE))
  }
  orfs <- read_tsv(result_files(prefixes[1])$tsv)
  expect_gt(nrow(orfs), 0)
  expect_true(all(orfs$pval_uniq < 0.05))
  expect_same_results(result_files(prefixes[2]), result_files(prefixes[1]))
  expect_same_tables(write_result_tables(prefixes[2]), write_result_tables(prefixes[1]))
})

test_that("ORFquant() with unique_reads = TRUE analyzes a region only with more than 4 positions of unique P-sites", {
  skip_if(is.null(runs$one_core))
  # Before, it analyzed a region with fewer too, if it had more than 4 positions of
  # P-sites of all reads. The region of DDT has 43 positions, 12 of unique P-sites
  loaded <- NULL
  quiet(loaded <- load_annotation(runs$annotation, envir = globalenv()))
  for_ORFquant <- get(load(paste0(runs$one_core, "_for_ORFquant")))
  region <- GenomicRanges::reduce(unlist(loaded$GTF_annotation$txs_gene["ENSG00000099977.16"]))
  genetic_code <- Biostrings::getGeneticCode("SGC0")
  attributes(genetic_code)$alt_init_codons <- names(which(genetic_code == "M"))
  # ORFquant() with the n positions of unique P-sites of the region that have the most reads
  with_unique <- function(n) {
    uniq <- for_ORFquant$P_sites_uniq
    in_region <- which(IRanges::overlapsAny(uniq, region))
    drop <- in_region[order(uniq$score[in_region], decreasing = TRUE)][-seq_len(n)]
    for_ORFquant$P_sites_uniq <- uniq[!seq_along(uniq) %in% drop]
    res <- NULL
    quiet(res <- ORFquant(region = region, for_ORFquant = for_ORFquant, genetic_code_region = genetic_code,
                          unique_reads = TRUE, annotation = loaded$GTF_annotation,
                          genome_sequence = loaded$genome_seq))
    res
  }
  expect_length(with_unique(4), 0)
  expect_gt(length(with_unique(5)$genomic_features), 0)
})

test_that("run_ORFquant() with canonical_start_only = TRUE leaves only the start codons of M, with 1 core and with 2 cores", {
  skip_if(is.null(runs$one_core))
  # Before, with 1 core, the alternative start codons of the genetic code stayed
  # start codons. The example gives the same results without them, so here
  # ORFquant() stops with the alternative start codons of the genetic code it gets
  local_mocked_bindings(ORFquant = function(region, for_ORFquant, genetic_code_region, ...) {
    stop("start codons: ", paste(attr(genetic_code_region, "alt_init_codons"), collapse = " "))
  })
  start_codons <- function(n_cores, canonical_start_only) {
    refusal(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                         annotation_file = runs$annotation, n_cores = n_cores,
                         prefix = file.path(work, "start_codons"), interactive = FALSE,
                         canonical_start_only = canonical_start_only))$message
  }
  for (n_cores in 1:2) {
    expect_match(start_codons(n_cores, TRUE), "start codons: ATG$", info = paste("n_cores =", n_cores))
    expect_match(start_codons(n_cores, FALSE), "start codons: TTG CTG$", info = paste("n_cores =", n_cores))
  }
})

test_that("with canonical_start_only = FALSE, run_ORFquant() compares the C-terminus of an ORF with that of the CDS without a start codon", {
  skip_if(is.null(runs$one_core))
  # Before, NC_protein_isoform read an alternative start codon at the start of the
  # last 10 codons of the CDS as M. Genetic code 11 has the amino acids of the
  # standard code and more alternative start codons, among them ATC, which starts
  # the last 10 codons of the CDS that the ORFs of PISD are compared with
  ann <- get(load(runs$annotation))
  ann$genetic_codes$genetic_code <- "11"
  annotation <- tempfile("code_11_Rannot_", work)
  save(ann, file = annotation)
  dir.create(file.path(work, "code_11"))
  prefix <- file.path(work, "code_11", "example")
  quiet(run_ORFquant(for_ORFquant_file = paste0(runs$one_core, "_for_ORFquant"),
                     annotation_file = annotation, n_cores = 1, prefix = prefix, gene_name = "PISD",
                     canonical_start_only = FALSE, interactive = FALSE))
  orfs <- orfs_table(prefix)
  all_orfs <- orfs_table(runs$one_core)
  all_orfs <- all_orfs[all_orfs$gene_name == "PISD", ]
  expect_identical(all_orfs$NC_protein_isoform, c("N", "N"))
  # The ORFs and values of the run with the standard code, except ORFs_pM, which sums
  # to a million over the ORFs of a run
  columns <- setdiff(colnames(all_orfs), "ORFs_pM")
  expect_equal(orfs[columns], all_orfs[columns], ignore_attr = TRUE)
})
