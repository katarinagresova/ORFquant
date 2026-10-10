# The annotation of GTF files that are not GENCODE's. The GTF of the example
# (GENCODE's attribute names, all biotypes, no gene on two strands) is changed in
# the ways that other GTFs differ from it. The annotation must be that of the
# example, or the one that NEWS.md describes.

work <- tempfile("orfquant_gtf_")
dir.create(work)
gtf_example <- readLines(gzfile(example_file("chr22_example.gtf.gz")))
fasta <- example_file("chr22_example.fa.gz")
# run_ORFquant() leaves its annotation in two global variables
withr::defer(rm(list = intersect(c("GTF_annotation", "genome_seq"), ls(globalenv())),
                envir = globalenv()), teardown_env())

# GenomicFeatures reads a start or stop codon as a CDS too
coding_types <- c("CDS", "start_codon", "stop_codon")
line_type <- function(lines) sub("^[^\t]*\t[^\t]*\t([^\t]*)\t.*", "\\1", lines)
line_attr <- function(lines, key) {
  has <- grepl(paste0(key, ' "'), lines, fixed = TRUE)
  ifelse(has, sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", lines), NA_character_)
}
# The lines without the attribute `key`, on the `lines` selected by `where`
drop_attr <- function(lines, key, where = rep(TRUE, length(lines))) {
  lines[where] <- sub(paste0(key, ' "[^"]*"; ?'), "", lines[where])
  lines
}

# Writes `lines` as a GTF file, prepares its annotation, and gives the path of the
# annotation file; `annotate()` gives the annotation
annotation_of <- function(lines, name) {
  dir <- file.path(work, name)
  dir.create(dir)
  gtf <- file.path(dir, "genes.gtf")
  writeLines(lines, gtf)
  suppressMessages(suppressWarnings(
    prepare_annotation_files(annotation_directory = dir, gtf_file = gtf,
                             genome_seq = fasta, forge_BSgenome = FALSE)))
  file.path(dir, "genes.gtf_Rannot")
}
annotate <- function(lines, name) get(load(annotation_of(lines, name)))

ann_genc <- annotate(gtf_example, "genc")
tx_ids <- unique(stats::na.omit(line_attr(gtf_example, "transcript_id")))

test_that("the annotation of the example has one row for each transcript", {
  expect_equal(ann_genc$trann$transcript_id, tx_ids)
  expect_true(all(ann_genc$trann$transcript_biotype %in%
                    c("protein_coding", "protein_coding_CDS_not_defined", "retained_intron",
                      "nonsense_mediated_decay", "non_stop_decay")))
})

# Names of attributes, as in the GTFs of Ensembl, NCBI, UCSC, gffread and StringTie
test_that("ids, names and biotypes are read by the name of their attribute", {
  # other names of the attributes, also mixed with GENCODE's: those of Ensembl, NCBI and StringTie GTFs
  renames <- list(c(gene_type = "gene_biotype", gene_name = "gene_symbol", transcript_type = "transcript_biotype"),
                  c(gene_name = "gene", transcript_type = "transcript_biotype"),
                  c(gene_name = "ref_gene_name"))
  for (i in seq_along(renames)) {
    gtf <- gtf_example
    for (old in names(renames[[i]])) gtf <- sub(paste0(old, ' "'), paste0(renames[[i]][[old]], ' "'), gtf)
    expect_equal(annotate(gtf, paste0("renamed", i))$trann, ann_genc$trann, info = paste(renames[[i]], collapse = ", "))
  }
})

test_that("a biotype or a name that only some lines of a gene or transcript have is found", {
  # as in the GTF of gffread, which writes them on the transcript lines only
  other <- !line_type(gtf_example) %in% c("gene", "transcript")
  gtf <- gtf_example
  for (key in c("gene_type", "gene_name", "transcript_type")) gtf <- drop_attr(gtf, key, other)
  expect_equal(annotate(gtf, "tx_lines")$trann, ann_genc$trann)
})

test_that("a GTF without gene names has the gene name no_name", {
  trann <- annotate(drop_attr(gtf_example, "gene_name"), "no_names")$trann
  expect_equal(trann$gene_name, rep("no_name", nrow(trann)))
  expect_equal(trann$gene_biotype, ann_genc$trann$gene_biotype)
})

test_that("a GTF without CDS lines is refused", {
  msg <- tryCatch(annotation_of(gtf_example[!line_type(gtf_example) %in% coding_types], "no_cds"),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "no CDS lines")
})

# Biotypes
test_that("the biotype mRNA of NCBI is protein_coding", {
  gtf <- sub('transcript_type "protein_coding"', 'transcript_type "mRNA"', gtf_example)
  trann <- annotate(gtf, "mrna")$trann
  expect_false("mRNA" %in% trann$transcript_biotype)
  expect_equal(trann, ann_genc$trann)
})

test_that("without biotypes, a transcript with a CDS and its gene are protein_coding", {
  # RANBP1 has no CDS lines here (and no start and stop codons): it is not coding
  ranbp1 <- "ENSG00000099901.18"
  gtf <- drop_attr(drop_attr(gtf_example, "gene_type"), "transcript_type")
  gtf <- gtf[!(line_type(gtf) %in% coding_types & line_attr(gtf, "gene_id") == ranbp1)]
  trann <- annotate(gtf, "no_types")$trann
  has_cds <- trann$transcript_id %in% line_attr(gtf, "transcript_id")[line_type(gtf) == "CDS"]
  expect_equal(trann$transcript_biotype, ifelse(has_cds, "protein_coding", "no_type"))
  expect_equal(trann$gene_biotype, ifelse(trann$gene_id == ranbp1, "no_type", "protein_coding"))
  # also those that GENCODE calls nonsense_mediated_decay, which have a CDS
  nmd <- ann_genc$trann$transcript_biotype == "nonsense_mediated_decay" & ann_genc$trann$gene_id != ranbp1
  expect_true(any(nmd))
  expect_true(all(trann$transcript_biotype[nmd] == "protein_coding"))
})

# Genes that genes() drops, as the PAR genes of UCSC's GTFs, which are on chrX and chrY.
# DDT with a copy of its transcript ENST00000350608.7 on the other strand: the gene
# is on both strands, so genes() drops it
ddt <- "ENSG00000099977.16"
ddt_tx <- "ENST00000350608.7"
gtf_two_strands <- local({
  copy <- gtf_example[line_attr(gtf_example, "transcript_id") %in% ddt_tx]
  copy <- sub("\t-\t", "\t+\t", sub(ddt_tx, paste0(ddt_tx, "_rev"), copy, fixed = TRUE))
  c(gtf_example, copy)
})

test_that("a gene on both strands is kept in the annotation, and is not in $genes", {
  ann <- annotate(gtf_two_strands, "two_strands")
  expect_equal(ann$genes, ann_genc$genes[names(ann_genc$genes) != ddt])
  expect_true(ddt %in% names(ann$txs_gene))
  expect_equal(ann$trann$transcript_id, c(tx_ids, paste0(ddt_tx, "_rev")))
  # the regions of the gene have no gene
  locus <- unlist(range(ann$txs_gene[ddt]))
  regions <- ann$introns[ann$introns %over% locus]
  expect_gt(length(regions), 0)
  expect_true(all(lengths(regions$gene_id) == 0))
  # the other genes have the regions of the example
  expect_equal(ann$ncRNAs[!ann$ncRNAs %over% locus], ann_genc$ncRNAs[!ann_genc$ncRNAs %over% locus])
})

test_that("run_ORFquant() runs on a gene that genes() dropped", {
  skip_on_os("windows")   # run_ORFquant() forks its workers
  quiet <- function(expr) invisible(utils::capture.output(suppressMessages(suppressWarnings(expr))))
  annotation <- file.path(work, "two_strands", "genes.gtf_Rannot")
  prefix <- file.path(work, "two_strands", "example")
  quiet(prepare_for_ORFquant(annotation_file = annotation,
                             bam_file = example_file("chr22_example.bam"),
                             path_to_rl_cutoff_file = example_file("chr22_example_cutoffs.tsv"),
                             dest_name = prefix, n_cores = 1))
  golden <- read_tsv(golden_files()$tsv)
  golden <- golden[golden$gene_name == "DDT", ]
  expect_gt(nrow(golden), 0)
  for (by in list(list(gene_name = "DDT"), list(gene_id = ddt))) {
    out <- paste0(prefix, "_", names(by))
    quiet(do.call(run_ORFquant, c(list(for_ORFquant_file = paste0(prefix, "_for_ORFquant"),
                                       annotation_file = annotation, n_cores = 1, prefix = out,
                                       interactive = FALSE), by)))
    orfs <- read_tsv(result_files(out)$tsv)
    expect_true(all(orfs$gene_name == "DDT"), info = names(by))
    expect_true(all(golden$ORF_id_tr %in% orfs$ORF_id_tr), info = names(by))
  }
})

# GENCODE's CDS lines end before the stop codon, and its stop_codon lines add it to the CDS.
# Without them, as in GTFs with CDS lines from factR, the codon after the CDS is the stop codon
test_that("without stop_codon lines, the ORFs get the results of the example", {
  skip_on_os("windows")   # run_ORFquant() forks its workers
  quiet <- function(expr) invisible(utils::capture.output(suppressMessages(suppressWarnings(expr))))
  annotation <- annotation_of(gtf_example[line_type(gtf_example) != "stop_codon"], "no_stop_codons")
  expect_identical(ann_genc$stop_in_gtf, "*")
  expect_true(is.na(get(load(annotation))$stop_in_gtf))
  prefix <- file.path(work, "no_stop_codons", "example")
  quiet(prepare_for_ORFquant(annotation_file = annotation,
                             bam_file = example_file("chr22_example.bam"),
                             path_to_rl_cutoff_file = example_file("chr22_example_cutoffs.tsv"),
                             dest_name = prefix, n_cores = 1))
  msgs <- character()
  quiet(withCallingHandlers(run_ORFquant(for_ORFquant_file = paste0(prefix, "_for_ORFquant"), annotation_file = annotation,
                                         n_cores = 1, prefix = prefix, interactive = FALSE),
                            message = function(m) msgs <<- c(msgs, conditionMessage(m))))
  expect_true(any(grepl("end before their stop codon", msgs)))
  expect_equal(read_tsv(result_files(prefix)$tsv), read_tsv(golden_files()$tsv), tolerance = 1e-6)
})

# The lines with an intron just before the stop codon of transcript `tx`, which is in its last exon: the
# exon ends before the stop codon, and a new exon goes from the next stop codon of the genome (at least
# 3 nt downstream) to the end of the old exon. The ORF does not change, and its stop codon is the first
# codon of an exon. The UTR lines of `tx` are removed (GenomicFeatures does not read them)
intron_before_stop <- function(lines, tx) {
  of_tx <- which(line_attr(lines, "transcript_id") == tx)
  f <- strsplit(lines[of_tx], "\t", fixed = TRUE)
  type <- sapply(f, `[`, 3)
  from <- as.integer(sapply(f, `[`, 4))
  to <- as.integer(sapply(f, `[`, 5))
  plus <- f[[1]][7] == "+"
  stop <- which(type == "stop_codon")
  ex <- which(type == "exon" & from <= from[stop] & to >= to[stop])
  # the stop codons in the exon's sequence (5' to 3'), after the annotated one
  exon_seq <- as.character(Biostrings::getSeq(Rsamtools::FaFile(fasta),
                                              GRanges(f[[1]][1], IRanges(from[ex], to[ex]), f[[1]][7])))
  at <- if (plus) from[stop] - from[ex] + 1 else to[ex] - to[stop] + 1
  stops <- gregexpr("(?=TAA|TAG|TGA)", exon_seq, perl = TRUE)[[1]]
  new <- min(stops[stops >= at + 3])
  new <- if (plus) from[ex] + new - 1 else to[ex] - new - 1
  set_range <- function(line, a, b) {
    x <- strsplit(line, "\t", fixed = TRUE)[[1]]
    x[4:5] <- c(a, b)
    paste(x, collapse = "\t")
  }
  exon_line <- lines[of_tx[ex]]
  if (plus) {
    lines[of_tx[ex]] <- set_range(exon_line, from[ex], from[stop] - 1)
    new_exon <- set_range(exon_line, new, to[ex])
  } else {
    lines[of_tx[ex]] <- set_range(exon_line, to[stop] + 1, to[ex])
    new_exon <- set_range(exon_line, from[ex], new + 2)
  }
  lines[of_tx[stop]] <- set_range(lines[of_tx[stop]], new, new + 2)
  lines <- append(lines, sub('exon_id "([^"]*)"', 'exon_id "\\1_2"', new_exon), after = of_tx[ex])
  lines[!(line_type(lines) == "UTR" & line_attr(lines, "transcript_id") %in% tx)]
}

# ORF_category_Gen compares the last nucleotide before the stop codon, also when an intron is before the stop codon
test_that("an ORF with the annotated stop codon after an intron has the annotated stop", {
  skip_on_os("windows")   # run_ORFquant() forks its workers
  quiet <- function(expr) invisible(utils::capture.output(suppressMessages(suppressWarnings(expr))))
  # TOP3B (- strand) and TCN2 (+ strand): the transcripts of their ORFs and of their ref_id
  txs <- c("ENST00000357179.10", "ENST00000698268.1", "ENST00000698271.1")
  lines <- gtf_example
  for (tx in txs) lines <- intron_before_stop(lines, tx)
  annotation <- annotation_of(lines, "intron_before_stop")
  ann <- get(load(annotation))
  # the CDS has one more exon, with the stop codon only
  expect_identical(lengths(ann$cds_txs[txs]), lengths(ann_genc$cds_txs[txs]) + 1L)
  expect_identical(sum(width(ann$cds_txs[txs])), sum(width(ann_genc$cds_txs[txs])))
  expect_identical(unname(sapply(width(ann$cds_txs[txs]), function(w) w[length(w)])), rep(3L, 3))
  prefix <- file.path(work, "intron_before_stop", "example")
  quiet(prepare_for_ORFquant(annotation_file = annotation,
                             bam_file = example_file("chr22_example.bam"),
                             path_to_rl_cutoff_file = example_file("chr22_example_cutoffs.tsv"),
                             dest_name = prefix, n_cores = 1))
  quiet(run_ORFquant(for_ORFquant_file = paste0(prefix, "_for_ORFquant"), annotation_file = annotation,
                     n_cores = 1, prefix = prefix, interactive = FALSE, gene_name = c("TOP3B", "TCN2")))
  golden <- read_tsv(golden_files()$tsv)
  golden <- golden[golden$transcript_id %in% txs, ]
  expect_identical(golden$ORF_category_Gen, c("exact_start_stop", "exact_start_stop"))
  orfs <- read_tsv(result_files(prefix)$tsv)
  orfs <- orfs[match(golden$ORF_id_tr, orfs$ORF_id_tr), ]
  expect_identical(orfs$ref_id, golden$ref_id)
  expect_identical(orfs$ORF_category_Tx, golden$ORF_category_Tx)
  expect_identical(orfs$ORF_category_Gen, golden$ORF_category_Gen)
})
