## Run the three steps of ORFquant from the command line, each in a new R
## process: prepare_annotation_files(), prepare_for_ORFquant() and run_ORFquant().
## Run it with Rscript; --help lists the options. The installed copy is at
## system.file("scripts", "run_orfquant.R", package = "ORFquant").

usage <- "Usage: Rscript run_orfquant.R (--gtf FILE --fasta FILE | --annotation FILE)
                             --bam FILE --offsets FILE --outdir DIR [options]

Runs prepare_annotation_files(), prepare_for_ORFquant() and run_ORFquant(),
with their default parameters, each in a new R process, and writes all output
files to DIR.

  --gtf FILE         annotation in GTF format (may be gzipped)
  --fasta FILE       genome sequence in FASTA format (may be bgzipped), indexed
                     with samtools faidx; the annotation refers to it by its
                     path, so keep it there
  --annotation FILE  the _Rannot file of an earlier run (in its
                     DIR/annotation), instead of --gtf and --fasta
  --bam FILE         Ribo-seq reads aligned to the genome with a spliced
                     aligner (e.g. STAR)
  --offsets FILE     tab-separated table of P-site offsets, with the columns
                     read_length, cutoff and compartment
  --outdir DIR       output directory, created if it doesn't exist
  --sample NAME      prefix of the output files (default: the BAM file's name
                     without .bam)
  --cores N          number of cores for run_ORFquant() (default: 1)
  --gene-names LIST  analyse only the genomic regions with these genes, a
                     comma-separated list of gene names
  --gene-ids LIST    the same with gene ids; with --gene-names, only the
                     regions with genes from both lists
  -h, --help         print this help and exit

With --gtf and --fasta, the annotation is written to DIR/annotation; give its
_Rannot file to --annotation to reuse it for other samples. At the end, the
output files are listed, one per line: their kind, a tab and their path.
"

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat(usage)
  quit(status = 0)
}
if (length(args) == 0) {
  cat(usage)
  quit(status = 1)
}

opts <- c(gtf = NA, fasta = NA, annotation = NA, bam = NA, offsets = NA, outdir = NA,
          sample = NA, cores = "1", "gene-names" = NA, "gene-ids" = NA)
i <- 1
while (i <= length(args)) {
  key <- sub("=.*", "", sub("^--", "", args[i]))
  if (!startsWith(args[i], "--") || !key %in% names(opts)) {
    stop("unknown option ", args[i], "; see --help")
  }
  if (grepl("=", args[i])) {
    opts[[key]] <- sub("^[^=]*=", "", args[i])
  } else {
    if (i == length(args)) {
      stop("no value for ", args[i])
    }
    i <- i + 1
    opts[[key]] <- args[i]
  }
  i <- i + 1
}

from_gtf <- !is.na(opts[["gtf"]]) || !is.na(opts[["fasta"]])
if (from_gtf == !is.na(opts[["annotation"]])) {
  stop("give either --gtf and --fasta, or --annotation")
}
needed <- c(if (from_gtf) c("gtf", "fasta") else "annotation", "bam", "offsets", "outdir")
if (any(is.na(opts[needed]))) {
  stop("missing ", paste0("--", needed[is.na(opts[needed])], collapse = ", "))
}
for (f in setdiff(needed, "outdir")) {
  if (!file.exists(opts[[f]])) {
    stop("--", f, ": file not found: ", opts[[f]])
  }
  opts[[f]] <- normalizePath(opts[[f]])
}
if (!grepl("^[1-9][0-9]*$", opts[["cores"]])) {
  stop("--cores must be a whole number above 0, not ", opts[["cores"]])
}
gene_list <- function(x) {
  if (is.na(x)) NA else trimws(strsplit(x, ",")[[1]])
}
if (!nzchar(system.file(package = "ORFquant"))) {
  stop("ORFquant is not installed in this R")
}
## Each step runs in a new R process: the worker processes of run_ORFquant()
## would copy the memory left by earlier steps in the same process (with a
## human annotation and 16 cores, 100 GB instead of 25 GB). library() first:
## loaded by ORFquant::f(), its dependency SparseArray makes R warn "stack
## imbalance".
run_step <- function(call) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote("suppressPackageStartupMessages(library(ORFquant))"),
                      "-e", shQuote(deparse1(call))))
  if (status != 0) {
    stop(deparse(call[[1]]), "() failed", call. = FALSE)
  }
}

dir.create(opts[["outdir"]], recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(opts[["outdir"]])
sample <- opts[["sample"]]
if (is.na(sample)) {
  sample <- sub("\\.bam$", "", basename(opts[["bam"]]))
}
prefix <- file.path(outdir, sample)

if (from_gtf) {
  annotation_dir <- file.path(outdir, "annotation")
  run_step(bquote(ORFquant::prepare_annotation_files(annotation_directory = .(annotation_dir),
                                                     gtf_file = .(opts[["gtf"]]),
                                                     genome_seq = .(opts[["fasta"]]),
                                                     forge_BSgenome = FALSE)))
  annotation_file <- file.path(annotation_dir, paste0(basename(opts[["gtf"]]), "_Rannot"))
} else {
  annotation_file <- opts[["annotation"]]
}

run_step(bquote(ORFquant::prepare_for_ORFquant(annotation_file = .(annotation_file),
                                               bam_file = .(opts[["bam"]]),
                                               path_to_rl_cutoff_file = .(opts[["offsets"]]),
                                               dest_name = .(prefix))))
for_ORFquant_file <- paste(prefix, "for_ORFquant", sep = "_")

run_step(bquote(ORFquant::run_ORFquant(for_ORFquant_file = .(for_ORFquant_file),
                                       annotation_file = .(annotation_file),
                                       n_cores = .(as.integer(opts[["cores"]])),
                                       prefix = .(prefix),
                                       gene_name = .(gene_list(opts[["gene-names"]])),
                                       gene_id = .(gene_list(opts[["gene-ids"]])),
                                       interactive = FALSE)))

outputs <- c(annotation = annotation_file,
             for_ORFquant = for_ORFquant_file,
             results = paste(prefix, "final_ORFquant_results", sep = "_"),
             orfs_tsv = paste(prefix, "Detected_ORFs.tsv", sep = "_"),
             orfs_gtf = paste(prefix, "Detected_ORFs.gtf", sep = "_"),
             proteins_fasta = paste(prefix, "Protein_sequences.fasta", sep = "_"),
             region_results = paste(prefix, "tmp_ORFquant_results", sep = "_"))
outputs <- outputs[file.exists(outputs)]
cat("\nORFquant output files:\n")
cat(paste(names(outputs), outputs, sep = "\t"), sep = "\n")
