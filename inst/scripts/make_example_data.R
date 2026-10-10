## How the example data in inst/extdata were made, for the vignette: 7 genes of
## human chr22 (GENCODE 47) with Ribo-seq reads from SRA run SRR15513199.
##
## The GTF and BAM are py-orfquant's chr22_uorf_region test fixture
## (scripts/build_uorf_fixtures.py), copied unchanged:
## - GTF: all lines of gencode.v47.annotation.gtf for the chr22 genes TXNRD2,
##   RANBP1, PPIL2, TOP3B, DDT, TCN2 and PISD (gzipped here).
## - BAM: the primary alignments within 500 bp of these genes, from the
##   alignment of SRR15513199 to GRCh38 with STAR 2.7.10a and the GENCODE 47
##   annotation; only the MD and NH tags kept, with a chr22-only header.
## The cutoff table has the P-site offsets that RiboseQC used for the P-sites
## of this sample in the gencode_pilot pipeline (RiboseQC_analysis() with
## rescue_all_rls = TRUE): the offsets of the read lengths it selected, and 12
## for all other read lengths of the sample. With this table,
## prepare_for_ORFquant() gives the P-sites of RiboseQC's _for_ORFquant file.
## RiboseQC's _P_sites_calcs file lists other offsets, which it did not use.
## The genome is chr22 of GENCODE's GRCh38.primary_assembly.genome.fa, with
## every base more than 1000 bp from these genes replaced by N, so that the
## bgzipped FASTA is small but keeps chr22's length and coordinates.
##
## Run from the package's root directory. The input paths are those of the
## machine the files were made on.

fixtures <- "/fast/home/k/kgresov/py-orfquant/tests/data/fixtures"
genome_fa <- "/fast/AG_Ohler/gabriel/gencode_pilot/resources/GRCh38.primary_assembly.genome.fa"
riboseqc_all <- "/fast/AG_Ohler/gabriel/gencode_pilot/results/riboseqc/data/SRR15513199.Aligned.sortedByCoord.out/_results_RiboseQC_all"
out <- "inst/extdata"
margin <- 1000

dir.create(out, recursive = TRUE, showWarnings = FALSE)
gtf <- file.path(fixtures, "chr22_uorf_region.gtf")
con <- gzfile(file.path(out, "chr22_example.gtf.gz"), "w")
writeLines(readLines(gtf), con)
close(con)
stopifnot(file.copy(file.path(fixtures, c("chr22_uorf_region.bam", "chr22_uorf_region.bam.bai")),
                    file.path(out, c("chr22_example.bam", "chr22_example.bam.bai")),
                    overwrite = TRUE))

# the read lengths of the sample (the columns of rld) and the offsets of the selected ones
riboseqc <- new.env()
load(riboseqc_all, envir = riboseqc)
read_lengths <- as.numeric(sub("^reads_", "", colnames(riboseqc$res_all$read_stats$rld)))
selected <- as.data.frame(riboseqc$res_all$selection_cutoffs$results_choice$nucl$final_choice)
cutoffs <- data.frame(read_length = read_lengths, cutoff = 12, compartment = "nucl")
cutoffs$cutoff[match(as.numeric(selected$read_length), read_lengths)] <- selected$cutoff
write.table(cutoffs, file.path(out, "chr22_example_cutoffs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

genes <- rtracklayer::import(gtf)
genes <- genes[genes$type == "gene"]
keep <- GenomicRanges::reduce(genes + margin, ignore.strand = TRUE)
fa_index <- Rsamtools::scanFaIndex(genome_fa)
chr22 <- fa_index[GenomicRanges::seqnames(fa_index) == "chr22"]
seq <- Rsamtools::scanFa(genome_fa, chr22)[[1]]
drop <- IRanges::gaps(IRanges::ranges(keep), start = 1, end = length(seq))
seq <- Biostrings::replaceAt(seq, drop, Biostrings::DNAStringSet(strrep("N", IRanges::width(drop))))

tmp <- tempfile(fileext = ".fa")
Biostrings::writeXStringSet(Biostrings::DNAStringSet(list(chr22 = seq)), tmp, width = 60)
Rsamtools::bgzip(tmp, file.path(out, "chr22_example.fa.gz"), overwrite = TRUE)
Rsamtools::indexFa(file.path(out, "chr22_example.fa.gz"))
unlink(tmp)
