# ORFquant 1.03.0

## Installation and dependencies

- ORFquant now installs and runs on current Bioconductor, and needs
  Bioconductor 3.22 or later (R 4.5 or later). Before, on Bioconductor 3.22,
  `prepare_annotation_files()` and `run_ORFquant()` failed, because
  GenomicFeatures' `makeTxDbFromGFF()` and `makeTxDb()` are defunct; they now
  come from txdbmaker. `prepare_for_ORFquant()` no longer warns that
  `cigarRangesAlongQuerySpace()` and `cigarRangesAlongReferenceSpace()` are
  deprecated: it uses cigarillo's functions instead.
- devtools and GenomicFiles are no longer needed.
  `prepare_annotation_files(forge_BSgenome = TRUE)` installs the forged
  BSgenome package with `install.packages()`.
- doMC and foreach are no longer needed: `run_ORFquant()` runs its worker
  processes with the parallel package, which comes with R. So with more than
  one core, `run_ORFquant()` no longer registers doMC as the backend of
  foreach, and no longer attaches GenomicRanges when ORFquant itself isn't
  attached.
- `library(ORFquant)` now attaches only GenomicRanges and the packages it
  attaches (IRanges, S4Vectors, Seqinfo, BiocGenerics, generics and stats4).
  Before, it attached the 18 packages of its Depends field and theirs.
  Scripts that relied on it to attach others, for example Biostrings,
  rtracklayer or ggplot2, need their own `library()` calls.
- ORFquant's functions now work when called as `ORFquant::f()` without
  `library(ORFquant)`, and when other attached packages mask names they use.
  Before, `prepare_annotation_files()` failed with `could not find function
  "scanFaIndex"` in the first case, and with `type 'S4' passed to shift()`
  with data.table attached after ORFquant.
- knitr and rmarkdown are no longer required. They, Gviz, lemon, dplyr and
  GenomeInfoDb are now suggested packages: `create_ORFquant_html_report()`
  and `plot_orfquant_locus()` stop with a message when the ones they use are
  missing. `plot_orfquant_locus()` no longer fails with `could not find
  function "keepSeqlevels"` when GenomeInfoDb isn't attached, or partway
  through the plot when lemon isn't installed, and no longer attaches Gviz.
- Loading ORFquant no longer warns "replacing previous import".
- The repository no longer holds the source tarballs of earlier versions
  (`ORFquant_0.99.0.tar.gz` to `ORFquant_1.02.0.tar.gz`, and
  `ORFquant_manuscript_version.tar.gz`); they are still in its history, for
  example at commit `01b2da8`. `ORFquant_1.02.0.tar.gz` held older code than
  the repository's version 1.02.0, which still called `disjointExons()`, no
  longer in GenomicFeatures.

None of these changes affect results. On our test data, all outputs are the
same as before them (run with txdbmaker attached), apart from the two
deprecation warnings.

## New features

- `run_ORFquant()` also writes `<prefix>_Detected_ORFs.tsv`, a tab-separated
  table with one row per ORF: the columns of `ORFs_tx`, starting with its
  transcript coordinates, so the results can be read without R. Columns with
  several values per ORF have them separated by commas, and ranges are written
  as `seqname:start-end:strand`. `write_TSV_file = FALSE` turns it off. The
  other output files don't change.
- The script `run_orfquant.R`, installed with the package
  (`system.file("scripts", "run_orfquant.R", package = "ORFquant")`), runs
  `prepare_annotation_files()`, `prepare_for_ORFquant()` and `run_ORFquant()`
  from the command line, each in a new R process, with their default
  parameters, for example `Rscript run_orfquant.R --gtf genes.gtf --fasta
  genome.fa --bam sample.bam --offsets cutoffs.tsv --outdir results`, and
  lists the files it writes. `--annotation` reuses the annotation of an
  earlier run; `--cores`, `--gene-names` and `--gene-ids` are passed to
  `run_ORFquant()`; `--help` lists all options.
- The script `run_orfquant.sbatch`, installed next to `run_orfquant.R`, runs
  it as a SLURM job: `sbatch run_orfquant.sbatch` followed by
  `run_orfquant.R`'s options, with `--cores` set to the job's CPUs. It asks
  for 16 CPUs, 64 GB of memory and 24 hours, and its comments give the time
  and memory measured on human samples.
- `prepare_annotation_files()`, `prepare_for_ORFquant()` and `run_ORFquant()`
  now check their input before the long steps, and stop with a message that
  says what is wrong. Before, some wrong input gave an error from deep inside
  a package, some after hours, and some no error at all. The runs on correct
  input are not affected. The checks are:
  - Files that don't exist (the genome sequence too, and all the inputs of
    `prepare_for_ORFquant()`) are listed together, before anything is written.
  - The directory of the output files (`dest_name`, `prefix`) must exist and
    be writable. Before, `prepare_for_ORFquant()` and `run_ORFquant()` failed
    only when they saved their results, after the BAM file was read or the
    ORFs were quantified.
  - `n_cores` must be a number, 1 or more, in `run_ORFquant()` too. Before,
    a missing or zero `n_cores` failed late, with `argument "n_cores" is
    missing` or `object 'ORFs_found' not found`.
  - The GTF file must have exon lines with `transcript_id` and `gene_id`
    (before, a GTF without `gene_id` or `transcript_id` failed in
    GenomicFeatures, and one without exon lines gave an annotation with 68 of
    the 140 transcripts of the example, rebuilt from its CDS and UTR lines).
    The error for chromosomes of the GTF that the genome sequence doesn't
    have now says so, and gives the chromosome names of the genome.
  - The cutoff table must have rows, numbers in `read_length` and `cutoff`,
    and at least one `compartment` that is `nucl` or a chromosome of the
    annotation. A compartment that is neither never gave P-sites, so one
    among others now gives a warning. Before, a table with no usable row gave a
    `for_ORFquant` file without P-sites, and no error.
  - The BAM file must have a chromosome that is in the annotation, and
    `prepare_for_ORFquant()` stops when it found no P-sites. Before, a BAM
    file with the chromosome names `1`, `2` for an annotation with `chr1`,
    `chr2` gave an empty `for_ORFquant` file without an error.
  - The P-sites in `for_ORFquant_file` must be on the same chromosomes, with
    the same lengths, as the annotation of `run_ORFquant()`, which did not
    check this. It catches a `for_ORFquant_file` made for another genome
    assembly or annotation.
- The progress messages of the three steps and of the plots, such as
  "Calculating P-sites positions and junctions ... ", now come from
  `message()`, not from `cat()`. So they go to the standard error, not to the
  standard output, and `suppressMessages()` or `--quiet`-style redirects can
  hide them. Their text is the same. The scripts `run_orfquant.R` and
  `run_orfquant.sbatch` show them as before.

## Changes in results

- `run_ORFquant()` now applies `canonical_start_only = TRUE` (the default)
  with `n_cores = 1` too; before, only runs with `n_cores > 1` did. Serial
  runs read the genetic code's alternative initiation codons (TTG and CTG in
  the standard code) as M when they were the first codon of a translated
  sequence. This affected two things:
  - A TTG or CTG starting at one of the first three bases of a transcript
    was a start codon. This changed `longest_ORF`, `compatible_with_longest`,
    `compatible_ORF_id_tr_longest`, `pct_fr_st` and `ave_pct_fr_st`.
  - `NC_protein_isoform` compares an ORF's last 10 amino acids with those of
    the annotated CDS. A TTG or CTG at the start of that stretch of the CDS
    was read as M, giving `C` instead of `same` (or `N_C` instead of `N`).

  On chr21 of a human sample this changed the annotation of 8 of 87 ORFs,
  but not which ORFs were found, their coordinates, P-sites or p-values.
  Serial runs now give the same results as parallel runs, which are
  unchanged.
- With `unique_reads_only = TRUE`, `run_ORFquant()` now skips genomic regions
  with too little signal from uniquely mapping reads. A region needs more
  than 4 P-site entries (ranges in `for_ORFquant`) to be analysed. The check
  on the unique P-sites was computed but its result was dropped, so only the
  check on all P-sites applied. In our test data and on chr21 of a human
  sample (30 such regions), no ORF was found in these regions before either:
  the ORFs are unchanged, and only the transcripts selected there (2 on
  chr21) are no longer in `selected_txs` and `ORFs_txs_feats`. The default,
  `unique_reads_only = FALSE`, is unaffected.
- `prepare_for_ORFquant()` now stores the bigWigs of uniquely mapping P-sites
  with mismatches (`path_to_P_sites_uniq_mm_plus_bw` and
  `path_to_P_sites_uniq_mm_minus_bw`) in `P_sites_uniq_mm`. Before, it stored
  them in `P_sites_uniq`, replacing the unique P-sites, and left
  `P_sites_uniq_mm` empty. ORFs called from such a file had `P_sites_raw_uniq`
  and `pval_uniq` computed on the mismatch track and `P_sites_raw_uniq_mm` 0,
  and `ORFs_feat` and `ORFs_txs_feats` had wrong `unique_reads`. With
  `unique_reads_only = TRUE`, ORFs were also detected and quantified on the
  mismatch track: on chr21 of a human sample, 5 ORFs were found instead of 76.
  With the default, `unique_reads_only = FALSE`, the same ORFs were found.
  Files made without these two bigWigs, from a BAM and a cutoff table, or by
  RiboseQC are not affected. Re-run `prepare_for_ORFquant()` for files made
  with them.
- With `canonical_start_only = FALSE`, `run_ORFquant()` no longer reads a TTG
  or CTG as M when it starts the stretch of the annotated CDS that
  `NC_protein_isoform` compares with the ORF's last 10 amino acids. These
  ORFs were labelled `C` instead of `same` (or `N_C` instead of `N`); on
  chr21 of a human sample, 5 of 87 ORFs change from `C` to `same`. Nothing
  else changes, and the default, `canonical_start_only = TRUE`, is
  unaffected.
- `calc_orf_pval()` now uses `tapers` for the degrees of freedom of the
  multitaper F-test (`2 * tapers - 2`). Before, it always used 46, the value
  for the default 24 tapers. `run_ORFquant()` always uses the default, so its
  results are unchanged; only direct calls to `calc_orf_pval()` with another
  `tapers` get different `pval` and `pval_uniq`.
- `prepare_annotation_files()` now gives the biotype `protein_coding` to
  transcripts without a biotype that have CDS lines, and to their genes.
  Before, they were `no_type`, as are all transcripts and genes of GTFs
  without biotypes, such as UCSC's RefSeq GTF (`hg38.ncbiRefSeq.gtf`).
  `run_ORFquant()` then prefers these transcripts as an ORF's compatible
  transcript, taking the first in sorted order, as with GENCODE's
  `protein_coding` transcripts (before, the first one found), and
  `plot_ORFquant_results()` labels their ORFs `protein_coding`, not
  "non-coding RNA". On chr21 of a human sample with UCSC's RefSeq GTF, all 89
  ORFs change their `gene_biotype`, `transcript_biotype` and
  `compatible_biotype` from `no_type` to `protein_coding`, and 6 their
  compatible transcript (`compatible_tx`, `compatible_ORF_id_tr` and
  `compatible_ORF_id_tr_longest`). Which ORFs are found, their categories,
  P-sites and p-values don't change. Without biotypes, transcripts with a CDS
  that GENCODE calls `nonsense_mediated_decay` are `protein_coding` too, so an
  ORF on one can be its own compatible transcript and get another
  `ORF_category_Tx_compatible`: on GENCODE's chr22 transcripts with their
  biotypes removed, 1 of 16 ORFs gets `overl_dORF` instead of GENCODE's
  `N_truncation`. GTFs in which every transcript has a biotype, such as
  GENCODE, Ensembl and NCBI RefSeq GTFs, are not affected.
- `run_orfquant.sbatch` runs R with one BLAS thread per process
  (`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`), and the README asks
  users of a multi-threaded BLAS to do the same (see Performance). The
  p-values of the longest ORFs depend on the number of BLAS threads: with
  more than one, `multitaper::dpss()` gives slightly different tapers for
  ORFs of about 20,000 nt or more. So with a multi-threaded BLAS, such as
  OpenBLAS in conda's R, `pval` and `pval_uniq` of these ORFs depended on the
  number of CPUs of the machine or job; with one thread they don't. On a
  human sample with GENCODE 47, 9 of 33,146 ORFs (in TTN, OBSCN, NEB and
  others) change, by at most 1e-11 of their value; nothing else changes.
  With R's own BLAS, which uses one thread, nothing changes.

## Bug fixes

- `run_ORFquant()` no longer fails with `NA/NaN argument` when exactly one
  genomic region has ORFs, for example in a run on one gene
  (lcalviell/ORFquant#26).
- `run_ORFquant()` no longer drops the names of `ORFs_tx` (the ORF ids) when
  every region with ORFs (or every one in a block of 1000 such regions) has
  the same number of ORFs, 2 or more, for example in a run on one gene with
  2 ORFs. The `ORF_id_tr` column was not affected.
- `prepare_annotation_files()` no longer fails with `NA/NaN argument` when
  the annotation has exactly one protein-coding transcript.
- `prepare_annotation_files()` now reads the GTF's ids, biotypes and gene
  names by attribute name, one row per transcript, taking each value from
  any line of the transcript (for genes, also from a line of the gene).
  Before, it kept one row per distinct set of values and named the columns
  it found by their position:
  - GTFs with biotypes on only some lines or transcripts failed with
    `subscript contains invalid names`, for example NCBI RefSeq GTFs (only
    their gene lines have `gene_biotype`) and GTFs written by
    `gffread -F -T` (biotypes on transcript lines only);
  - GTFs with both `gene_type` and `gene_biotype` (or both
    `transcript_type` and `transcript_biotype`) gave an annotation with its
    columns mixed up, for example gene names as transcript ids, without an
    error.

  Missing biotypes are now `no_type`, as all were when no line had one (or
  `protein_coding`, see "Changes in results"). The
  transcript biotype `mRNA`, NCBI's name for coding transcripts, is read as
  `protein_coding`. `run_ORFquant()` then prefers these transcripts as an
  ORF's compatible transcript (`compatible_tx`), and `plot_ORFquant_results()`
  no longer labels their ORFs "non-coding isoform", as with GENCODE. Gene
  names are also read from `gene` (NCBI) and `ref_gene_name` (StringTie)
  values, so `run_ORFquant(gene_name = ...)` works with these GTFs.
  `?prepare_annotation_files` lists the attributes read. GTFs that worked
  before give the same annotation, apart from these gene names, the `mRNA`
  biotype, and biotypes that were missing (`NA`) for some transcripts, now
  `no_type` or `protein_coding`.
- `prepare_annotation_files()` now stops with an error saying that the GTF
  file has no CDS lines when it has none, for example StringTie's output.
  Before, it failed with `wrong sign in 'by' argument`.
- `prepare_annotation_files()` no longer fails with `<n> elements in value to
  replace <m> elements` when the GTF has genes with exons on both strands or
  on more than one chromosome, for example UCSC's RefSeq GTF
  (`hg38.ncbiRefSeq.gtf`), whose genes in the pseudoautosomal regions have
  the same id on chrX and chrY. The annotation's `genes` leaves these genes
  out, as before, and their UTR, intron and non-coding exon regions get no
  gene id. With such genes, `run_ORFquant()` with `gene_name` or `gene_id`,
  `plot_ORFquant_results()` and `plot_orfquant_locus()` no longer fail with
  `subscript contains invalid names`. Other GTFs give the same annotation
  and results.
- `run_ORFquant()` no longer fails with `missing value where TRUE/FALSE
  needed` (with more than one core, `task 1 failed - "missing value where
  TRUE/FALSE needed"`) when the transcript of an ORF has no biotype in the
  annotation (lcalviell/ORFquant#25), as in annotations made by RiboseQC or
  earlier versions from GTFs where some transcripts have none. These
  transcripts are taken as not protein-coding.
- `prepare_for_ORFquant()` now reads the cutoff table's columns by name. It
  needs the columns `read_length`, `cutoff` and `compartment`, in any order,
  and ignores any others. Before, it checked only that the table had 3
  columns and took the compartment from the third. A table headed `rl` (as
  the error message said) gave no P-sites and no error, and so did one with
  `compartment` first. Tables with the documented header give the same
  results as before. Some headers with other names, for example
  `read_lengths` or a third column `comp`, worked before and now stop with an
  error. The error message now names `read_length`, not `rl`.
- `prepare_for_ORFquant()` now stops with an error when only one bigWig of a
  pair is given, for example `path_to_P_sites_plus_bw` without
  `path_to_P_sites_minus_bw` (likewise for the `uniq` and `uniq_mm` pairs), or
  when a cutoff table is given with `path_to_P_sites_plus_bw` or
  `path_to_P_sites_minus_bw`. Before, these calls ran without a warning. The
  track then held only the plus strand when only the plus bigWig was given,
  and nothing when only the minus bigWig was given; with a cutoff table, it
  replaced the P-sites computed from the BAM. Other calls are not affected.
- `plot_ORFquant_results()` no longer fails with `'breaks' are not unique`
  when the largest number of selected transcripts per gene is 3, 6 or 9
  (lcalviell/ORFquant#18). All other plots are unchanged.
- `plot_orfquant_locus()`:
  - no longer fails for genes none of whose discarded transcripts is coding,
    or without junction reads;
  - no longer fails with `Too many stacks to draw` for genes with many
    transcripts: a track with more than 21 transcripts gets more height, and
    the page with it. Plots with up to 21 transcripts in every track keep
    their 7-inch page;
  - closes its PDF when plotting fails. Before, the device stayed open, and
    the session's later plots went into that file;
  - labels each discarded ORF with its own transcript. Before, in most genes,
    many labels named another discarded transcript (with the right
    coordinates);
  - colours the selected ORFs by `ORFs_pM`, from dark green (lowest) to
    bright green (highest), as in the legend. Before, they were all white.

  In our test data, 17 of 63 genes failed; all 63 plot now.
- `create_ORFquant_html_report()` now renders the report from a copy of the
  template in a temporary directory, not in the installed package. It now
  works when the R library is read-only (containers, shared conda
  environments), and reports rendered at the same time from one installation
  no longer share intermediate files, which could make a report list another
  report's input or fail.
- `create_ORFquant_html_report()` now ends its `sink()` when rendering fails.
  Before, the rest of the session's console output went to
  `<output>_ORFquant_report_output.txt`.
- The `getSeq()` error for ranges that wrap twice around a circular
  chromosome now reads "Ranges wrapping twice isn't implemented yet...",
  without the stray quote, line break and spaces it had.
- `run_ORFquant()` and `prepare_annotation_files()` now close the TxDb
  databases they make when they no longer need them. Their logs no longer
  show `Error in x$.self$finalize() : attempt to apply non-function` (often
  hundreds of lines in a whole-genome run) and `call dbDisconnect() when
  finished working with a connection` (lcalviell/ORFquant#3). These messages
  came from R closing the databases later, and did not change the results.

## Performance

- On a human sample with a 1.8 GB BAM and GENCODE 47, the changes below
  together make a run take 3 h 11 min with 16 cores instead of 7 h 46 min
  (59% less), and at most 35 GB of memory instead of 45 GB. With 64 cores, a
  run takes 57 min and at most 115 GB. Each entry below gives what its change
  saved when it was made.
- `run_ORFquant()` is faster. On a human sample with a 1.8 GB BAM and
  GENCODE 47, with 16 cores, a run takes 4 h 3 min instead of 7 h 46 min
  (48% less) and at most 37 GB of memory instead of 45 GB. With one core, a
  run on chr21 of a human sample takes 431 s instead of 704 s (39% less),
  and with 8 cores 90 s instead of 145 s; these chr21 times don't include
  the fourth change below, or `annotate_splicing()` setting its exon columns
  once per ORF. Results don't change: on our test data and on chr21, all
  outputs are the same as before, and on the whole human sample the
  exported files are. Six changes make it faster:
  - It translated transcripts and ORFs with
    `translate(if.fuzzy.codon = "solve")`, which spends 50 to 80 ms per
    call building a table of the codons with ambiguous bases, such as N. It
    now skips that for sequences with only A, C, G and T, which give the
    same protein without it.
  - It set columns of GRanges objects many times per ORF or per exon, for
    example each ORF's p-values in `calc_orf_pval()` and each exon's splice
    type in `annotate_splicing()`, and each assignment takes 6 to 37 ms. It
    now keeps the values in plain vectors and sets each column once.
  - `annotate_ORFs()` mapped each ORF to the annotated transcripts one
    transcript at a time; it now maps it to all of them in one call.
    `annotate_splicing()` looked for the annotated CDS exons overlapping
    each exon of an ORF, set the exon's columns, and added the exon to its
    result with `c()` and `sort()`, one exon at a time; it now does each
    once per ORF.
  - Each genomic region looked for its P-sites and junctions among those of
    the whole genome, which takes 12 to 32 ms per region for each of the
    three P-site tracks and the junctions on a human sample.
    `run_ORFquant()` now finds them for all regions at once and gives each
    region only its own. Likewise, `detect_translated_orfs()` looked up
    each transcript's exons and introns by name among all annotated
    transcripts, about 20 ms per lookup with a human annotation; it now
    looks them up among the region's selected transcripts.
  - The last step, which combines the results of all genomic regions,
    combined them with `GRangesList()`, which handles the columns and the
    sequence information (Seqinfo) of its elements one element at a time.
    It now gathers each column from all elements and binds it once, and
    merges the Seqinfo pairwise. On the results of a whole human genome,
    this step takes 108 s instead of 490 s.
  - Setting a column of a GRanges object with `x$name <- value` takes about
    20 ms even when it is done once, because it also checks and updates the
    whole object. Many steps set several columns of the same object in a
    row; they now take the object's columns, set them, and put them back in
    one step. On chr21, a run made about 13,000 such assignments and now
    makes about 900. The times above don't include this change; on chr21 it
    makes a run 10% faster with one core and 6% faster with 8 cores.
- `run_orfquant.sbatch` is faster with many CPUs, because it runs R with one
  BLAS thread per process (see Changes in results). With a multi-threaded
  BLAS, such as OpenBLAS in conda's R, each process that `run_ORFquant()`
  forks started as many BLAS threads as the job had CPUs, and these threads
  took CPU time from the other processes. On the human sample above, with 64
  CPUs, a run takes 1 h 11 min instead of 1 h 27 min, and half the CPU time
  (53 h instead of 102 h). Scripts that call `run_ORFquant()` or
  `run_orfquant.R` get the same by setting `OPENBLAS_NUM_THREADS=1` and
  `OMP_NUM_THREADS=1` before they start R (see the README).
- With more than one core, `run_ORFquant()` gives the genomic regions to its
  worker processes one at a time, in order, each to the next worker that is
  free. Before, each worker got an equal share of the regions at the start.
  The time per region varies a lot, so at the end most workers were idle
  while the last ones finished their shares. On the human sample above, with
  64 cores and one BLAS thread per process, the regions take 53 min instead
  of 56 min, and the whole run 68 min instead of 72 min. With 16 cores, the
  run takes the same time as before (3 h 20 min) and at most 34 GB of
  memory instead of 37 GB. On chr21 with 8 cores, the regions take 75 s
  instead of 86 s. Results don't change. The
  "% completed" lines now show the progress of the whole run, not of one
  worker's share. When the regions start, each worker prints a line
  "starting worker pid=...". If a region fails, the error message starts
  with "one node produced an error" instead of "task N failed".
- `prepare_for_ORFquant()` has a new argument `n_cores` (default 1), and
  `run_orfquant.R` gives it the value of `--cores`. With more than one core,
  forked processes calculate the P-sites and junctions of the chunks of
  `chunk_size` alignments while the next chunks are read. With any number of
  cores, the P-sites of all chunks are now added up once, after the last
  chunk. Before, the P-sites of each chunk were added to those of all the
  earlier chunks, which took longer with each chunk. On the human sample
  above (68 million alignments), with 16 cores, "Calculating P-sites
  positions and junctions" takes 108 s instead of 8 min 26 s, and at most
  19 GB of memory instead of about 5 GB. With one core, on a test BAM file
  read in chunks of 1,000 alignments, it takes 78 s instead of 139 s.
  Results don't change: the `for_ORFquant` file is the same as before, with
  any `chunk_size` and number of cores. One exception, which doesn't change
  the results of `run_ORFquant()`: in a BAM file with more than `chunk_size`
  alignments, a P-site track without P-sites (for example `P_sites_uniq_mm`
  for a BAM file without MD tags) no longer has a `score` column, the same
  as in a BAM file with fewer alignments.
- With more than one core, the export at the end of `run_ORFquant()` is
  faster. Forked processes save `<prefix>_tmp_ORFquant_results` and build
  the result tables, at most 3 tables at a time. Then a forked process saves
  `<prefix>_final_ORFquant_results` while the TSV, FASTA and GTF files are
  written. On the human sample above, with 16 cores, "Exporting ORFquant
  results" takes 143 s instead of 6 min 13 s, and uses about as much memory
  as the genomic regions before it, instead of 10 GB. The warnings of the
  forked processes are shown as before. Results don't change. With this
  change and the one above, a run on the human sample with 16 cores takes
  3 h 11 min instead of 3 h 20 min, and at most 35 GB of memory instead of
  34 GB. With 64 cores, a run takes 57 min instead of 68 min, and at most
  115 GB of memory, as before.

## Documentation

- The new help page `?ORFquant_output` (also `?ORFs_tx`) describes each
  column of `ORFs_tx` and of `<prefix>_Detected_ORFs.tsv`, how ORFquant
  computes it, and the values of `ORF_category_Tx`,
  `ORF_category_Tx_compatible` and `ORF_category_Gen`.
- Help pages that contradicted the code are corrected:
  - `?select_quantify_ORFs` gave the defaults of `cutoff_cums`, `cutoff_pct`
    and `cutoff_P_sites` as 99, 1 and 10; they are `NA` (not applied), 2 and
    `NA`. It also said that `P_sites` is `P_sites_raw` divided by the scaling
    factor, but it is multiplied by it.
  - `?annotate_ORFs` called `N_truncation` an N-terminal extension. It said
    that `compatible_with` holds transcript ids, but it holds ORF ids. Its
    return value was described as that of `annotate_splicing()`.
  - `?detect_translated_orfs` called `pct_fr` a percentage, but it is a
    fraction.
  - `?calc_orf_pval` said that `cutoff` applies to the average in-frame
    signal per codon, but it applies to `pct_fr`. Its return value was
    described as that of `select_start()`.
- The vignette, `vignette("ORFquant")`, now runs the whole analysis on example
  data included in the package (`inst/extdata`): 7 genes of human chr22, with
  their GENCODE 47 annotation, Ribo-seq reads from SRA run SRR15513199 and
  their P-site offsets, and the sequence of chr22, replaced by N away from
  these genes. `inst/scripts/make_example_data.R` describes how they were
  made, and `inst/extdata/README.md` gives their sources, terms of use and the
  papers to cite. Before, the vignette was not installed with the package, and
  it downloaded its data from links that no longer work.
- The repository no longer holds `ORFquant-manual.pdf`. It was the manual of
  version 1.02.0 from June 2020, without the changes to the help pages since
  then. Read the help pages with `help(package = "ORFquant")`, or make a PDF
  of them with `R CMD Rd2pdf` (this needs LaTeX). The old PDF is still in the
  repository's history, for example at commit `5fc1d86`.
