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

## Performance

- `run_ORFquant()` is faster. It translated transcripts and ORFs with
  `translate(if.fuzzy.codon = "solve")`, which spends 50 to 80 ms per call
  building a table of the codons with ambiguous bases, such as N. It now
  skips that for sequences with only A, C, G and T, which give the same
  protein without it. With one core, a run on chr21 of a human sample took
  659 s instead of 717 s (8% less); with 8 cores, the expected saving (about
  7 s) was smaller than the variation between runs. Results don't change: on
  our test data and on chr21, all outputs are the same as before.

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
  made. Before, the vignette was not installed with the package, and it
  downloaded its data from links that no longer work.
