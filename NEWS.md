# ORFquant 1.4.0

ORFquant's own code does not change. The changes below come from RiboseQC
1.5.0, whose functions ORFquant uses for the first two steps and exports
again: `prepare_annotation_files()`, `load_annotation()`, `FaFile_Circ`,
`prepare_for_ORFquant()`, `get_ps_fromspliceplus()` and
`get_ps_fromsplicemin()`. The NEWS of RiboseQC 1.5.0 gives all its changes.

## Installation and dependencies

- ORFquant now needs RiboseQC 1.5.0 or later (`Remotes` installs the tag
  `v1.5.0`). RiboseQC 1.5.0 no longer has the methods of `unlist()` for
  `GRanges` and `GAlignments`, which changed `unlist()` in the whole R
  session. ORFquant does not need them since 1.3.2. RiboseQC 1.5.0 needs
  *ggplot2* 3.4.0 or newer.

## New arguments

- `prepare_for_ORFquant()` has the new argument `stranded`, with the values
  of `RiboseQC_analysis()`: `TRUE` (the default), `FALSE` or `"inverse"`.
  Use `"inverse"` for libraries with the reads on the strand opposite to the
  RNA: the function then inverts the strand of each read before it
  calculates the P-sites and the junction reads. Before, for these libraries
  each P-site was on the other strand, and the reads on the junctions were
  not counted. `TRUE` and `FALSE` give the same results as before.
- `prepare_annotation_files()` has the new argument `genetic_codes`: the
  genetic codes of some chromosomes, e.g. `genetic_codes = c(chrM = "3")`,
  with the ids of `Biostrings::GENETIC_CODE_TABLE`. They replace the codes
  that the chromosome names give. So the mitochondrion of a yeast or
  Drosophila genome with UCSC names (`chrM`), which gets code 2 from its
  name, can get its own code, and `run_ORFquant()` finds and translates its
  ORFs with that code. The default, `NULL`, gives the same codes as before.

## Input checks

- `prepare_annotation_files()` warns when the lines of the GTF file do not
  agree: when the lines of a transcript give different values of an
  attribute (for example two biotypes), or the transcripts of a gene give
  different gene names or biotypes. Then the order of the lines chooses the
  value (lcalviell/ORFquant#10). The annotation does not change.
- `prepare_annotation_files()` stops at once for a GTF file without CDS
  lines also when the file has `start_codon` or `stop_codon` lines. Before,
  it made the annotation, and `run_ORFquant()` stopped later with `subscript
  contains out-of-bounds indices`.
- `prepare_annotation_files()` stops with an error that names a transcript
  when a CDS ends after the end of its chromosome (one that is not circular)
  in the genome sequence, for example when the GTF file and the FASTA file
  are of different genome assemblies.
- `load_annotation()`, and so `prepare_for_ORFquant()` and `run_ORFquant()`,
  stop with an error that says what is wrong when the annotation file is not
  a file of `prepare_annotation_files()`, for example the GTF file. Before,
  the error was `bad restore file magic number`.
- `prepare_for_ORFquant()` stops at once when a compartment of the cutoff
  file is a chromosome that is not circular, for example `chr22`. Before, a
  read length of `nucl` and of that chromosome gave each read of the
  chromosome two P-sites. A compartment is `nucl` (all chromosomes that are
  not circular) or the name of a circular chromosome.
- The error of `prepare_for_ORFquant()` for a BAM file without a chromosome
  of the annotation now names the two files, and says that the reads must be
  aligned to the genome.

## Changes in results

- `prepare_for_ORFquant()`: the P-sites of the uniquely mapping reads
  (`P_sites_uniq`) and of the uniquely mapping reads with mismatches
  (`P_sites_uniq_mm`) now come from the correct reads. Before, a read with 2
  or more introns and its P-site in a middle block could get the mapping
  quality and the MD tag of another such read. In ORFquant, these P-sites
  give `P_sites_raw_uniq`, `P_sites_raw_uniq_mm` and `pval_uniq`, and with
  `unique_reads_only = TRUE` also the ORFs. This is rare: in a human sample,
  37 of 71.8 million P-sites of `P_sites_uniq` move. The P-sites of all reads
  (`P_sites_all`) and the results of our test data do not change. For this
  fix, `get_ps_fromspliceplus()` and `get_ps_fromsplicemin()` now give the
  P-sites in the order of the reads, not sorted by position.

## Bug fixes

- `prepare_annotation_files()` no longer stops on a GTF file without
  introns, for example a GTF file of chrM only, and `prepare_for_ORFquant()`
  no longer stops on an annotation without splice junctions.
- `prepare_annotation_files()` no longer stops when the chromosomes that are
  not circular have no genetic code or more than one, for example a GTF file
  of chrM only with chrM circular. The annotations that were made before do
  not change.
- `getSeq()` of a `FaFile_Circ` now reads a range that starts after the end
  of a circular chromosome. Before, such a range was dropped, so the
  sequences did not match the ranges, or R crashed.
- The messages of `load_annotation()` now say where it puts
  `GTF_annotation` and `genome_seq`: "the global workspace" or "the
  environment envir".

## Performance

- `prepare_for_ORFquant()` finds the P-sites of each chunk of the BAM file
  faster. The results don't change.

## Documentation

- The help of `prepare_annotation_files()` says that the order of the lines
  in the GTF file does not change the results if the lines agree, as in
  GENCODE and Ensembl files (lcalviell/ORFquant#10), and which annotations
  give wrong ORF categories in ORFquant: since ORFquant 1.3.1, only an
  annotation with `stop_in_gtf` `NA`.

# ORFquant 1.3.2

## Installation and dependencies

- ORFquant no longer needs the method of `unlist()` for `GRanges` that
  RiboseQC 1.4.0 and older have. This method changes `unlist()` in the whole
  R session, and a later RiboseQC removes it. `run_ORFquant()` used it to
  join the results of the regions, and `annotate_ORFs()` for an ORF in a
  region with the CDS of two or more genes. ORFquant now runs with and
  without it. The results do not change.

# ORFquant 1.3.1

## Changes in results

- `ORFs_spl_feat_longest` now removes the last codon of the CDS of `ref_id`
  only when it is a stop codon. A CDS can end without its stop codon (for
  example with the GENCODE tag `cds_end_NF`) also in a GTF whose other CDS
  include it. Since 1.03.0 (see "Changes in results" there), such a CDS lost
  its last 3 coding nucleotides, so an ORF that ends at its end, or a splice
  site in its last 3 nucleotides, got the wrong label. On a human sample
  (SRR15513199 with GENCODE 47), 29 of the 6,673 CDS of `ref_id` end without
  a stop codon. They are the reference of 33 ORFs. In the rows of these
  ORFs, 58 `ref` and 6 `missing_CDS` ranges get their last codon back, and 2
  labels change: `down_lastCDS` back to `same_lastCDS` (ADAMTSL1) and
  `down_3ss` back to `same_3ss`. The example data don't change.
- `ORF_category_Tx`, `ORF_category_Tx_compatible`, `ORF_category_Gen` and
  `NC_protein_isoform` now also keep the last codon of a CDS that ends
  without a stop codon, as `ORFs_spl_feat_longest` does (see above). Before,
  they compared the ORF with an end 3 nucleotides before the end of such a
  CDS, so an ORF that ends there was, for example, `C_extension` and
  `Alt3_stop`, not `ORF_annotated` and `exact_start_stop`. On the human
  sample, this changes 1 ORF (ADAMTSL1, which ends at the end of the CDS of
  its `ref_id`): its `ORF_category_Gen` changes from `Alt3_start_Alt3_stop`
  to `Alt3_start`, and its `NC_protein_isoform` from `N_C` to `N`. 6 of the
  8,830 CDS of the ORFs' transcripts and compatible transcripts end without
  a stop codon; the transcript categories of their 6 ORFs don't change. The
  example data don't change.

## Bug fixes

- `plot_ORFquant_results()` with coverage files no longer stops with "The
  `facets` argument of `facet_grid()` was deprecated in ggplot2 2.2.0 and
  is now defunct" when the results have a `CDS_spanning` splice event (an
  ORF exon over an intron of the CDS of its reference). ggplot2 4.0.0,
  which ORFquant needs, no longer accepts this argument (ggplot2 3 gave a
  warning), so the function stopped before it wrote any plot. The example
  data have no such event, but a whole human sample (SRR15513199 with
  GENCODE 47) has 15; it now gives all plots.
- `plot_ORFquant_results()` now stops at once with "The results in ... have
  no ORFs, so there is nothing to plot" when the results file has no ORFs.
  Before, it stopped later with `no rows to aggregate`. `run_ORFquant()`
  writes no results without ORFs, but a filtered results file can have none.
- `plot_orfquant_locus()` now names only the missing packages when it stops
  because packages it suggests are not installed. Before, it named all four.

## Documentation

- `?plot_ORFquant_results` now tells that with few ORFs, for example the
  results of a few genes, some plots are empty, and ggplot2 warns about
  groups with too few data points.
- `?select_quantify_ORFs` and `?run_ORFquant` now tell that
  `cutoff_P_sites` and `cutoff_cums` filter the values of each round of the
  selection, before the removed ORFs give their signal to the other ORFs. So
  a cutoff can remove an ORF whose value in the results of a run without
  this cutoff passes it.
- `?ORFquant_output` (section "P-site runs") now tells what a minimum of 5
  P-sites, not 5 ranges, would change in a human sample.

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
- `prepare_annotation_files()` no longer forges and installs a BSgenome
  package from a 2bit file: this failed with current BSgenomeForge versions
  (lcalviell/ORFquant#17, #19, #22 and #27). Give the genome as a FASTA file
  (`genome_seq`); `twoBitToFa` of the UCSC tools makes one from a 2bit file.
  `forge_BSgenome` now defaults to `FALSE`, and without `genome_seq` the
  function stops at once with this advice. `twobit_file`, `scientific_name`
  and `annotation_name` are no longer used. Annotations that earlier versions
  made with a forged BSgenome package still load.
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

- ORFquant now needs RiboseQC (1.4.0 or later), which it imports for the
  annotation step and `prepare_for_ORFquant()` (see "New features").
  RiboseQC is not on Bioconductor: `BiocManager::install("katarinagresova/ORFquant")`
  and `remotes::install_github()` install it from GitHub (the `Remotes` field
  of `DESCRIPTION`), and `environment.yml` has its dependencies. Through
  RiboseQC, ORFquant needs DT, ggpubr, scales, viridis and withr, and knitr
  and rmarkdown again. ORFquant itself no longer imports AnnotationDbi and
  txdbmaker, but RiboseQC does. `library(ORFquant)` attaches no new package.
- ORFquant now needs ggplot2 4.0.0 or later. `plot_ORFquant_results()` with
  coverage files no longer warns that the `size` aesthetic for lines is
  deprecated (ggplot2 3.4.0). Its violin plots use the arguments of ggplot2
  4.0.0 for the line at the median, not `draw_quantiles`, which is
  deprecated. With ggplot2 4.0.0 or later, the plots do not change.
- ORFquant no longer imports BSgenome: it used none of its functions.
  RiboseQC imports it, so it is still installed with ORFquant.
  `elementNROWS()` came into ORFquant through BSgenome; ORFquant now imports
  it from S4Vectors.

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
  - The chromosomes of the BAM file that are in the annotation must have the
    same lengths as in the annotation. Before, a BAM file of another genome
    assembly failed in the first chunk of reads with the unclear error
    `sequence chr22 has incompatible seqlengths`.
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
- `load_annotation()` now returns the annotation and the genome sequence,
  invisibly, as a list with the elements `GTF_annotation` and `genome_seq`,
  and `ORFquant()` takes them as its new arguments `annotation` and
  `genome_sequence`. So `ORFquant()` can run without global variables.
  Before, `load_annotation()` returned the genome sequence, invisibly, which
  was not documented. The global variables stay as they were, so old scripts
  work unchanged: `load_annotation()` still assigns `GTF_annotation` and
  `genome_seq`, `ORFquant()` without the new arguments reads them, and
  `run_ORFquant(interactive = TRUE)` still assigns `ORFquant_results` and
  `for_ORFquant`. `run_ORFquant()` and `prepare_for_ORFquant()` now use the
  annotation that they load, not the global variables. `R CMD check` no
  longer gives a note for these assignments.

- The annotation step is now the code of RiboseQC. `prepare_annotation_files()`,
  `load_annotation()` and the class `FaFile_Circ` (with its `getSeq()` and
  `seqinfo()` methods) come from RiboseQC, and ORFquant exports them
  again, so scripts that call them after `library(ORFquant)` still work.
  Both packages had a copy of this code, and the copies had drifted apart;
  the copy of ORFquant, with its fixes, went to RiboseQC, so each fix is now
  in one place. The annotation of a GTF is the same as before, apart from
  the changes listed here and under "Changes in results" and "Bug fixes".
  `_Rannot` files that earlier versions of ORFquant or RiboseQC made still
  load. New files hold the class `FaFile_Circ` of RiboseQC, and R needs
  RiboseQC to load them. What users can see:
  - The `*_similbed.bed` files in the annotation directory have a header line,
    as those of RiboseQC always had; for a set of regions that is empty, the
    file has only the header line. Before, ORFquant wrote no header, and an
    empty file for an empty set.
  - `prepare_annotation_files()` has the new argument `annot_file`, the path
    of the `_Rannot` file. Its default gives the path of before,
    `<annotation_directory>/<name of the GTF file>_Rannot`.
  - `load_annotation()` has the new argument `envir`, where it assigns
    `GTF_annotation` and `genome_seq`. Its default is the environment of the
    caller, not always the global environment. A script that calls it at the
    top level doesn't change. A function that calls `load_annotation()` and
    then `ORFquant()`, which by default reads the global variables, must give
    `envir = globalenv()`. `run_ORFquant()` and `prepare_for_ORFquant()` do,
    so they still assign the two global variables.
  - `load_annotation()` gives two messages: "assigning FaFile_Circ object
    <path> to the global workspace as genome_seq" and "assigning
    GTF_annotation object GTF_annotation to parent workspace".
  - Four messages of `prepare_annotation_files()` no longer say "ORFquant".
    They are the stop without `genome_seq` ("prepare_annotation_files() no
    longer forges a BSgenome package"), the stop for a GTF without CDS lines
    ("the annotation needs the coding sequences (CDS) of the transcripts"),
    and the two stops for a GTF without exon lines or without ids ("the
    annotation needs both on every line", "the annotation needs the exons and
    the coding sequences (CDS) of the transcripts").
  - `prepare_annotation_files()` no longer changes the random numbers of the
    session. Before, it called `set.seed(666)` when the TxDb had no exons for
    some transcripts (the message with 3 examples).
  - When `genome_seq` is a `FaFile_Circ`, the TxDb also marks the chromosomes
    of `circ_chroms` as circular. For a FASTA file, nothing changes.
- `get_ps_fromspliceplus()` and `get_ps_fromsplicemin()`, which give the
  P-sites of spliced reads, now come from RiboseQC too, and ORFquant exports
  them again. The two packages had the same code, so the P-sites don't
  change. Their help pages are now those of RiboseQC.
- `prepare_for_ORFquant()` now comes from RiboseQC (1.3.0), and ORFquant
  exports it again, so scripts that call it after `library(ORFquant)` still
  work. Its code went to RiboseQC unchanged, so the `for_ORFquant` files
  don't change. Its help page is now that of RiboseQC.

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
- With GTFs whose CDS lines end before the stop codon and that have no
  `stop_codon` lines, `run_ORFquant()` now takes the codon after each
  annotated CDS as its stop codon. These are GTFs with CDS lines made by
  factR (lcalviell/ORFquant#26), and GENCODE, Ensembl or RefSeq GTFs without
  their `stop_codon` lines. Before, ORFs were compared with annotated stop
  codons 3 nt too early: the annotated ORF itself was `C_extension` and
  `Alt3_stop`, not `ORF_annotated` and `exact_start_stop`, and
  `NC_protein_isoform` compared the wrong C-terminus.
  `prepare_annotation_files()` already found these GTFs: their `stop_in_gtf`
  is `NA`, because most CDS don't end with a stop codon. `annotate_ORFs()`
  now uses it, and `run_ORFquant()` says so in a message; annotations made
  before don't need to be made again. One rule is used for the whole
  annotation, so in a GTF with both kinds of CDS, ORFs at the stop codons of
  the smaller part still get wrong categories. On the data of
  lcalviell/ORFquant#26 (a StringTie and factR GTF, 7,399 ORFs), 6,918 ORFs
  had another annotation than with `stop_codon` lines added to the GTF, and
  no ORF was `ORF_annotated`; now all but 9 have the same. 5 of these 9 have
  a stop codon after an intron, and with `stop_codon` lines they got the
  wrong `ORF_category_Gen` (fixed too, see below). Which ORFs are found,
  their P-sites and p-values don't change. GTFs
  whose CDS lines include the stop codon or that have `stop_codon` lines,
  such as GENCODE, Ensembl, RefSeq, gffread and ORFanage GTFs, are not
  affected.
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
- `run_ORFquant(stn.orf_find.nostarts = TRUE)` now also finds ORFs without a
  start codon. Before, the option had no effect: `detect_translated_orfs()`
  called `get_orfs()` with `Stop_Stop = FALSE`. These ORFs are the parts of a
  frame that no ORF with a start codon covers: from the start of the
  transcript or a stop codon to the next stop codon, the first start codon or
  the end of the transcript, in frames that have an ORF with a start codon.
  So an N-terminal extension gives two ORFs: the part before the start codon
  and the ORF from the start codon. No column marks these ORFs. Their
  `Protein` does not start with M, except with `canonical_start_only = FALSE`
  when their first codon is TTG or CTG. They are selected and quantified
  with the other ORFs, so other ORFs of the same genes can change too. On the
  example data, one PPIL2 ORF without a start codon, before the start codon
  of another ORF, replaces another PPIL2 ORF. On a human sample (SRR15513199
  with GENCODE 47), 428 ORFs without a start codon join the 8,139 ORFs of the
  default, 30 of these 8,139 ORFs are gone (21 of them `ORF_annotated`), and
  448 others get other `P_sites`. 21 of the 25 genes of the 30 ORFs still
  have an ORF with a start codon that overlaps them. In NUCB2, CNOT2, TMEM98
  and EMC10, only ORFs without a start codon are left. There is no
  readthrough after an ORF without a stop codon after it; before, this
  stopped the run with `replacement has length zero`. With the default
  (`FALSE`), results don't change.
- In `ORFs_spl_feat_longest`, an ORF that overlaps no annotated CDS (its
  `ref_id` is `NA`), such as a uORF, is now compared with the longest
  annotated CDS of its gene, or with no CDS when its gene has none. Before,
  it was compared with the CDS used for the last ORF before it in its region
  that overlapped a CDS, which could be of another gene, or with no CDS when
  there was no such ORF: the result depended on the order of the ORFs. On a
  human sample (SRR15513199 with GENCODE 47), 379 of the 8,139 ORFs overlap
  no CDS, and 115 of them change their rows of `ORFs_spl_feat_longest`: 103
  in genes with a CDS, mostly uORFs, compared before with no CDS (102) or
  with the CDS of another gene (1), and 12 in genes without a CDS, such as
  pseudogenes, compared before with the CDS of another gene. `ref_id` stays
  `NA`, and all other outputs don't change.
- In `ORFs_spl_feat_longest` and `ORFs_spl_feat_maxORF`, ORFs on the minus
  strand now get the `spl_type` suffixes of the plus strand. New and missing
  exons at their 5' end now end with `_5prime` (for example
  `missing_CDS_5prime`, not `missing_CDS`), and when no exon of the ORF
  overlaps the reference, all its rows end with `_notoverl`, not `_3prime`.
  `annotate_splicing()` computed the `_5prime` labels of the minus strand but
  did not keep them. On a human sample (SRR15513199 with GENCODE 47), 4,298
  rows of 1,093 of the 8,139 ORFs change in `ORFs_spl_feat_longest`, and
  1,747 rows of 206 ORFs in `ORFs_spl_feat_maxORF`. The plus strand, the
  other outputs and the report don't change. `?annotate_splicing` now
  describes these suffixes.
- `ORF_category_Gen` now compares the last nucleotide before the stop codon
  of the ORF with that of the annotated CDS, as `ORF_category_Tx` does.
  Before, it compared the end of the CDS with the position 3 nt after the
  ORF's last exon, which is in an intron when the ORF's stop codon is split
  by the intron or is the first codon of an exon. Such an ORF with the
  annotated stop codon was `Alt5_stop` or `Alt3_start_Alt5_stop`, not
  `exact_start_stop` or `Alt3_start`. On a human sample (SRR15513199 with
  GENCODE 47), 22 of the 8,139 ORFs have such a stop codon, and 8 of them
  change: 5 from `Alt5_stop` to `exact_start_stop`, and 3 from
  `Alt3_start_Alt5_stop` to `Alt3_start`. The other 14 end at another place
  than the CDS of their `ref_id`, so their category stays the same. All other
  outputs don't change. GTFs whose CDS lines end before the stop codon
  (`stop_in_gtf` `NA`, see above) already got this rule and are not affected,
  so a GTF with and without its `stop_codon` lines now gives these ORFs the
  same category.
- With the default `stn.orf_quant.scaling = "total_Psites"`, the final
  scaling factor (`optim_feats`, the third value of `scaling_factors`) now
  puts all ORFs that share exonic bins or splice junctions, directly or
  through other ORFs, in one group. So the P-sites assigned to the ORFs of a
  group add up to the P-sites on their exonic bins, as `?ORFquant_output`
  says. Before, the groups were wrong: an ORF could be alone in its group, or
  in two groups. An ORF alone in its group got the final factor 1, so its
  `P_sites` was its `P_sites_raw`, and the P-sites that it shared with other
  ORFs counted two times. On the example data, the TCN2 uORF
  `ENST00000450638.5_250_336` had 22 P-sites (12.9 % of the gene), almost all
  of them those of the annotated ORF. Now it has less than 2 % and is no
  longer found. On a human sample (SRR15513199 with GENCODE 47), 8 of the
  8,139 ORFs are no longer found (3 `ORF_annotated`, 3 `novel` and 2 `uORF`,
  all ending before the annotated stop codon, such as `Alt5_stop`), 3 ORFs
  are found that were not before, and 131 ORFs in 51 genes get other
  `P_sites` (128 of them lower, 90 by more than 10 %). `ORFs_pM`, which is
  scaled over all ORFs, changes for all ORFs. In results of earlier
  versions, an ORF that was alone in a wrong group has a final factor of 1
  and a second factor (`adj_feats`) below 1, and the other ORFs of its gene
  can be wrong too: on this sample, 57 ORFs have this sign, and the fix
  removes or changes 44 of them. With `stn.orf_quant.scaling =
  "average_coverage"`, results don't change.
- `compatible_ORF_id_tr` is now always the ORF on `compatible_tx`. When an
  ORF is on more than one protein-coding transcript, ORFquant takes the
  first of them in sorted order. Before, it sorted the ORF ids and the
  transcript ids apart, so when a transcript id was the start of another
  (for example `MSTRG.1.1` and `MSTRG.1.10` in the C locale), the two
  columns could name two transcripts. A later step corrected the pair,
  except when the ORF had a protein-coding compatible transcript and its
  `longest_ORF` had none. Then `ORF_category_Tx_compatible` (also in the
  FASTA headers) mixed the two transcripts. `compatible_ORF_id_tr_longest`
  is chosen in the same way. Such ORFs are rare: GTFs from GENCODE, RefSeq,
  StringTie, IsoQuant and factR gave none.
- With `stn.orf_quant.scaling = "average_coverage"` (not the default), the
  coverage of a splice junction in the final scaling factor (`optim_feats`)
  is now its reads divided by 60, as in the steps before it. Before, its
  reads were divided by its width and then by 60, so the junctions added
  almost nothing to the coverage of the region, and the factor was too low.
  On the example data, `P_sites` of all 15 ORFs is now higher (for example
  98.7, not 77.9, for the annotated TCN2 ORF), and so `ORFs_pM` changes too.
  `ORF_pct_P_sites` and the ORFs that ORFquant finds don't change, as the
  factor is the same for all ORFs of a region. `?ORFquant_output` now
  describes this factor. With the default `"total_Psites"`, results don't
  change.
- `ORFs_spl_feat_longest` now compares the ORF's exons with the CDS without
  its stop codon, as the ORF ends before its stop codon. Before, with a GTF
  that has `stop_codon` lines (GENCODE, Ensembl and RefSeq GTFs do), the
  CDS included its stop codon. So the last exon of an ORF with the
  annotated stop codon was `up_lastCDS`, not `same_lastCDS`, and when an
  intron split the stop codon or came before it, the ORF also had a
  `missing_CDS_3prime` row of 1 to 3 nt. On a human sample (SRR15513199
  with GENCODE 47), 3,924 of the 4,125 `exact_start_stop` ORFs had
  `up_lastCDS`, and only 11 of the 8,134 ORFs had `same_lastCDS`. On the
  example data, the last exons of the annotated TCN2 and TOP3B ORFs are now
  `same_lastCDS`. GTFs whose CDS lines end before the stop codon
  (`stop_in_gtf` `NA`, see above) are not affected. `ORFs_spl_feat_maxORF`,
  the plots and the other outputs don't change.
- `NMD_candidate` is now `TRUE` only when the last exon-exon junction is
  downstream of the ORF's stop codon: `Distance_to_lastExEx` more than 3.
  Before, it was `TRUE` for each positive distance, but the distance is
  measured from the ORF's end, without its stop codon. So it was also `TRUE`
  when the stop codon was the first codon of the last exon (distance 1) or
  went across the last junction (distances 2 and 3). The same for
  `NMD_candidate_compatible_txs`. `Distance_to_lastExEx` doesn't change, so
  a filter on it for the usual minimum distance (often 50 nucleotides) gives
  the same ORFs. On a human sample (SRR15513199 with GENCODE 47), 16 of the
  946 ORFs with `NMD_candidate` `TRUE` had a distance of 1 to 3, and 19 of
  the 15,207 values of `NMD_candidate_compatible_txs`. On the example data,
  nothing changes.

- `prepare_annotation_files()` now chooses the representative start and stop
  codons of each gene from that gene's own codons. Before, when the start
  and stop codons of two genes alternated along a chromosome (for example in
  overlapping genes), the flags for a gene's most common and outermost
  codons could go to the other gene's codons. So `upstr_stasto`,
  `mostfreq_stasto` and the columns `reprentative_mostcommon`,
  `reprentative_boundaries` and `reprentative_5len` of `cds_txs_coords` (and
  of `cds_txs_coords_similbed.bed`) were wrong for some transcripts and
  genes: on GENCODE v47, 339 of 20,503 genes get another
  `reprentative_mostcommon` transcript. ORFquant doesn't read these columns,
  so the ORFs don't change. RiboseQC reads them for its metagene, so with a
  new annotation its results can change. Annotations that earlier versions
  made keep their values.
- `prepare_for_ORFquant()` (RiboseQC 1.4.0) now finds the uniquely mapping
  reads with mismatches (`P_sites_uniq_mm` in `<dest_name>_for_ORFquant`) by
  the mismatches in their MD tag: a number and then a base. Before, a read
  had mismatches when its MD tag had more than 3 characters. So a read
  shorter than 20 nt with 1 mismatch (for example `5A9`) did not count, and
  a read of 1,000 nt or more without mismatches did. Only these reads
  change, and they get a P-site only when the cutoff file has their read
  length. In ORFquant, these P-sites give only the column
  `P_sites_raw_uniq_mm`. The cutoff files of our test data have read lengths
  from 20 to 51 nt, so nothing changes there.

## Bug fixes

- `run_ORFquant()` no longer fails with `NA/NaN argument` when exactly one
  genomic region has ORFs, for example in a run on one gene
  (lcalviell/ORFquant#26).
- `run_ORFquant()` no longer drops the names of `ORFs_tx` (the ORF ids) when
  every region with ORFs (or every one in a block of 1000 such regions) has
  the same number of ORFs, 2 or more, for example in a run on one gene with
  2 ORFs. The `ORF_id_tr` column was not affected.
- `run_ORFquant()` no longer fails with `1 elements in value to replace 0
  elements` when the cutoffs of the quantification remove all ORFs of a
  genomic region, for example `stn.orf_quant.cutoff_P_sites = 10` on a region
  where all ORFs have fewer than 10 P-sites. The region now has no ORFs, as a
  region where no ORF is detected. When no region has ORFs, the run stops
  with `No ORFs found!`. Runs that did not fail give the same results.
- `run_ORFquant(write_GTF_file = FALSE)` now writes
  `<prefix>_Protein_sequences.fasta`, as `write_protein_fasta` (default
  `TRUE`) asks. Before, the FASTA file was written only together with the GTF
  file. Runs with the default values write the same files as before.
- `run_ORFquant()` with a `gene_id` value that is not in the annotation now
  stops with an error that names the value. Before, it stopped with an unclear
  error from the subsetting of `txs_gene`. A `gene_name` value that is not in
  the annotation now gives a warning that names it; as before, the other values
  are used. Runs with known genes give the same results.
- In `<prefix>_Detected_ORFs.gtf`, the `ORF_id` of each `CDS` line is now the
  ORF of the line. Before, all `CDS` lines of a transcript with more than one
  ORF had the `ORF_id` of one of them, for example a uORF and the annotated
  CDS of its transcript (lcalviell/ORFquant#13); their other values
  (`P_sites`, `ORF_pct_P_sites`, `ORF_pct_P_sites_pN`, `ORFs_pM`) were
  already right. On a human sample (SRR15513199 with GENCODE 47), 141 of
  7,994 transcripts with ORFs have more than one, and 1,232 `CDS` lines of 145
  ORFs had the id of another ORF. The ORFs, the other output files and the
  coordinates in the GTF don't change. All ORFs of a transcript still have
  its `transcript_id`, so tools that take one CDS for each transcript, such
  as a conversion to genePred, join them into one CDS: `?run_ORFquant` now
  describes the GTF and says to group its `CDS` lines by `ORF_id`.
- `annotate_splicing()` now gives the right `spl_type` when the exons of
  `orf_gen` are not in ascending order, for example minus-strand exons in
  transcript order. Before, it took each exon's overlap with the reference
  from another exon: it gave wrong labels, or failed with `failed to rbind
  column 'ref' across DataFrame objects`. `run_ORFquant()` always gives the
  exons in ascending order, so its results don't change.
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
- `prepare_for_ORFquant()` now sorts the P-sites from bigWig files in the
  order of the chromosomes of the annotation, as it sorts the P-sites from a
  BAM file. Before, they were not in this order when the bigWig files list
  the chromosomes in another order, or when a chromosome has P-sites on one
  strand only. The results of `run_ORFquant()` do not change.
- `plot_ORFquant_results()` no longer fails with `'breaks' are not unique`
  when the largest number of selected transcripts per gene is 3, 6 or 9
  (lcalviell/ORFquant#18). All other plots are unchanged.
- `plot_ORFquant_results()` with coverage files no longer fails when all
  regions of a type of splice site have the same value, for example when the
  type has one region (lcalviell/ORFquant#7). These regions are now one group,
  as when the type has two different values. Before, the function failed with
  `invalid number of intervals` when the value was below 2. When the value was
  2 or more, it failed with `Insufficient values in manual scale` if the type
  also had regions without a mixture of ORFs, and else it drew the line in a
  different red (`#B60048`, now `#6D0091`). All other plots are unchanged.
- `plot_ORFquant_results()` with coverage files no longer fails in the plots
  of the splice sites when a type of splice site has no regions (`1 elements
  in value to replace 0 elements`), has regions on one strand only (`subscript
  out of bounds`), or has one region to plot (`attempt to set 'rownames' on an
  object with no dimensions`). A type without regions gets the panel "Not
  enough distinct regions". On our 4 test data sets and on chromosome 21 of
  SRR15513199, the function failed; it now makes all plots.
- `plot_ORFquant_results()` with bedGraph coverage files no longer fails
  with `subscript contains out-of-bounds ranges` when a window around a
  splice site goes past the last covered base of its chromosome and strand
  (lcalviell/ORFquant#7). A bedGraph file has no chromosome lengths, so its
  coverage ended at its last covered base. The coverage now gets the
  chromosomes and their lengths of the P-sites file. BigWig files, and
  bedGraph files that did not fail, give the same plots as before.
- In the `CDS_spanning` plot of `plot_ORFquant_results()`, all `same_5ss`
  regions (the regions without a mixture of ORFs at 5' splice sites) are now
  in the downstream panel. Before, when there were more of them than of
  `same_3ss` regions, some were in the upstream panel; when there were fewer,
  the function failed with `attempt to set too many names`.
- In the `ORFs_maxiso_tpm` plot of `plot_ORFquant_results()` (section 3.4
  of the HTML report), each bin of `ORF_pct_P_sites` of the major ORF now
  has its violin of gene TPMs, with the bins of the `ORFs_maxiso` plot.
  Before, since R 4.0, all genes were in one violin labelled `NA`.
- The `Sel_txs_genes`, `Sel_txs_genes_tpm` and `Sel_txs_genes_pct` plots of
  `plot_ORFquant_results()` (sections 1.1 to 1.3 of the HTML report) now put
  the genes in the bins `1-3`, `4-6`, `7-9` and `>9` of their numbers of
  selected transcripts, in this order. Before, the labels of the bins were
  right only when the largest number of a gene was 49 (the last bin was
  `10-49`). For example, when it was 4, the genes with 4 selected
  transcripts had the bin `3-4`, and when it was 10, the bin `9-10`. The bars
  were in the order of their numbers of genes.
- The file `<prefix>_ORFquant_plots_RData` of `plot_ORFquant_results()` no
  longer holds the annotation, the results and the P-sites. Before, they
  were in the environment of each plot, which `save()` writes with it: in
  the file of a human sample (SRR15513199, GENCODE 47), they were 180 of its
  181 MB. `create_ORFquant_html_report()` loads one such file for each
  sample. On the example data, the file is now 3.5 MB, not 6.6 MB. The plots
  do not change.
- The log axes of the plots of `plot_ORFquant_results()` now have labels up
  to 1,000,000 (`ORFs_pM` and gene TPM) and 100,000 (ORF length). Before,
  they had none above 1000 and 10,000, so the upper part of the axis had no
  labels when the values were higher, for example on the example data. The
  y axis of `ORFs_found_ORFs_pM` now has the label `ORFs_pM`, not
  `P-sites_pNpM`.
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
    bright green (highest), as in the legend. Before, they were all white;
  - no longer warns that `qplot()` is deprecated (ggplot2 3.4.0), or that the
    font width of character 0x09 is unknown (the tabs in the name of the
    P-site track). The plots are unchanged;
  - warns when a call gives `bam_files` or `col`, because they were never
    used: the P-sites and junction reads come from the P-site file of the
    results, and these tracks are always forestgreen. `bam_files` is again
    optional;
  - no longer fails with `object 'selgene' not found` on Bioconductor 3.23.
    S4Vectors 0.50 evaluates the condition of `subset()` in another
    environment when the call is in a `%>%` pipe, so the condition did not
    find the variables of the function;
  - sets the option `ucscChromosomeNames` back to its value before the call
    when it ends. Before, it left it at `FALSE`, so Gviz no longer checked the
    chromosome names of the later plots of the session.

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
- `create_ORFquant_html_report()` now stops at once when an input file does
  not exist, or when `input_files` and `input_sample_names` have different
  lengths. Before, a missing file stopped the rendering with `error in
  evaluating the argument 'x' in selecting a method for function 'get'`. In
  section 3.5 of the report, a sample without the plot of the splice sites
  (`plot_ORFquant_results()` without coverage files) now has a note, not an
  empty tab. `?create_ORFquant_html_report` now says that the function needs
  pandoc, and that it writes the console output of the rendering to
  `<output_file>_ORFquant_report_output.txt`.
- The `getSeq()` error for ranges that wrap twice around a circular
  chromosome now reads "Ranges wrapping twice isn't implemented yet...",
  without the stray quote, line break and spaces it had.
- `run_ORFquant()` and `prepare_annotation_files()` now close the TxDb
  databases they make when they no longer need them. Their logs no longer
  show `Error in x$.self$finalize() : attempt to apply non-function` (often
  hundreds of lines in a whole-genome run) and `call dbDisconnect() when
  finished working with a connection` (lcalviell/ORFquant#3). These messages
  came from R closing the databases later, and did not change the results.
- `prepare_annotation_files()` now gives the message "N genes were dropped
  because they have exons located on both strands..." once, not twice.
- `from_tx_togen()` now stops with an error naming the ORFs that are not
  inside exactly one transcript of `exons`. Before, an ORF outside its
  transcript gave `subscript contains out-of-bounds indices`, and with a
  transcript name twice in `exons`, ORFs got the genomic coordinates of
  another ORF, without an error. `?from_tx_togen` now says that it maps the
  ORFs of one transcript: `exons` has one element and `introns` is a
  `GRanges`. `run_ORFquant()` calls it this way, with ORFs inside the
  transcript, so results don't change.
- `prepare_annotation_files()` no longer stops with `solving row <n>:
  'allow.nonnarrowing' is FALSE and the solved start (0) is < 1` when a
  transcript has a CDS shorter than 3 nt. It removes these transcripts from
  the annotation.
- `load_annotation()` also reads `_Rannot` files that early versions of
  RiboseQC made, with the name of a BSgenome package in `genome`. Before, it
  stopped with `'package' must be of length 1`.
- `prepare_annotation_files()` (RiboseQC 1.4.0) now stores the path of the
  FASTA file (`genome_seq`) as an absolute path. Before, it stored the path
  as given, so with a relative path, `load_annotation()` and
  `run_ORFquant()` failed in another working directory: the FASTA file did
  not exist. Annotations made before keep the path that they have.

## Performance

- On a human sample with a 1.8 GB BAM and GENCODE 47, the changes below
  together make a run take 3 h 11 min with 16 cores instead of 7 h 46 min
  (59% less; measured before the last change below), and at most 35 GB of
  memory instead of 45 GB. With 64 cores, a run takes 47 min and at most
  116 GB. Each entry below gives what its change saved when it was made.
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
- The genomic regions take less time again. With one core (timed with R's
  profiler on), a run on chr21 of a human sample takes 350 s instead of
  415 s (16% less), and a run on 160 random regions of the other chromosomes
  484 s instead of 563 s (14% less). With 64 cores, a run on the human sample
  with the 1.8 GB BAM takes 47 min instead of 57 min (17% less), and at most
  116 GB of memory instead of 115 GB. Results don't change: on our test data
  and on the whole sample with the 1.8 GB BAM, all outputs and the results of
  each region are the same as before, and on chr21 and the 160 regions, each
  changed function gives the same result as before in every call. Five
  changes make it faster:
  - `select_start()` and `calc_orf_pval()` took the P-sites of each ORF from
    the coverage of its transcript, with S4 operations for each ORF. They now
    take the coverage once as a plain vector, and a part of it for each ORF.
    `select_start()` chose the start codon of each stop codon with
    `endoapply()` on the ORFs split by stop codon; it now works on the ORFs'
    positions and takes the chosen ORFs once.
  - `select_quantify_ORFs()` built an SQLite database of the region's ORFs
    (a TxDb) to get their exonic parts with `exonicParts()`. It now gets the
    same parts directly, with `disjoin()`.
  - `from_tx_togen()` removed the introns from the ORFs one ORF at a time; it
    now does it for all ORFs of a transcript at once, with `psetdiff()`.
  - `annotate_splicing()` looked for the overlaps of an ORF's exons and the
    CDS exons three times; it now does it once. It also builds its list
    columns in one step.

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
  - `?get_orfs` said that the function loads the annotation. It finds the
    ORFs of a transcript sequence.
  - `?plot_orfquant_locus` said that `locus` is a gene name, but it is a
    `gene_id`. The page now also says that the function takes the annotation
    from the global variable `GTF_annotation`.
  - `?detect_readthrough` said that `P_sites_uniq_mm` is an `Rle`, but it
    is a `GRanges`, as `P_sites`.
  - `?select_txs` said that the function returns a `GRanges`, but a region
    with fewer than 4 P-sites and junction reads on its exonic bins and
    junctions gives an empty `GRangesList()`.
  - `?select_quantify_ORFs` described `optimiz`, but did not say that
    `run_ORFquant()` never sets it.
- The HTML report of `create_ORFquant_html_report()` said that the category
  "other" holds the overlapping uORFs, but the plots have them in the
  category "uORF". The report now says this, and names all categories in
  "other".
- `?select_start` and `?ORFquant_output` now describe the first filter of
  `select_start()`. Before it chooses the start codon of each stop codon, it
  removes the candidate ORFs below `cutoff_ave`
  (`stn.orf_find.start_sel_cutoff_ave`) and `cutoff`
  (`stn.orf_find.start_sel_cutoff`), but it compares them with the ORF's
  in-frame percentages, not fractions. So the default 0.5 removes only the
  candidates with an `ave_pct_fr` below 0.5, not 50: those with almost no
  P-sites in frame. ORFquant has always worked this way, and it still does,
  so that results don't change: with 50, 13 ORFs of a human sample change.
  4 of its 8,139 ORFs (3 annotated CDSs), with an `ave_pct_fr` of 49 to 52,
  are lost or replaced by an ORF that starts further downstream, and 7 more
  such N-terminally truncated ORFs are new. When
  `cutoff` selects a start codon, its `ave_pct_fr_st` and `pct_fr_st` are
  fractions, not percentages. `?select_start` also said that the start codon
  needs more than `cutoff_ave`; at least `cutoff_ave` is enough.
- `?calc_orf_pval`, `?take_Fvals_spect` and `?ORFquant_output` now say why
  ORFquant finds fewer short ORFs, such as uORFs (lcalviell/ORFquant#14).
  The 3-nucleotide periodicity test has less power for short ORFs. There is
  no minimum length, but an ORF shorter than 25 nt (8 codons or fewer) is
  padded with zeros to 50 nt before the test. In simulations (P-sites from a
  Poisson distribution, 70% in frame), with 30 P-sites in the ORF, the test
  gives a p-value below 0.05 for 95% of ORFs of 30 codons, 93% of ORFs of
  10 codons and 57% of ORFs of 3 codons. In a human sample (SRR15513199,
  8,139 ORFs), 32 of the 131 uORFs have 10 codons or fewer (the shortest
  has 3), but only 1 of the 49 ORFs of 12 codons or fewer has a `pval`
  below 0.001, against 51% of the longer ORFs. ORFquant has always worked
  this way, and it still does. The quantification cutoff
  `stn.orf_quant.cutoff_pct` does not take the ORF length into account, but
  on this sample it removed few uORFs: without it, there are 6 more (3 of
  them of 20 codons or fewer).
- `?prepare_for_ORFquant` and the README now say which alignments ORFquant
  uses. Alignments with an insertion or a deletion (`I` or `D` in the CIGAR)
  give no P-sites and no junction reads. Junction reads are counted only on
  the exon-exon junctions of the annotation; a read on another junction, for
  example of an isoform that is not in the GTF, still gives a P-site.
  ORFquant uses only the transcripts of the GTF, so isoforms to quantify, for
  example from long-read sequencing, must be in the GTF. ORFquant has always
  worked this way, and it still does. In two human samples aligned with STAR,
  0.7% (SRR15513199) and 0.03% (HCT116 cells) of the alignments have an
  insertion or a deletion, and 128 and 317 of them are on a CDS. 0.1% and
  0.5% of the junction reads are on junctions that are not in GENCODE 47;
  with a long-read annotation of the HCT116 cells, 0.4%.
- `?prepare_for_ORFquant` now says that soft-clipped bases (`S` in the CIGAR)
  are not counted: the read length of the cutoff table is the aligned length,
  and the cutoff counts from the first aligned base at the 5' end. So
  alignments with soft clipping (e.g. STAR's default) can give other read
  lengths and P-sites than end-to-end alignments, and the cutoff table must
  come from the same alignments (lcalviell/ORFquant#10). ORFquant has always
  worked this way.
- `?run_ORFquant` now says that a run on a few genes (`gene_name`, `gene_id`
  or `genomic_region`) gives the same ORFs and values as a run on all genes,
  except `ORFs_pM`, which sums to a million over the ORFs of the run.
  ORFquant analyzes each region of overlapping genes whole, with all its
  P-sites. On a human sample (SRR15513199), a run on 51 genes, and a run on
  each of 9 of them, gave the ORFs of the whole-genome run, with the same
  values in all other columns of `ORFs_tx`. A GTF or BAM file cut to a few
  genes can give other results, because a region then loses genes or their
  reads. In our test, 28 of 52 genes shared their region with genes that the
  cut removed. For 1 of them, the `P_sites` of its 2 ORFs changed by 0.01.
- `?ORFquant_output` and `?annotate_ORFs` now say how `ORF_category_Gen`
  labels an ORF that overlaps no CDS, but has exons on both sides of the CDS
  of the nearest gene and none between them, that is, an intron that holds
  this whole CDS. It is `novel_Upstream` on the plus strand, but
  `novel_Downstream` on the minus strand: on the plus strand, the label comes
  from the start of the ORF; on the minus strand, an ORF whose end is
  downstream of the CDS is `novel_Downstream`. For all other ORFs, the two
  rules give the same label. ORFquant has always worked this way, and it
  still does. No ORF is like this in a human sample (SRR15513199, with 148
  `novel_Upstream`, `novel_Downstream` and `novel_Internal` ORFs) or in the
  data of lcalviell/ORFquant#26 (139). `?ORFquant_output` also says that the
  nearest gene is the one with the CDS exon nearest to the ORF's first exon
  in genomic coordinates, and that an ORF is `novel_Internal` when one of its
  exons is between the first and the last CDS nucleotide.
- `?detect_readthrough` now describes the columns `compatible_id` and
  `compatible_original_id` of `ORFs_readthroughs` (beta). When ORFquant finds
  the same readthrough after the ORFs of several transcripts of a region, it
  usually keeps one of them, and these columns have the `ORF_id_tr` and the
  `ORF_orig_tr` of the others. Because of an indexing error in
  `detect_readthrough()`, a readthrough found after a readthrough with
  another protein in the same region has its own ids in these columns, and
  can lack the id of one other readthrough with the same protein. ORFquant
  has always worked this way, and it still does. No other column or output
  of ORFquant uses these two columns. In 22 human samples (GENCODE 47), 66 of
  2,470 readthroughs have their own id in `compatible_id`. A fix changes
  only these 66 values: 60 become empty, and 6 get the id of another
  readthrough.
- `?detect_readthrough` and `?get_reathr_seq` now describe the regions and
  the tests of `ORFs_readthroughs` (beta). A readthrough starts at the stop
  codon of its ORF, so its `Protein` starts with `*`, and ends one codon
  before the next stop codon: the amino acid of that codon is not in
  `Protein`, and its P-sites do not count. After the first region,
  `detect_readthrough()` extends a readthrough over each next region that
  passes the test alone, also over regions that fail, and does not test the
  readthrough after its last extension. So a readthrough over more than one
  stop codon can have a `pval` above 0.05. The tests need a `pct_fr` of 0.5,
  also when `stn.orf_find.cutoff_fr_ave` is lower. The last duplicate test
  compares `Protein` and `P_sites_raw` each on its own, not as pairs.
  ORFquant has always worked this way, and it still does. In 22 human samples
  (GENCODE 47), all 2,415 readthroughs over one stop codon pass the test, and
  16 of the 55 over more do not (one goes over 148 stop codons, with a `pval`
  of 0.50). With this version, one of these samples (SRR15513199) has 11
  readthroughs: 6 end one codon before a stop codon, and the other 5 at the
  end of the transcript.
- `?run_ORFquant` and `?get_orfs` now say what `canonical_start_only = FALSE`
  does. The alternative initiation codons of the genetic code (TTG and CTG in
  the standard code) are then start codons only at the first codon of a
  frame, that is, at one of the first three bases of a transcript: ORFquant
  translates each frame from its first codon, and Biostrings' `translate()`
  reads these codons as M only at the first codon. Inside a transcript, only
  the codons of M (ATG in the standard code) are start codons. So `FALSE`
  does not find ORFs with a near-cognate start codon, such as N-terminal
  extensions that start at a CTG. ORFquant has always worked this way, and it
  still does.
- `?ORFquant_output` now says how ORFquant counts the P-sites of runs. The
  P-site tracks of `prepare_for_ORFquant()` join neighbouring positions with
  the same number of P-sites into one range, and three steps use such a run
  as one unit. `detect_translated_orfs()` maps a run to a transcript only
  when the run is inside one exon, so a transcript loses the P-sites of a run
  that crosses the edge of one of its exons, or, depending on the other
  selected transcripts, gets them at shifted positions. `select_txs()` and
  `select_quantify_ORFs()` count all the P-sites of a run on each exonic bin
  that it overlaps. `ORFquant()` analyses a region only when more than 4
  ranges overlap it. ORFquant has always worked this way, and it still does.
  In a human sample (SRR15513199), 5% of the ranges are runs of 2 to 6
  positions, 14 of 8,139 ORFs have a `P_sites_raw` 1 or 2 too low, and 270
  of 939,002 exonic bins count too many P-sites (0.038% more in total). On
  43 genes of this sample, single positions give one more ORF (RPL19, with
  an `ORF_pct_P_sites` of 2.04), and lower the `ORF_pct_P_sites` of the two
  other RPL19 ORFs (from 3.4 to 2.6 and from 96.6 to 95.3). The effect of
  single positions on all ORFs of the sample is not known.
- `?select_txs` now says how it removes nested transcripts when transcripts
  have the same features (exonic bins and junctions), such as transcripts
  with the same exons. Of 2 such transcripts it keeps one, but of 3 it keeps
  2, of 2 that are nested in a third transcript it keeps one with the third,
  and of 2 that have all features of the region it keeps both. The selection
  can then take all of them, so `selected_txs` and the transcripts of
  `<prefix>_Detected_ORFs.gtf` can have more than one transcript with the
  same exons. An ORF on them is in the results once, with each of them in
  `compatible_with`. ORFquant has always worked this way, and it still does.
  In a human sample (SRR15513199, GENCODE 47), 5 of the 156 groups of
  transcripts with the same exons have 2 selected transcripts, each a gene
  with only these 2 transcripts and one ORF.
- `?annotate_splicing` now explains the labels `same`, `up` and `down` of
  `spl_type` in `ORFs_spl_feat_longest` and `ORFs_spl_feat_maxORF`: each
  end of an ORF exon is at, upstream of or downstream of the same end of
  the reference exon, on the strand of the ORF. It also describes the
  column `ref`, and says that the reference is the CDS of `ref_id` without
  its stop codon, or the ORF of `ref_id_maxORF`.
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
