# ORFquant 1.03.0

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
