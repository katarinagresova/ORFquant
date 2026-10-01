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
