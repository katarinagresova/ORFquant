# Example data

A small data set for the vignette, the tests and the command-line script
`inst/scripts/run_orfquant.R`. It covers 7 genes of human chr22 (GRCh38):
*TXNRD2*, *RANBP1*, *PPIL2*, *TOP3B*, *DDT*, *TCN2* and *PISD*.
`inst/scripts/make_example_data.R` describes how the files were made.

| File | Content | Source |
|---|---|---|
| `chr22_example.gtf.gz` | All GENCODE release 47 lines of the 7 genes, unchanged. | GENCODE [1] |
| `chr22_example.fa.gz`, `.fai`, `.gzi` | The sequence of chr22 from the GRCh38 primary assembly (GENCODE's `GRCh38.primary_assembly.genome.fa`). Every base more than 1 kb from the 7 genes is replaced by N. Compressed with `bgzip` and indexed. | GRCh38, Genome Reference Consortium [2] |
| `chr22_example.bam`, `.bai` | The primary alignments of the Ribo-seq reads of SRA run SRR15513199 within 500 bp of the 7 genes. The alignment used STAR 2.7.10a and the GENCODE 47 annotation. Only the MD and NH tags are kept. | Chothani et al. [3] |
| `chr22_example_cutoffs.tsv` | The P-site offset for each read length of this sample: the offsets that Ribo-seQC used for its P-sites of the whole sample (with `rescue_all_rls = TRUE`), that is, those of the read lengths it selected, and 12 for all others. | Ribo-seQC [4] |

SRR15513199 is GEO sample GSM5527724: untreated human aortic endothelial cells
(GEO series GSE182371, BioProject PRJNA756018).

## Attribution and terms of use

These files are not covered by the GPL licence of the ORFquant code. They stay
under the terms of their sources:

- **Reads.** The data are open access in the NCBI SRA and GEO. NCBI places no
  restriction on their use or distribution. The submitters can claim rights in
  the data, and NCBI cannot grant permission for them. The BAM file holds a
  small subset of the reads of the sample, with attribution to the authors [3].
- **GENCODE annotation and genome sequence.** GENCODE data are open access.
  EMBL-EBI places no restriction on their use or redistribution other than
  those of the original data owners, and expects attribution [1].

None of the sources states a licence. If you use these files in a publication,
cite the papers below.

## References

1. Mudge JM et al. (2025). GENCODE 2025: reference gene annotation for human
   and mouse. *Nucleic Acids Research* 53, D966–D975.
   doi:10.1093/nar/gkae1078
2. Schneider VA et al. (2017). Evaluation of GRCh38 and de novo haploid genome
   assemblies demonstrates the enduring quality of the reference assembly.
   *Genome Research* 27, 849–864. doi:10.1101/gr.213611.116
3. Chothani SP et al. (2022). A high-resolution map of human RNA translation.
   *Molecular Cell* 82, 2885–2899. doi:10.1016/j.molcel.2022.06.023
4. Calviello L, Sydow D, Harnett D, Ohler U (2019). Ribo-seQC: comprehensive
   analysis of cytoplasmic and organellar ribosome profiling data. *bioRxiv*
   601468. doi:10.1101/601468
