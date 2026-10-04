# ORFquant
An R package for Splice-aware quantification of translation using Ribo-seq data


*ORFquant* is an R package that aims at detecting and quantifiying ORF translation on complex transcriptomes using Ribo-seq data.
This package uses syntax and functions present in Bioconductor packages like *GenomicFeatures*, *rtracklayer* or *BSgenome*. 
*ORFquant* aims at quantifying translation at the single ORF level taking into account the presence of multiple transcripts expressed by each gene.
To do so, the *ORFquant* pipeline consists of transcript filtering, *de-novo* ORF finding, ORF quantification and ORF annotation.
A variety of annotation methods, both in transcript and genomic space, is performed for each ORF, to yield a more complete picture of alternative splice sites usage, uORF translation, translation on NMD candidates and more.

This repository is a fork of [lcalviell/ORFquant](https://github.com/lcalviell/ORFquant) that installs and runs on current Bioconductor. [NEWS.md](NEWS.md) lists what changed, including fixes that change some results.

More details can be found in our manuscript:

### Quantification of translation uncovers the functions of the alternative transcriptome ###

*Lorenzo Calviello^, Antje Hirsekorn, Uwe Ohler^*

**Nature Structural and Molecular Biology** (2020), doi: https://doi.org/10.1038/s41594-020-0450-4

Preprint: **bioRxiv** (2019), doi: https://doi.org/10.1101/608794


## Installation

*ORFquant* needs R 4.5 or later and Bioconductor 3.22 or later.

With *BiocManager*:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("katarinagresova/ORFquant", dependencies = TRUE)
```

`dependencies = TRUE` also installs the packages that `plot_orfquant_locus()` and `create_ORFquant_html_report()` need (*Gviz*, *lemon*, *dplyr*, *GenomeInfoDb*, *knitr*, *rmarkdown*). The HTML report also needs [pandoc](https://pandoc.org). On Linux, the packages are usually built from source, which needs the development files of system libraries such as libcurl and libjpeg (`libcurl-devel` and `libjpeg-turbo-devel` on RHEL or Fedora, `libcurl4-openssl-dev` and `libjpeg-dev` on Debian or Ubuntu).

Or with conda, using [environment.yml](environment.yml) from this repository (R 4.5, Bioconductor 3.22, all dependencies and pandoc):

```bash
conda env create -f environment.yml
conda activate orfquant
Rscript -e 'remotes::install_github("katarinagresova/ORFquant", upgrade = "never")'
```

To also install the vignette, add `build_vignettes = TRUE` to either install command (this needs pandoc), and read it with `vignette("ORFquant")`.


## Usage

Three steps are required to use *ORFquant* on your data. Run each one in a new R session (e.g. with `Rscript`): they load the annotation into global variables.

```r
library(ORFquant)

# 1. Once per annotation and genome: writes annotation/genes.gtf_Rannot
prepare_annotation_files(annotation_directory = "annotation",
                         gtf_file = "genes.gtf",
                         genome_seq = "genome.fa",
                         forge_BSgenome = FALSE)

# 2. Once per sample: P-sites and junction reads, written to sample_for_ORFquant
prepare_for_ORFquant(annotation_file = "annotation/genes.gtf_Rannot",
                     bam_file = "sample.bam",
                     path_to_rl_cutoff_file = "sample_cutoffs.tsv",
                     dest_name = "sample")

# 3. Find and quantify ORFs: writes sample_final_ORFquant_results,
#    sample_Detected_ORFs.gtf, sample_Detected_ORFs.tsv and
#    sample_Protein_sequences.fasta
run_ORFquant(for_ORFquant_file = "sample_for_ORFquant",
             annotation_file = "annotation/genes.gtf_Rannot",
             n_cores = 4,
             prefix = "sample")
```

The [Ribo-seQC](https://github.com/lcalviell/Ribo-seQC) package can also create the input of step 3 from a Ribo-seq BAM file.

`sample_final_ORFquant_results` holds a list, which `get(load("sample_final_ORFquant_results"))` returns. Its `ORFs_tx` has one range per ORF, in transcript coordinates, with the ORF's P-sites, p-values, categories and `ORFs_pM` (P-sites per ORF length, scaled to sum to a million, akin to TPM); `ORFs_gen` has the ORFs' genomic coordinates; `selected_txs` lists the transcripts selected for quantification. `sample_Detected_ORFs.tsv` has the `ORFs_tx` table, one row per ORF, for use without R. `?ORFquant_output` describes each of its columns and the ORF categories; see also `?run_ORFquant` and `?ORFquant`.

Plots and an HTML report of the results:

```r
plot_ORFquant_results(for_ORFquant_file = "sample_for_ORFquant",
                      ORFquant_output_file = "sample_final_ORFquant_results",
                      annotation_file = "annotation/genes.gtf_Rannot")
create_ORFquant_html_report(input_files = "sample_final_ORFquant_results_plots/sample_ORFquant_plots_RData",
                            input_sample_names = "sample",
                            output_file = "sample_ORFquant_report.html")
```

The script `run_orfquant.R`, installed with the package, runs the three steps from the command line, each in a new R process, with their default parameters:

```sh
Rscript $(Rscript -e 'cat(system.file("scripts", "run_orfquant.R", package = "ORFquant"))') \
  --gtf genes.gtf --fasta genome.fa --bam sample.bam --offsets sample_cutoffs.tsv \
  --outdir results --cores 4
```

It writes the annotation to `results/annotation` and the other files to `results/sample_*`, as above, and lists them at the end. For other samples, `--annotation results/annotation/genes.gtf_Rannot` instead of `--gtf` and `--fasta` reuses the annotation. `--gene-names` and `--gene-ids` restrict the analysis to the genomic regions of some genes; `--help` lists all options.

On a SLURM cluster, `run_orfquant.sbatch`, installed next to it, runs it as a job, with `--cores` set to the job's CPUs:

```sh
sbatch $(Rscript -e 'cat(system.file("scripts", "run_orfquant.sbatch", package = "ORFquant"))') \
  --gtf genes.gtf --fasta genome.fa --bam sample.bam --offsets sample_cutoffs.tsv \
  --outdir results
```

It asks for 16 CPUs, 64 GB of memory and 24 hours: a human sample with a 1.8 GB BAM took 4 h 3 min and 37 GB. With 64 CPUs (`sbatch -c 64 --mem=160G`, before the script's path) it took 1 h 27 min and 114 GB, so ask for more CPUs if your cluster has them. Its comments give more measurements and say how to change these resources.

The [vignette](vignettes/ORFquant.Rmd) runs all these steps on example data included in the package, 7 genes of human chr22, and shows the results.


### Input files

- **GTF**: `exon` and `CDS` lines, each with `transcript_id` and `gene_id`. `gene_name`, `gene_biotype` (or `gene_type`) and `transcript_biotype` (or `transcript_type`) are used when present.
- **Genome sequence**: a FASTA file indexed with `samtools faidx`, with the same chromosome names as the GTF and the BAM. The `_Rannot` file refers to it by its path, so keep it there. Chromosomes named as in `circ_chroms` (`chrM`, `MT` and other mitochondrial names by default) are treated as circular.
- **BAM**: Ribo-seq reads aligned to the genome with a spliced aligner (e.g. STAR), sorted by coordinate and indexed with `samtools index`. Reads flagged as duplicates and secondary alignments are skipped. Reads with a mapping quality above 50 count as uniquely mapping (STAR gives them 255).
- **Cutoff table**: the P-site offset for each read length, as a tab-separated table with the columns `read_length`, `cutoff` and `compartment`, in any order. `cutoff` is the distance of the P-site from the read's 5' end; `compartment` is `nucl` for all chromosomes that are not circular, or the name of a circular one. For example:

  ```
  read_length	cutoff	compartment
  28	12	nucl
  29	12	nucl
  ```

  Instead of a cutoff table, `prepare_for_ORFquant()` also takes P-site bigWig files, one per strand (`path_to_P_sites_plus_bw` and `path_to_P_sites_minus_bw`). It still reads the BAM file, for the junction reads.


### Notes

- `prepare_annotation_files()` can also forge and install a *BSgenome* package from a 2bit file (`twobit_file`, `forge_BSgenome = TRUE`, the default). This fails on Bioconductor 3.22 (see lcalviell/ORFquant#17, #19, #22 and #27), so use `genome_seq` and `forge_BSgenome = FALSE` as above.
- `library(ORFquant)` attaches only *GenomicRanges* (and the packages it attaches). Scripts that use, for example, *Biostrings*, *rtracklayer* or *ggplot2* need their own `library()` calls.
- `n_cores` above 1 uses forked processes (`parallel::makeForkCluster()`), which need Linux or macOS.


For any question, please email:

calviello.l.bio@gmail.com or uwe.ohler@mdc-berlin.de

For problems with this fork, please open an issue at https://github.com/katarinagresova/ORFquant/issues.


Enjoy!
