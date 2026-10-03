#!/bin/bash
## SLURM job running ORFquant on one sample with run_orfquant.R, installed with
## the package next to this file. Its arguments are passed to run_orfquant.R,
## with --cores set to the job's CPUs:
##
##   sbatch run_orfquant.sbatch --gtf genes.gtf --fasta genome.fa \
##     --bam sample.bam --offsets sample_cutoffs.tsv --outdir results
##
## Submit it from a shell in which Rscript is the R with ORFquant (e.g. after
## conda activate orfquant): the job inherits the shell's environment. Change
## the resources below on the sbatch command line (e.g. sbatch -c 32 --mem=96G)
## and add --partition or --account if your cluster needs them.
##
## Resources, measured on human samples with GENCODE 47 and 16 CPUs: a sample
## with a 69 MB BAM took 1 h 52 min and at most 29 GB of memory (8 min and 5 GB
## of it for the annotation); one with a 1.8 GB BAM, reusing the annotation,
## 7 h 46 min and 45 GB. Each CPU works in a forked copy of R, so memory grows
## with the CPUs: with 32, the first sample took 59 min and 53 GB (reusing the
## annotation). The annotation also writes up to 0.25 GB of temporary files, in
## TMPDIR (by default /tmp): if the job's /tmp is smaller, ask for more
## temporary space (--tmp, or what your cluster uses) or set TMPDIR to a larger
## directory.
##
## For more samples, reuse the first one's annotation once its job is done:
##
##   sbatch --dependency=afterok:<its job id> run_orfquant.sbatch \
##     --annotation results/annotation/genes.gtf_Rannot \
##     --bam sample2.bam --offsets sample2_cutoffs.tsv --outdir results
#SBATCH --job-name=orfquant
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
script=$(Rscript -e 'cat(system.file("scripts", "run_orfquant.R", package = "ORFquant"))')
if [ -z "$script" ]; then
  echo "Error: ORFquant is not installed in the R of $(command -v Rscript)" >&2
  exit 1
fi
Rscript "$script" --cores "$SLURM_CPUS_PER_TASK" "$@"
