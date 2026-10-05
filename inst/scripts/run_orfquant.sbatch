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
## the resources below on the sbatch command line (e.g. sbatch -c 64 --mem=160G)
## and add --partition or --account if your cluster needs them.
##
## Resources, measured on human samples with GENCODE 47: a sample with a 1.8 GB
## BAM, reusing the annotation, took 3 h 11 min and at most 35 GB of memory with
## 16 CPUs (the default below), and 57 min and 115 GB with 64 CPUs. So if
## your cluster has them, ask for more CPUs, as in the example above. Each CPU
## works in a forked copy of R, so memory grows with the CPUs. Before the
## speed-ups listed in NEWS.md, with 16 CPUs, that sample took 7 h 46 min and
## 45 GB, and one with a 69 MB BAM 1 h 52 min and 29 GB (8 min and 5 GB of it
## for the annotation), or 59 min and 53 GB with 32 CPUs (reusing the
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
# One BLAS thread per R process: run_ORFquant() forks one process per CPU, and
# a multi-threaded BLAS (such as conda's OpenBLAS) would start as many threads
# as the job has CPUs in each of them. It reads these when R starts.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
script=$(Rscript -e 'cat(system.file("scripts", "run_orfquant.R", package = "ORFquant"))')
if [ -z "$script" ]; then
  echo "Error: ORFquant is not installed in the R of $(command -v Rscript)" >&2
  exit 1
fi
Rscript "$script" --cores "$SLURM_CPUS_PER_TASK" "$@"
