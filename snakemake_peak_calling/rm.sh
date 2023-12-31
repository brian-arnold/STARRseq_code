#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

source /Genomics/argo/users/bjarnold/miniconda3/etc/profile.d/conda.sh

rm -r FASTQ_trimmed
rm -r fastp_output
rm -r SAMs
rm -r BAMs_sorted
rm -r BAMs_nsorted
rm -r MACS2
rm -r deeptools
rm -r genrich_single
rm -r genrich_multi
