#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=8G         # memory per cpu-core (4G is default)
#SBATCH --time 0-01:00:00        # DAYS-HOURS:MINUTES:SECONDS

source /Genomics/argo/users/bjarnold/miniforge3/etc/profile.d/conda.sh
conda activate bedtools

REF=/Genomics/kocherlab/bjarnold/STARRseq/data/alternate_references/renamed/Nmel.fasta
BED=./kmers_test.txt

bedtools getfasta -fi ${REF} -bed ${BED} -fo ./bedtools_tmp.fa

