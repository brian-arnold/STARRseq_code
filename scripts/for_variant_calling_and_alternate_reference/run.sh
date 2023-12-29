#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

source /Genomics/argo/users/bjarnold/miniforge3/etc/profile.d/conda.sh

# conda activate bioinformatics
# python3 compute_depth_per_sample.py
# python3 filter_vcfs.py

# conda activate cyvcf
# python3 03_filter_vcfs_within_sample_allele_frequency.py Aaur
# python3 03_filter_vcfs_within_sample_allele_frequency.py Apur
# python3 03_filter_vcfs_within_sample_allele_frequency.py Avir
# python3 03_filter_vcfs_within_sample_allele_frequency.py Nmel
# python3 03_filter_vcfs_within_sample_allele_frequency.py Amel
# python3 03_filter_vcfs_within_sample_allele_frequency.py Bimp

conda activate bioinformatics
python3 04_make_alternate_reference.py