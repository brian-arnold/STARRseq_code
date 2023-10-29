#!/bin/bash
#SBATCH -J template
#SBATCH -o out_Nmel
#SBATCH -e err_Nmel
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=20        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=3G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

######
# CHANGE OUT/ERR FILES TOO
######

DIR=/Genomics/kocherlab/bjarnold/Nmel
CONFIGFILE=${PWD}/config_files/config_Nmel.yml

# DIR=/Genomics/kocherlab/bjarnold/Aaur
# CONFIGFILE=${PWD}/config_files/config_Aaur.yml

# DIR=/Genomics/kocherlab/bjarnold/Apur
# CONFIGFILE=${PWD}/config_files/config_Apur.yml

# DIR=/Genomics/kocherlab/bjarnold/Avir
# CONFIGFILE=${PWD}/config_files/config_Avir.yml

cd ${DIR}
source /Genomics/argo/users/bjarnold/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
# for --singularity-args, list all parent directories containing subdirectories that have data or files the snakemake workflow needs to access
snakemake --snakefile /Genomics/kocherlab/bjarnold/snakemake/Snakefile \
--configfile ${CONFIGFILE} \
-p --use-singularity \
--singularity-args "--bind /Genomics/kocherlab/bjarnold --bind /Genomics/kocherlab/bmjones --bind /Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1" \
--cores 20 \
--rerun-triggers mtime 