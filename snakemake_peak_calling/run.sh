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
SPECIES=Nmel

CONFIGFILE=${PWD}/config_files/config_${SPECIES}.yml
BASEDIR=/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_snakemake_output
DIR=${BASEDIR}/${SPECIES}
mkdir -p ${DIR}

cd ${DIR}
source /Genomics/argo/users/bjarnold/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
# for --singularity-args, list all parent directories containing subdirectories that have data or files the snakemake workflow needs to access
snakemake --snakefile /Genomics/kocherlab/bjarnold/STARRseq/code/snakemake_peak_calling/Snakefile \
--configfile ${CONFIGFILE} \
-p --use-singularity \
--singularity-args "--bind /Genomics/kocherlab/bjarnold --bind /Genomics/kocherlab/lab/data/GenomeReleases/official_release_v2.1.1 --bind /Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1" \
--cores 20 \
--rerun-triggers mtime