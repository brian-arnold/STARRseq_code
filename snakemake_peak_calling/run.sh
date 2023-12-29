#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=60        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G         # memory per cpu-core (4G is default)
#SBATCH --time 3-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

######
# CHANGE OUT/ERR FILES TOO
######
# single-pop:
# DONE: Nmel, Avir, Apur, Aaur
# TO-DO: Amel, Hlig, Lzep, Hqua, Lbal, Lvie, Bimp

# multi-pop: 
# TO-DO: Lalb, Sinv


SPECIES=Dmel

CONFIGFILE=${PWD}/config_files/config_${SPECIES}.yml
RUNDIR=${PWD}
RESULTSDIR=/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_snakemake_output
DIR=${RESULTSDIR}/${SPECIES}
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
--rerun-triggers mtime \
2> ${RUNDIR}/err_${SPECIES} > ${RUNDIR}/out_${SPECIES}

# --rerun-triggers mtime \
# --rerun-incomplete \