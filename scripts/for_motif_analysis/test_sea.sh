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

SPECIES=Nmel
TYPE=dedup

MEME_SUITE_DIR=/Genomics/argo/users/bjarnold/meme/bin
DATA_DIR=/Genomics/kocherlab/bjarnold/STARRseq/data/meme_suite/fastas/${SPECIES}
DB=/Genomics/kocherlab/bjarnold/STARRseq/data/meme_suite/B10_Matrices_MEME.txt
OUT_DIR=/Genomics/kocherlab/bjarnold/STARRseq/data/meme_suite/sea_output/${SPECIES}/${TYPE}

${MEME_SUITE_DIR}/sea \
--p ${DATA_DIR}/MACS2_triplet_peaks_${TYPE}.fa \
--m ${DB} \
--oc ${OUT_DIR}
