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

conda activate deepstarr

SPECIES=Apur
DATA_DIR=/Genomics/kocherlab/bjarnold/STARRseq/data/deepSTARR_input/min_biorep_support_3/qval_0_folddiff_2/${SPECIES}
FASTA=MACS2_peaks_dedup.fa # MACS2_peaks_raw.fa, MACS2_peaks_dedup.fa, random_seqs.fa
# FASTA=MACS2_peaks_raw.fa # MACS2_peaks_raw.fa, MACS2_peaks_dedup.fa, random_seqs.fa
# FASTA=random_seqs.fa # MACS2_peaks_raw.fa, MACS2_peaks_dedup.fa, random_seqs.fa


DEEPSTARR_DIR=/Genomics/kocherlab/bjarnold/DeepSTARR_github_repo
MODEL=${DEEPSTARR_DIR}/DeepSTARR.model
DEEPSTARR_PYTHON_SCRIPT=${DEEPSTARR_DIR}/DeepSTARR/DeepSTARR/DeepSTARR_pred_new_sequence.py

python ${DEEPSTARR_PYTHON_SCRIPT} \
--model ${MODEL} \
--seq ${DATA_DIR}/${FASTA}