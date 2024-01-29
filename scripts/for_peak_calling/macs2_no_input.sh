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

conda activate macs2

LIBRARY=Amel
BIOREP=F3
BASEDIR=/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_snakemake_output/${LIBRARY}

INPUT_RNA=${BASEDIR}/BAMs_sorted/Amel-${BIOREP}-RNA_csorted_dedup.bam
NAME=/Genomics/kocherlab/bjarnold/STARRseq/data/MACS2_no_input/${LIBRARY}-${BIOREP}
GSIZE=223937270

macs2 callpeak \
-t ${INPUT_RNA} \
-f BAMPE \
-g ${GSIZE} \
-n ${NAME} \
--bdg \
-q 1 \
--min-length 100