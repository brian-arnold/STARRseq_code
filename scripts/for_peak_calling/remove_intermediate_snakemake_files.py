#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

d = "/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_snakemake_output/old"

species = ['Nmel', 'Apur', 'Avir' , 'Aaur']

for s in species:
    to_remove = glob.glob(f'{d}/{s}/FASTQ_trimmed/*')
    to_remove.extend(glob.glob(f'{d}/{s}/SAMs/*'))
    for r in to_remove:
        # print(r)
        os.remove(r)

