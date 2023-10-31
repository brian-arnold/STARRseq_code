#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict

def main():

  base_dir = '/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_snakemake_output'
  species = ['Nmel', 'Aaur', 'Apur', 'Avir']

  out_dir = '/Genomics/kocherlab/bjarnold/STARRseq/data/peak_calling_final'
  
  for s in species:
    cmd1 = f'cp -r {base_dir}/{s}/MACS2 {out_dir}/{s}'
    cmd2 = f'cp -r {base_dir}/{s}/genrich_single {out_dir}/{s}'

    os.system(cmd1)
    os.system(cmd2)

if __name__ == '__main__':
  main()
