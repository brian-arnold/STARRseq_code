#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

def main():

  new_assembly_dir="/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1"
  old_assembly_dir="/Genomics/kocherlab/bmjones/genomes"

  # get list of species, which correspond to the subirectories in the new assembly directory
  subdirs = glob.glob(f'{new_assembly_dir}/*')
  species = [s.split('/')[-1] for s in subdirs]
  # print(species)

  for s in species:
    # s is a 4 character string, but convert the last 3 characters to lower case
    s_lower = s[:-3] + s[-3:].lower()
    ref_fastas_paths = glob.glob(f'{new_assembly_dir}/{s}/*.fasta.gz')
    ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
    ref_fasta = [f for f in ref_fastas if '.alt.' not in f][0]

    ref = f'{s_lower}\t{new_assembly_dir}/{s}/{ref_fasta}'
    print(ref)

if __name__ == '__main__':
  main()
