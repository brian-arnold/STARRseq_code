#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

def main():

  assembly_dir = "/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1"
  # get list of all subdirectories in assembly_dir
  species_list = os.listdir(assembly_dir)
  print(species_list)
  for i,spec in enumerate(species_list):
    # if i > 0:
    #   continue
    dir = f'{assembly_dir}/{spec}'
    fasta_file = f'{dir}/{spec}_genome_v3.fasta.gz'
    # check if fasta file exists
    if os.path.isfile(fasta_file):
      if not os.path.isfile(f'{dir}/{spec}_genome_v3.fasta.gz.bwt'):
        print(f'Indexing {fasta_file}')
        os.system(f'bwa index {fasta_file}')
      else:
        print(f'skipping {spec}, fasta file already indexed')
    else:
      print(f"fasta file doesn't exist for {spec}")


if __name__ == '__main__':
  main()
