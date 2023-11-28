#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

def main():

  # species_list = ["Nmel", "Aaur", "Apur", "Avir"]
  species_list = ["Amel", "Lzep", "Sinv", "Lvie", "Hqua", "Hlig", "Bimp", "Lbal"]
  
  beryl_dir = "/Genomics/kocherlab/bmjones/STARRseq"
  beryl_subdirs = glob.glob(f'{beryl_dir}/*')

  dest_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/FASTQ"

  beryl_subdirs_starrseq = defaultdict(list)
  for subdir in beryl_subdirs:
    # STARRseq data is organized by species designated as 4 letter string
    # the first letter is the first letter of the genus, capitalized
    # the last 3 letters are the first 3 letters of the species, lower case
    # get these species subdirecotries
    subdir_name = subdir.split('/')[-1]
    if len(subdir_name) == 4:
      print(subdir)
      beryl_subdirs_starrseq[subdir_name] = subdir


  for species in beryl_subdirs_starrseq:
    if species in species_list:
      fastq = glob.glob(f'{beryl_subdirs_starrseq[species]}/FASTQ/*.fastq.gz')
      print(species)
      print(f'{beryl_subdirs_starrseq[species]}/FASTQ')
      for old in fastq:
        new = old.split('/')[-1]
        new = new.replace("2611__", "")
        new = new.replace("__CrossSpp_STARR_pool_19May23", "")
        new = new.replace("2597__", "")

        print(old)
        print(new)

        os.makedirs(f'{dest_dir}/{species}', exist_ok=True)
        os.system(f'cp {old} {dest_dir}/{species}/{new}')

if __name__ == '__main__':
  main()
