#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

def main():

  # use lowercase for species to match syntax of naming convention for all other files/directories

  rm_dirs = {'Aaur':	'/scratch/tmp/aewebb/Annotations/AAUR/AAUR_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Apur':	'/scratch/tmp/aewebb/Annotations/APUR/APUR_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Avir':	'/scratch/tmp/aewebb/Annotations/AVIR/AVIR_v3/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Hlig':	'/scratch/tmp/aewebb/Annotations/HLIG/HLIG_v3/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Lalb':	'/scratch/tmp/aewebb/Annotations/LALB/LALB_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Lbal':	'/scratch/tmp/aewebb/Annotations/LBAL/LBAL_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Lvie':	'/scratch/tmp/aewebb/Annotations/LVIE/LVIE_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  # 'LVIE':	'/scratch/tmp/aewebb/Annotations/LVIE/LVIE_v3.alt.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Lzep':	'/scratch/tmp/aewebb/Annotations/LZEP/LZEP_v3.1/Assembly/RepeatModeler/R1/MaskedAssembly/',
  'Hqua':	'/Genomics/kocherlab/aewebb/RM/HQUA/',
  'Nmel':	'/Genomics/kocherlab/aewebb/RM/NMEL/' }

  out_dir = '/Genomics/kocherlab/bjarnold/STARRseq/data/repeat_modeler_intervals'

  for i,species in enumerate(rm_dirs):
    # if i > 0:
    #   continue
    rm_dir = rm_dirs[species]
    out_files = glob.glob(rm_dir + '*.out')
    assert len(out_files) == 1
    in_file = out_files[0]
    command = f'rmsk2bed < {in_file} > {out_dir}/{species}_repeat_elements.bed'
    os.system(command)


if __name__ == '__main__':
  main()
