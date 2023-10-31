#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
from cyvcf2 import VCF, Writer
import numpy as np

def main():

  species = sys.argv[1]
  d_in = '/Genomics/kocherlab/bjarnold/STARRseq/data/VCFs/alternate_reference/00_filtered'
  d_out = '/Genomics/kocherlab/bjarnold/STARRseq/data/VCFs/alternate_reference/01_filtered_within_sample_allele_frequency'
  vcf_in_file = f'{d_in}/{species}.vcf.gz'
  vcf_out_file = f'{d_out}/{species}.vcf.gz'

  vcf_in = VCF(vcf_in_file, gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles
  vcf_out = Writer(vcf_out_file, vcf_in)
  
  for i,variant in enumerate(vcf_in):
    if np.sum(variant.gt_alt_freqs > 0.5) >= 2:
      vcf_out.write_record(variant)
    # if i <= 15:
    #   print(variant.gt_alt_depths)
    #   print(variant.gt_ref_depths)
    #   print(variant.gt_alt_freqs)
    #   print(np.sum(variant.gt_alt_freqs > 0.5))
    #   print("####")

  vcf_in.close()
  vcf_out.close()

  cmd = f"tabix {vcf_out_file}"
  os.system(cmd)

if __name__ == '__main__':
  main()
