#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import subprocess

def main():

  # GET ASSEMBLY LOCATIONS FROM reference_genome_locations.py

  genomes = ['/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/AVIR/AVIR_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/AAUR/AAUR_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/APUR/APUR_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v2.1.1/NMEL/NMEL_genome_v2.1.0.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/HLIG/HLIG_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/LZEP/LZEP_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/LBAL/LBAL_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/LALB/LALB_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1/LVIE/LVIE_genome_v3.fasta.gz',
  '/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v2.1.1/HQUA/HQUA_genome_v2.1.0.fasta.gz',
  '/Genomics/kocherlab/bjarnold/STARRseq/data/genomes/AMEL/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz',
  '/Genomics/kocherlab/bjarnold/STARRseq/data/genomes/BIMP/BIMP_genome_vNP3.fna.gz',
  '/Genomics/kocherlab/bjarnold/STARRseq/data/genomes/SINV/GCF_016802725.1_UNIL_Sinv_3.0_genomic.fna.gz'
  ]

  for g in genomes:
    print(g)
    cmd = ['zcat', g, '|', 'grep', '-v', '"^>"', '|', 'tr', '-cd', '"ACGTacgt"', '|', 'wc', '-c']
    os.system(" ".join(cmd))
    print("#### including N's ####")
    cmd = ['zcat', g, '|', 'grep', '-v', '"^>"', '|', 'tr', '-cd', '"ACGTNacgtn"', '|', 'wc', '-c']
    os.system(" ".join(cmd))
    print("\n")

    # result = subprocess.run(cmd, capture_output=True, text=False)
    # print(result)
    # print("#### including N's ####")
    # cmd = ['zcat', g, '|', 'grep', '-v', '"^>"', '|', 'tr', '-cd', '"ACGTNacgtn"', '|', 'wc', '-c']
    # result = subprocess.run(cmd, capture_output=True, text=False)
    # print(result)   
    # print("\n\n\n\n\n\n")


if __name__ == '__main__':
  main()
