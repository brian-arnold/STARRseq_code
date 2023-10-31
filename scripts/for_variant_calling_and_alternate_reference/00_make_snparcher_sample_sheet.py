#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import glob

def main():

  species_for_variant_calling = ['Nmel', 'Aaur', 'Apur', 'Avir']
  fastq_dir = "/Genomics/kocherlab/bjarnold/STARRseq/data/FASTQ"

  new_assembly_dir="/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v3.1"
  old_assembly_dir="/Genomics/kocherlab/lab/data/GenomeReleases/official_release_v2.1.1"

  dest_assembly_dir="/Genomics/kocherlab/bjarnold/STARRseq/code/snakemake_snparcher/ref_genomes_uncompressed"

  # GET READ LOCATIONS
  reads = defaultdict(lambda: defaultdict(list))
  for datatype in ['input', 'RNA']:
    for s in species_for_variant_calling:
      read1 = glob.glob(f'{fastq_dir}/{s}/*{datatype}-read-1*')
      read3 = glob.glob(f'{fastq_dir}/{s}/*{datatype}-read-3*')
      # sort files in read1 and read3
      read1 = sorted(read1)
      read3 = sorted(read3)
      for i,j in zip(read1, read3):
        assert i.replace("read-1", "read-3") == j
        reads[datatype][s].append((i,j))

  # for s in reads:
  #   print(s, reads[s])



  # GET REFERENCE GENOME LOCATIONS
  new_assemblies = defaultdict(str)
  # get list of species, which correspond to the subirectories in the new assembly directory
  subdirs = glob.glob(f'{new_assembly_dir}/*')
  species = [s.split('/')[-1] for s in subdirs]
  for s in species:
    # s is a 4 character string, but convert the last 3 characters to lower case
    s_lower = s[:-3] + s[-3:].lower()
    if s_lower in species_for_variant_calling:
      ref_fastas_paths = glob.glob(f'{new_assembly_dir}/{s}/*.fasta.gz')
      ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
      ref_fasta = [f for f in ref_fastas if '.alt.' not in f][0]

      # ref = f'{s_lower}\t{new_assembly_dir}/{s}/{ref_fasta}'

      # check if {dest_assembly_dir}/{ref_fasta}.gz exists
      if not os.path.isfile(f'{dest_assembly_dir}/{ref_fasta}'):
        ref = f'{new_assembly_dir}/{s}/{ref_fasta}'
        os.system(f'cp {ref} {dest_assembly_dir}')
        os.system(f'gunzip {dest_assembly_dir}/{ref_fasta}')

      new_assemblies[s_lower] = f'{dest_assembly_dir}/{ref_fasta}'.replace('.gz','')

  for s in new_assemblies:
    print(s, new_assemblies[s])

  print('###')
  old_assemblies = defaultdict(str)
  subdirs = glob.glob(f'{old_assembly_dir}/*')
  species = [s.split('/')[-1] for s in subdirs]
  for s in species:
    # s is a 4 character string, but convert the last 3 characters to lower case
    s_lower = s[:-3] + s[-3:].lower()
    if s_lower in species_for_variant_calling and s_lower not in new_assemblies:
      ref_fastas_paths = glob.glob(f'{old_assembly_dir}/{s}/*_v2.1.0.fasta.gz')
      ref_fastas = [f.split('/')[-1] for f in ref_fastas_paths]
      assert len(ref_fastas) == 1, print(s)
      ref_fasta = ref_fastas[0]

      if not os.path.isfile(f'{dest_assembly_dir}/{ref_fasta}'):
        ref = f'{old_assembly_dir}/{s}/{ref_fasta}'
        os.system(f'cp {ref} {dest_assembly_dir}')
        os.system(f'gunzip {dest_assembly_dir}/{ref_fasta}')
        old_assemblies[s_lower] = f'{dest_assembly_dir}/{ref_fasta}'.replace('.gz','')
  
  for s in old_assemblies:
    print(s, old_assemblies[s])


  print('BioSample,LibraryName,Run,refGenome,refPath,fq1,fq2')
  run_count = 0
  for datatype in reads:
    for s in reads[datatype]:
      for i,r in enumerate(reads[datatype][s]):
        run_count += 1
        biosamp = f'{s}_{datatype}_F{i+1}'
        libname = biosamp
        run = run_count
        ref = s
        if s in new_assemblies:
          refpath = new_assemblies[s]
        elif s in old_assemblies:
          refpath = old_assemblies[s]
        else:
          sys.exit('ERROR')
        fq1 = r[0]
        fq2 = r[1]
        print(','.join([biosamp, libname, str(run), ref, refpath, fq1, fq2]))
      



if __name__ == '__main__':
  main()
