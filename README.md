# STARRseq_code

Calling peaks from DNA input, RNA output FASTQ files

Before beginning
- In `STARRseq/code/scripts/`
    - Use `bwa_index.py` to create indices of all the genomes you’ll use for peak calling
    - Use `make_repeatmasker_BED.py` to get BED files from repeat modeler output to potentially exclude repeat regions
    - Use `collect_fastqs.py` to copy fastq files for each species into a single location. This script also renames files that have unnecessary strings

Fastq to peak data
- In `STARRseq/code/scripts/for_peak_calling/`
    - Use `compute_effective_genome_size.py` to compute the number of AGCT (non N’s) for effective genome size, this info will get used in the config file for the snakelike workflow. You can use the `reference_genome_locations.py` script to locate the appropriate reference genome to use
    - After you run the snakelike workflow (below), use the `reorganize_files.py` script to copy peak-calling files for all species to a single directory that you can then copy to your local machine for running analyses there
- In ‘/STARRseq/code/snakemake_peak_calling’
    - Put all output in `STARRseq/data/peak_calling_snakemake_output`
    - Run snakelike with `run.sh` script, making sure all relevant parent directories are mounted so the singularity container has access to necessary files


Calling short germline variants, from Fastq to VCF

- In `STARRseq/code/scripts/for_variant_calling_and_alternate_reference/`
    - Use `00_make_snparcher_sample_sheet.py` to help make a sample.csv file. NOTE: this also makes a copy of the reference genome that is not compressed, a requirement for snpArcher. These are stored in `ref_genomes_uncompressed` in `snakemake_snparcher`
    - After running snakelike pipeline (below)
        - Run `01_compute_depth_per_sample.py` to compute the depths for each species as a sanity check
        - Run `02_filter_vcfs.py` to filter VCFs for minimum depth, missing genotypes, and balletic sites
        - Run `03_filter_vcfs_within_sample_allele_frequency.py` to create another VCF file, filtered for sites in which 2+ bores have allele frequencies greater than 50%
        - Run `04_make_alternate_reference.py` to create a new reference sequence using the VCF file from the previous step
        - Run `05_rename_scaffolds.py` to revert scaffold names to the original, since these get renamed from the previous step
- In ‘STARRseq/code/snakemake_snparcher’
    - You changed HaplotypeCaller command so that it ignores duplicates flag in BAM
    - In `config/config.yaml`, keep db_scatter_factor set to 1; scattering DB intervals led to an error in which SNPs were duplicated such that merging them back together caused a GATK error
    - Raw VCF files will be stored in the `results` subdirectory


Analyses

- use notebooks/02_compare_peak_callers_across_bioreps_REFACTORED.ipynb to analyze peak data