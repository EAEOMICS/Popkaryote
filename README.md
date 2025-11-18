# Popkaryote
Population Genomics Pipeline for Bacteria
<img width="1100" height="304" alt="imagen" src="https://github.com/user-attachments/assets/c3aff12b-01d0-4baa-9af2-dde71ed792fc" />

# POPKARYOTE: A novel tool to infer population genomics statistics from prokaryotic organisms

---

## Usage Parameters

| Parameter | Description |
|----------|-------------|
| **vcf** | Path to input VCF file *(optional)* |
| **output** | Directory to save results |
| **max_clusters** | Maximum number of clusters to test *(default = 8)* |
| **window** | Window size for sliding-window statistics *(default = 10000)* |
| **step** | Distance between start positions of windows. Defaults to window size (non-overlapping) *(default = None)* |
| **min_samples_per_pop** | Minimum number of samples per population to compute statistics *(default = 3)* |
| **threads** | Number of threads for running Snippy *(default = 4)* |
| **reference** | Reference genome in FASTA or GENBANK format (can be multi-contig) |
| **gbff** | Path to GBFF annotation file of the reference genome |
| **sequence_type** | Type of sequence file for Snippy:<br>• `fastq` (paired-end)<br>• `fastq-single` (single-end)<br>• `fasta` *(default = fastq)* |
| **sequence_file_location** | Location of sequence files for VCF analysis |
| **dates_file** | TSV file with sample name (column 1) and decimal sampling date (column 2) |
| **unknown_month** | Use default for unknown sampling month *(default = True)* |
| **assemblies** | Location of assembly files to be annotated |
| **type_organism** | Type of organism:<br>• `0` = Bacteria<br>• `1` = Virus *(default = 0)* |

---

## Usage Example

```bash
snakemake --cores {NUM_OF_CORES} --sdm conda \
  --config output='{output_directory_full_path}' \
          step=100 \
          min_samples_per_pop=30 \
          threads={NUM_OF_CORES} \
          reference='{REF_GENOME_FULL_PATH}'

