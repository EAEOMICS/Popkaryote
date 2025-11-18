# POPKARYOTE: A novel tool to infer population genomics statistics from prokaryotic organisms

<img width="1100" height="304" alt="imagen" src="https://github.com/user-attachments/assets/c3aff12b-01d0-4baa-9af2-dde71ed792fc" />

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
          reference='{REF_GENOME_FULL_PATH}' \
          sequence_file_location='{SEQUENCE_DIRECTORY_FULL_PATH}'
```

---

## 📊 Pipeline Phases

| Phase | Description |
|-------|-------------|
| **PHASE_1_Automatic_population_inference** | Inference of populations based on FASTA/FASTQ files provided by the user |
| **PHASE_1_Automatic_population_inference_with_provided_VCF** | Inference of populations based on VCF file provided by the user |
| **PHASE_2_Global_statistics** | Inference of nucleotide diversity, Watterson’s Theta, Tajima’s D, Linkage Disequilibrium, Maximum Likelihood Tree, Recombination detection, Ancestral Reconstruction |
| **PHASE_3_Population_specific_statistics** | Specific population statistics for every inferred population (nucleotide diversity, Watterson’s Theta, Tajima’s D, Linkage Disequilibrium) |
| **PHASE_4_Pairwise_FST_comparisons_using_Hudson_method** | FST calculation for every population using Hudson’s method |
| **PHASE_5_Plot_population_comparisons** | Generation of plots for comparisons of population statistics results |
| **PHASE_6_Annotation** | Annotation of all the provided genomes |
| **PHASE_7_Pangenome** | Estimation of the pangenome |
| **PHASE_8_PhaME** | Execution of PhaME module |
| **PHASE_9_AMRfinder** | Inference of Antibiotic Resistance Genes (ARGs) and Virulence Factors (VFs) |
| **PHASE_10_SNPEff** | Annotation of the effects of SNPs in the genes |

---

## Usage Example for computing Pangenome

```bash
snakemake PHASE_7_Pangenome --cores 10 --sdm conda \
  --config output='{OUTPUT_DIR}' step=100 min_samples_per_pop=2 \
  threads={NUM_OF_CORES} reference='{REF_GENOME_FULL_PATH}' \
  sequence_type='fasta' sequence_file_location='{SEQUENCE_DIRECTORY_FULL_PATH}'\
  assemblies='{ASSEMBLIES_DIRECTORY_FULL_PATH}'

```

POPKARYOTE is a pipeline developed entirely by the EGG group of the Universitat Autonoma de Barcelona, all rights reserved" +"\n"+"https://github.com/EAEOMICS
