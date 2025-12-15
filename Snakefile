### POPKARYOTE

# ANSI escape codes para verde y reset
GREEN = "\033[1;32m"
RESET = "\033[0m"

banner = r"""
            ++     +                             
     ++     ++++++++    ++                      
       ++ +++       ++++                        
        ++   +++++++   ++                       
       ++  ++       ++  ++                                               
  ++   +  ++         ++  ++               _ (`-.                  _ (`-.  .-. .-')   ('-.       _  .-')                              .-') _       ('-. _ 
     +++  +   ++ ++   +  ++++++          ( (OO  )                ( (OO  ) \  ( OO )  ( OO ).-. ( \( OO )                            (  OO) )    _(  OO) )
      ++  +  ++   +   ++  +             _.`     \  .-'),-----.  _.`     \ ,--. ,--.  / . --. /  ,------.    ,--.   ,--. .-'),-----. /     '._  (,------. 
       +  +   ++++     +  ++           (__...--'' ( OO' .-.  '(__...--''  |  .'   /  | \-.  \   |   /`. '    \  `.'  / ( OO'  .-.  '|'--...__)  |  .---' 
       ++  +           ++  ++*+++      |  /  | | /   |  | |  | |  /  | |  |      /   '-'  |  |  |  /  | |  .-')     /  /   |  | |  |'--.  .--' (|  |/
        +   +           +++  +         |  |_.' | \_) |  |\|  | |  |_.' |  |     '    | |_.'  |  |  |_.' | (OO  \   /   \_) |  |\|  |   |  |    (|  '--.
   ++++  +  ++            +  ++        |  .___.'   \ |  | |  | |  .___.'  |  .   \   |  .-.  |  |  .  '.'  |   /  /\_    \ |  | |  |   |  |     |  .--'
          +  ++            +  ++       |  |         `'  '-'  ' |  |       |  |\   \  |  | |  |  |  |\  \   `-./  /.__)    `'  '-'  '   |  |     |  `---.
           ++   +          ++  ++++    `--'           `-----'  `--'       `--' '--'  `--' `--'  `--' '--'    `--'           `-----'    `--'     `------'
            +++  ++       ++  ++           
         +++  ++   +++++++   ++            
               +++++     +++++            
            ++    +  ++++     ++          
                      +         +                      
  """
print(f"{GREEN}{banner}{RESET}")


configfile: "utils/config.yml"
OUTPUT_DIR = config["output"]
WINDOW_SIZE = config['window']
MIN_SAMPLES = config['min_samples_per_pop']
MAX_CLUSTERS = config['max_clusters']
STEP= config['step']
FASTA= config['assemblies']
TYPE=config['type_organism']
if config["vcf"] != "":
  vcf = config["vcf"]
else:
  vcf = f"{OUTPUT_DIR}/output.vcf.gz"

if config['unknown_month']!= "False":
  UM='True'
else: UM='False'

if config['dates_file']!= "":
  DF=config['dates_file']
else: DF=None


#Define rules
rule all:
  input:
    vcf,
    f"{OUTPUT_DIR}/global",
    f"{OUTPUT_DIR}/by_population",
    f"{OUTPUT_DIR}/fst_pairwise",
    f"{OUTPUT_DIR}/LDWeaver_Results",
    f"{OUTPUT_DIR}/comparison_plots",
    f"{OUTPUT_DIR}/annotation",
    f"{OUTPUT_DIR}/Pangenome",
    f"{OUTPUT_DIR}/PhaME",
    f"{OUTPUT_DIR}/PhaME/files",
    f"{OUTPUT_DIR}/Pangenome/coevolution_results",
    f"{OUTPUT_DIR}/comparison_plots/pi_comparison.png",
    f"{OUTPUT_DIR}/comparison_plots/tajima_comparison.png",
    f"{OUTPUT_DIR}/comparison_plots/WT_comparison.png",
    f"{OUTPUT_DIR}/explained_variance.npy",
    f"{OUTPUT_DIR}/pca_coords.csv",
    f"{OUTPUT_DIR}/pca_clusters.png",
    f"{OUTPUT_DIR}/auto_populations.csv",
    f"{OUTPUT_DIR}/elbow_clusters.png",
    f"{OUTPUT_DIR}/global/pi_window_global.png",
    f"{OUTPUT_DIR}/global/pi_raw_global.png",
    f"{OUTPUT_DIR}/global/pi_global.png",
    f"{OUTPUT_DIR}/global/pi_window_global_barplot.png",
    f"{OUTPUT_DIR}/global/pi_raw_global_barplot.png",
    f"{OUTPUT_DIR}/global/pi_global_barplot.png",
    f"{OUTPUT_DIR}/global/pi_global.csv",
    f"{OUTPUT_DIR}/global/tajima_d_global.csv",
    f"{OUTPUT_DIR}/global/tajima_raw_global.png",
    f"{OUTPUT_DIR}/global/tajima_global.png",
#    f"{OUTPUT_DIR}/global/LD_heatmap.png",
    f"{OUTPUT_DIR}/global/hd_global.txt",
    f"{OUTPUT_DIR}/global/Hd_50_100_200.png",
    f"{OUTPUT_DIR}/global/WT_raw_global.png",
    f"{OUTPUT_DIR}/global/WT_global.png",
    f"{OUTPUT_DIR}/global/iqtree_output.contree",
    f"{OUTPUT_DIR}/core.ref.fa",
    f"{OUTPUT_DIR}/core.ref.fa.fai",
    f"{OUTPUT_DIR}/bam_list.txt",
    f"{OUTPUT_DIR}/sites.bed",
    f"{OUTPUT_DIR}/output.vcf.gz.csi",
    f"{OUTPUT_DIR}/results/windows.bed",
    f"{OUTPUT_DIR}/output.vcf",
    f"{OUTPUT_DIR}/PhaME/phame.ctl",

rule help:
  run:
    print("\n" + "="*100)
    print("POPKARYOTE: A novel tool to infer population genomics statistics from prokaryotic organisms".center(80))
    print("="*100 + "\n")

    print("Usage Parameters:\n")

    print("  vcf                     Path to input VCF file (optional)")
    print("  output                  Directory to save results")
    print("  max_clusters            Maximum number of clusters to test (default = 8)")
    print("  window                  Window size for sliding-window statistics (default = 10000)")
    print("  step                    Distance between start positions of windows")
    print("                          If not given, defaults to window size (non-overlapping) (default = None)")
    print("  min_samples_per_pop     Minimum number of samples per population to compute statistics (default = 3)")
    print("  threads                 Number of threads for running Snippy (default = 4)")
    print("  reference               Reference genome in FASTA or GENBANK format (can be multi-contig)")
    print("  gbff                    Path to gbff annotation file of the refernce genome:")
    print("  sequence_type           Type of sequence file for Snippy:")
    print("                            - fastq (paired-end)")
    print("                            - fastq-single (single-end)")
    print("                            - fasta (default = fastq)")
    print("  sequence_file_location  Location of sequence files for VCF analysis")
    print("  dates_file              TSV file with sample name (column 1) and decimal sampling date (column 2)")
    print("  unknown_month           Use default for unknown sampling month (default = True)")
    print("  assemblies              Location of assembly files to be annotated")
    print("  type_organism           Type of organism:")
    print("                            - 0 = Bacteria")
    print("                            - 1 = Virus (default = 0)")

    print("\n" + "="*100 + "\n")
    print("usage example:")
    print("snakemake --cores {NUM_OF_CORES} --sdm conda --config output='{output_directory_full_path}' step=100 min_samples_per_pop=30 threads={NUM_OF_CORES} reference='{REF_GENOME_FULL_PATH}' ")
    print("\n" + "="*100 + "\n")
    print("POPKARYOTE is a pipeline developed entirely by the EGG group of the Universitat Autonoma de Barcelona, all rights reserved" +"\n"+"https://github.com/EAEOMICS")
    print("\n" + "="*100 + "\n")

    
rule vcf_from_sequences:
  input:
    seq=config['sequence_file_location'],
    ref=config['reference'],
  output:
    ref=f"{OUTPUT_DIR}/core.ref.fa",
    fai=f"{OUTPUT_DIR}/core.ref.fa.fai",
    bamlist=f"{OUTPUT_DIR}/bam_list.txt",
    sites=f"{OUTPUT_DIR}/sites.bed",
  conda:
    "envs/vcfmaker.yml"
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/vcf_from_sequences.log"
  shell:
    """
    python -c "from utils.vcf_maker import vcf_maker; \
      vcf_maker(fastq_location='{input.seq}', output_dir='{OUTPUT_DIR}', \
      reference='{input.ref}', threads={threads}, \
      file_type='{config[sequence_type]}')"
    cp "$(find '{OUTPUT_DIR}' -type f -name 'ref.fa' | head -n1)" '{output.ref}'
    cp "$(find '{OUTPUT_DIR}' -type f -name 'ref.fa.fai' | head -n1)" '{output.fai}'
    """

rule variant_recalling:
  input:
    ref=f"{OUTPUT_DIR}/core.ref.fa",
    fai=f"{OUTPUT_DIR}/core.ref.fa.fai",
    bamlist=f"{OUTPUT_DIR}/bam_list.txt",
    sites=f"{OUTPUT_DIR}/sites.bed",
  output:
    vcf=vcf,
    tbi=f"{OUTPUT_DIR}/output.vcf.gz.csi",
    windows=f"{OUTPUT_DIR}/results/windows.bed",
    o4=f"{OUTPUT_DIR}/output.vcf",
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/variant_recalling.log"
  shell:
    """
    [ $(ulimit -Sn) -lt $(ls {OUTPUT_DIR}/*/snps.vcf.gz | wc -l) ] && ulimit -Sn $(( $(ls {OUTPUT_DIR}/*/snps.vcf.gz | wc -l) + 100 ))
    echo -e "\e[32mVariant re-calling Started: $(date)\e[0m"
    bedtools makewindows -g {input.fai} -w 100000 > {output.windows}
    awk '{{print $1":"$2+1"-"$3}}' {output.windows} | \
    parallel -j {threads} "
      bcftools mpileup -f {input.ref} -b {input.bamlist} -Ou -r {{}} |
      bcftools call -mv -Oz -o {OUTPUT_DIR}/chunk_{{#}}.vcf.gz
    "
    ls {OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | xargs bcftools concat -Oz | bcftools sort -Oz -o {output.vcf}
    bcftools index {output.vcf}
    gzip -d -k {output.vcf}
    echo ""
    echo -e "\e[32mVariant re-calling Ended: $(date) Output file: output.vcf.gz\e[0m"
    """


rule PHASE_1_Automatic_population_inference:
  input:
    i1=f"{OUTPUT_DIR}/output.vcf",
  output:
    o1=f"{OUTPUT_DIR}/explained_variance.npy",
    o2=f"{OUTPUT_DIR}/pca_coords.csv",
    o3=f"{OUTPUT_DIR}/pca_clusters.png",
    o4=f"{OUTPUT_DIR}/auto_populations.csv",
    o5=f"{OUTPUT_DIR}/elbow_clusters.png",
#  conda:
#    "envs/popkaryote.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_1_Automatic_population_inference.log"
  shell:
    """
    python -c "from utils.pca_and_coords import pca_and_save_coords; \
      pca_and_save_coords(vcf_path = '{input}', output_dir = '{OUTPUT_DIR}'); \
      from utils.kmeans_and_plot import cluster_and_plot; \
      cluster_and_plot(output_dir = '{OUTPUT_DIR}', \
      coords_csv = '{output.o2}',\
      var_exp_path = '{output.o1}',\
      max_clusters = {MAX_CLUSTERS},\
      min_cluster_size = {MIN_SAMPLES})"
    """
rule PHASE_1_Automatic_population_inference_with_provided_VCF:
  input:
    i1=config["vcf"],
  output:
    o1=f"{OUTPUT_DIR}/explained_variance.npy",
    o2=f"{OUTPUT_DIR}/pca_coords.csv",
    o3=f"{OUTPUT_DIR}/pca_clusters.png",
    o4=f"{OUTPUT_DIR}/auto_populations.csv",
    o5=f"{OUTPUT_DIR}/elbow_clusters.png",
#  conda:
#    "envs/popkaryote.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_1_Automatic_population_inference_with_provided_VCF.log"
  shell:
    """
    python -c "from utils.pca_and_coords import pca_and_save_coords; \
      pca_and_save_coords(vcf_path = '{input}', output_dir = '{OUTPUT_DIR}'); \
      from utils.kmeans_and_plot import cluster_and_plot; \
      cluster_and_plot(output_dir = '{OUTPUT_DIR}', \
      coords_csv = '{output.o2}',\
      var_exp_path = '{output.o2}',\
      max_clusters = {MAX_CLUSTERS},\
      min_cluster_size = {MIN_SAMPLES})"
    """
rule PHASE_2_Global_statistics:
  input:
    i1=vcf,
    i2=config['reference'],
  output:
    o1=directory(f"{OUTPUT_DIR}/global"),
    o2=f"{OUTPUT_DIR}/global/pi_window_global.png",
    o3=f"{OUTPUT_DIR}/global/pi_raw_global.png",
    o4=f"{OUTPUT_DIR}/global/pi_global.png",
    o5=f"{OUTPUT_DIR}/global/pi_window_global_barplot.png",
    o6=f"{OUTPUT_DIR}/global/pi_raw_global_barplot.png",
    o7=f"{OUTPUT_DIR}/global/pi_global_barplot.png",
    o9=f"{OUTPUT_DIR}/global/pi_global.csv",
    o10=f"{OUTPUT_DIR}/global/tajima_d_global.csv",
    o11=f"{OUTPUT_DIR}/global/tajima_raw_global.png",
    o12=f"{OUTPUT_DIR}/global/tajima_global.png",
#    o13=f"{OUTPUT_DIR}/global/LD_heatmap.png",
    o14=f"{OUTPUT_DIR}/global/hd_global.txt",
    o15=f"{OUTPUT_DIR}/global/Hd_50_100_200.png",
    o16=f"{OUTPUT_DIR}/global/Watterson_theta_global.csv",
    o17=f"{OUTPUT_DIR}/global/WT_raw_global.png",
    o18=f"{OUTPUT_DIR}/global/WT_global.png",
    o19=directory(f"{OUTPUT_DIR}/LDWeaver_Results"),
    o20=f"{OUTPUT_DIR}/global/iqtree_output.contree",
#  conda:
#    "envs/popkaryote.yml"
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/PHASE_2_Global_statistics.log"
  shell:
    """
    python -c "from utils.phase2_global_statistics import compute_global_statistics; \
      compute_global_statistics(\
      vcf_path='{input.i1}',\
      output_dir='{output.o1}',\
      window_size={WINDOW_SIZE},\
      step={STEP},\
      threads={threads},\
      dates='{DF}',\
      month={UM},\
      reference='{input.i1}')"
    """
        
rule PHASE_3_Population_specific_statistics:
  input:
    i1=vcf,
    i2=f"{OUTPUT_DIR}/auto_populations.csv",
  output:
    o1=directory(f"{OUTPUT_DIR}/by_population"),
#  conda:
#    "envs/popkaryote.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_3_Population_specific_statistics.log"
  shell:
    """
    python -c "from utils.phase3_population_statistics import compute_population_statistics; \
      compute_population_statistics( \
      vcf_path='{input.i1}', \
      populations_csv='{input.i2}', \
      output_dir='{output.o1}', \
      window_size={WINDOW_SIZE},)"
    """

rule PHASE_4_Pairwise_FST_comparisons_using_Hudson_method:
  input:
    i1=vcf,
    i2=f"{OUTPUT_DIR}/auto_populations.csv",
  output:
    o1=directory(f"{OUTPUT_DIR}/fst_pairwise"),
#  conda:
#    "envs/popkaryote.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_4_Pairwise_FST_comparisons_using_Hudson_method.log"
  shell:
    """
    python -c "from utils.phase4_pairwise_fst import compute_pairwise_fst; \
      compute_pairwise_fst( \
      vcf_path='{input.i1}', \
      populations_csv='{input.i2}', \
      output_dir='{output.o1}', \
      window_size={WINDOW_SIZE}, \
      min_samples={MIN_SAMPLES})"
    """
    
rule PHASE_5_Plot_population_comparisons:
  input:
    i1=f"{OUTPUT_DIR}/by_population",
  output:
    o1=directory(f"{OUTPUT_DIR}/comparison_plots"),
    o2=f"{OUTPUT_DIR}/comparison_plots/pi_comparison.png",
    o3=f"{OUTPUT_DIR}/comparison_plots/tajima_comparison.png",
    o4=f"{OUTPUT_DIR}/comparison_plots/WT_comparison.png",
#  conda:
#    "envs/popkaryote.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_5_Plot_population_comparisons.log"
  shell:
    """
    python -c "from utils.plot_global_comparisons import plot_population_comparisons; \
      plot_population_comparisons( \
      input_dir='{input.i1}', \
      output_dir='{output.o1}')"
    """

rule PHASE_6_Annotation:
  input:
    i1=config['assemblies']
  output:
    o1=directory(f"{OUTPUT_DIR}/annotation"),
  conda:
    "envs/annotation.yml"
  threads: 4
  log:
    f"{OUTPUT_DIR}/log/PHASE_6_Annotation.log"
  shell:
    """
    mkdir {output.o1}
    for file in {input.i1}/*.f*; do
    basename=$(basename "$file")
    basename="${{basename%%.f*}}"
    cleaned="{output.o1}/$basename.cleaned.fasta"
    awk '/^>/{{print ">contig_" ++i; next}}{{print}}' "$file" > "$cleaned"  
    prokka --outdir {output.o1}/"$basename" --prefix "$basename" "$cleaned"
    done
    REF=$(find {OUTPUT_DIR} -type f -name "ref.fa" -print -quit)
    prokka --outdir {output.o1}/Reference --prefix "Reference"  "$REF"
    """

rule PHASE_7_Pangenome:
  input:
    i1=f"{OUTPUT_DIR}/annotation",
    i2=f"{OUTPUT_DIR}/global/iqtree_output.contree"
  output:
    o1=directory(f"{OUTPUT_DIR}/Pangenome"),
    o2=directory(f"{OUTPUT_DIR}/Pangenome/coevolution_results"),
  conda:
    "envs/annotation.yml"
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/PHASE_7_Pangenome.log"
  shell:
    """
    mkdir -p {output.o2}
    panaroo -i {input.i1}/*/*.gff -o {output.o1} --clean-mode strict -a core --core_threshold 0.95 -t {threads}
    panaroo-spydrpick -i {output.o1}/gene_presence_absence.Rtab -o {output.o2} --tree {input.i2}
    """
rule PHASE_8_PhaME:
  input:
    i1 = f"{OUTPUT_DIR}/annotation",
    ref = config['reference']
  output:
    o1 = directory(f"{OUTPUT_DIR}/PhaME"),
    o2 = directory(f"{OUTPUT_DIR}/PhaME/files"),
    phame = f"{OUTPUT_DIR}/PhaME/phame.ctl"
  conda:
    "envs/annotation.yml"
  threads: config['threads']
  params:
    fasta_source = FASTA,
    code_type = TYPE
  log:
    f"{OUTPUT_DIR}/log/PHASE_8_PhaME.log"
  shell:
    """
    mkdir -p {output.o2}
    cp {input.i1}/*/*.cleaned.fasta {output.o2} 2>/dev/null || true
    
    if [ -d "{params.fasta_source}" ]; then
      cp {params.fasta_source}/* {output.o2}
    fi
    cp {input.ref} {output.o2}
    REF_NAME=$(basename "{input.ref}")
    cat > {output.phame} <<EOF
    refdir = {output.o2}
    workdir = {output.o1}
    reference = 1
    reffile = $REF_NAME
    project = Phame_alignment
    cdsSNPS = 1
    buildSNPdb = 1
    FirstTime = 1
    data = 3
    reads = 2
    tree = 2
    bootstrap = 1
    N = 10000
    PosSelect = 3
    code = {params.code_type}
    clean = 1
    threads = {threads}
    cutoff = 0.1
    EOF
    
    phame {output.phame}
    """
rule PHASE_9_AMRfinder:
  input:
    i1=config['assemblies'],
  output:
    o1=directory(f"{OUTPUT_DIR}/AMR_and_VF"),
  conda:
    "envs/amrfinder.yml"
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/PHASE_9_AMRfinder.log"
  shell:
    """
    mkdir -p {output.o1}
    amrfinder -u
    for f in {input.i1}/*;do
      basename=$(basename "$f")
      basename="${{basename%%.f*}}"
      amrfinder --nucleotide $f --plus\
      --threads {threads} \
      --name $basename \
      -o {output.o1}/"$basename".tsv
    done
    """

rule PHASE_10_SNPEff:
  input:
    i1=config['assemblies'],
    ref=config['gbff'],
    i2=f"{OUTPUT_DIR}/output.vcf"
  output:
    o1=directory(f"{OUTPUT_DIR}/SNPEff"),
  conda:
    "envs/amrfinder.yml"
  threads: config['threads']
  log:
    f"{OUTPUT_DIR}/log/PHASE_10_SNPEff.log"
  shell:
    """
    mkdir -p {output.o1}
    mkdir -p ./.snakemake/conda/*/share/snpeff-5.3.0a-0/data/SNPeff_db
    cp {input.ref} ./.snakemake/conda/*/share/snpeff-5.3.0a-0/data/SNPeff_db
    mv ./.snakemake/conda/*/share/snpeff-5.3.0a-0/data/SNPeff_db/*.gbff > ./.snakemake/conda/*/share/snpeff-5.3.0a-0/data/SNPEff_db/genes.gbk
    echo "SNPEff_db.genome : SNPEff" >> ./.snakemake/conda/*/share/snpeff-5.3.0a-0/data/snpEff.config 
    java -Xmx1G -jar ./.snakemake/conda/*/share/snpeff-5.3.0a-0/snpEff.jar build -genbank -v SNPEff_db
    java -Xmx40G -jar ./.snakemake/conda/*/share/snpeff-5.3.0a-0/snpEff.jar -v -stats {output.o1}/snpEff.html SNPEff_db {input.i2} > {output.o1}/variants.ann.vcf

    """







