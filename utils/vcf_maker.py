import os
import subprocess
import re

def build_tab_with_bash(directory, output_file):
    # Ajusta el patrón según la convención de tus archivos
    bash_command = f"""
    find "{directory}" -type f -name "*1.fastq*" | while read R1; do
    R2=$(echo "$R1" | sed 's/1.fastq/2.fastq/')
    if [ -f "$R2" ]; then
        SAMPLE=$(basename "$R1" | sed 's/_1\.fastq.*$//')
        echo -e "$SAMPLE\t$R1\t$R2"
    fi
    done > "{output_file}"
    """

    subprocess.run(bash_command, shell=True, executable='/bin/bash')


def build_tab_with_bash_fasta(directory, output_file):
    # Para archivos FASTA: un solo archivo por muestra
    bash_command = f"""
    find {directory} -type f -name "*.f*" | while read FASTA; do
        SAMPLE=$(basename $FASTA | sed 's/\\.f.*//')
        echo -e \"$SAMPLE\\t$FASTA\"
    done > {output_file}
    """
    subprocess.run(bash_command, shell=True, executable='/bin/bash')


def vcf_maker(fastq_location, output_dir,reference,threads=4,file_type='fastq'):
    print("\n[Pre-PHASE VCF] Computing VCF file with Snippy...")
    os.makedirs(output_dir,exist_ok=True)
	
    output_file=output_dir+'/snippy_input.tab'

    if file_type=='fastq':
      build_tab_with_bash(directory=fastq_location, output_file=output_file)
    else:
      build_tab_with_bash_fasta(directory=fastq_location, output_file=output_file)
    
    snippy_multi_command = f"snippy-multi {output_file} --ref {reference} --cpus {threads} > {os.path.join(output_dir, 'run_snippy.sh')}"
    snippy_change_bash= f"sed -i \"s|--outdir '|--outdir '{output_dir}/|g\" {os.path.join(output_dir, 'run_snippy.sh')}"
    snippy_run=f"bash {os.path.join(output_dir, 'run_snippy.sh')}"
    snippy_change_bash2=f"sed -i \"s|snippy-core |cd {output_dir}; snippy-core |g\" {os.path.join(output_dir, 'run_snippy.sh')}; cd -"
    vcf_merge=f"bcftools merge {output_dir}/*/snps.vcf.gz -o {output_dir}/merged.vcf"
    subprocess.run(snippy_multi_command, shell=True, executable='/bin/bash')
    subprocess.run(snippy_change_bash, shell=True, executable='/bin/bash')
    subprocess.run(snippy_change_bash2, shell=True, executable='/bin/bash')
    subprocess.run(snippy_run, shell=True, executable='/bin/bash')
    subprocess.run(vcf_merge, shell=True, executable='/bin/bash')

    reindexing=f"""
    echo "pairwise VCF merging started: $(date)"
    ls {output_dir}/*/snps.vcf.gz > {output_dir}/vcf_list.txt
    ls {output_dir}/*/snps.bam > {output_dir}/bam_list.txt
    [ $(ulimit -Sn) -lt $(ls {output_dir}/*/snps.vcf.gz | wc -l) ] && ulimit -Sn $(( $(ls {output_dir}/*/snps.vcf.gz | wc -l) + 100 ))
    bcftools merge -l {output_dir}/vcf_list.txt -Oz --threads {threads} -o {output_dir}/merged_sites.vcf.gz
    
    bcftools index {output_dir}/merged_sites.vcf.gz
    bcftools query -f'%CHROM\t%POS0\t%POS\n' {output_dir}/merged_sites.vcf.gz > {output_dir}/sites.bed
    """
    subprocess.run(reindexing, shell=True, executable='/bin/bash')