import gzip
import os

def clean_vcf_for_scikit_allel(vcf_path, replace_missing=False):
    """
    Cleans a VCF file to make it compatible with scikit-allel.
    - If replace_missing=True, replaces fields like ./., ./.:.:.:... with a consistent FORMAT-compatible string.
    - Also fixes rare cases like '0' → '0/0', '1' → '1/1'
    - Saves the result as a compressed .vcf.gz file
    """
    vcf_out = vcf_path.replace('.vcf.gz', '_preprocessed.vcf.gz').replace('.vcf', '_preprocessed.vcf.gz')
    opener = gzip.open if vcf_path.endswith('.gz') else open
    read_mode = 'rt' if vcf_path.endswith('.gz') else 'r'

    with opener(vcf_path, read_mode) as infile, gzip.open(vcf_out, 'wt') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                cols = line.strip().split('\t')
                format_fields = cols[8].split(':')

                # Dynamically build the replacement string
                replacement = '0/0'
                for field in format_fields[1:]:
                    if field in ['GQ', 'PL']:
                        replacement += ':30'
                    else:
                        replacement += ':1'

                fixed = cols[:9]
                for sample in cols[9:]:
                    # Missing data: ./., ./.:.:.:.:...
                    if replace_missing and (sample == './.' or sample.startswith('./.:')):
                        fixed.append(replacement)
                    # Rare cases: just '0' or '1'
                    elif sample == '0':
                        fixed.append('0/0')
                    elif sample == '1':
                        fixed.append('1/1')
                    else:
                        fixed.append(sample)

                outfile.write('\t'.join(fixed) + '\n')

    print(f"✅ Preprocessed file saved as: {vcf_out}")
    return vcf_out
