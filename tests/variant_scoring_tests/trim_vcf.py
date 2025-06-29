def extract_vcf_sample(input_file, output_file, max_data_rows=40, max_columns=13):
    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')

    data_row_count = 0
    found_header = False

    for line in infile:
        if not found_header:
            if line.startswith('#CHROM'):
                header_cols = line.strip().split('\t')[:13]
                outfile.write('\t'.join(header_cols) + '\n')
                found_header = True
            continue  
        if data_row_count < 10:
            cols = line.strip().split('\t')[:13]
            outfile.write('\t'.join(cols) + '\n')
            data_row_count += 1
        else:
            break

    infile.close()
    outfile.close()
# Example usage
input_vcf = '/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/enroll_hd/regulatory_vcfs/gwa12345.mis4.9064.hg38.cisreg0.5.mad0.01.chr18.vcf'
output_vcf = 'toy_trimmed2.vcf'
extract_vcf_sample(input_vcf, output_vcf)
