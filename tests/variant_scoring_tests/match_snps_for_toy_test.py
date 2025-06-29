def match_vcf_with_predictions(vcf_file, prediction_file, output_file):
 # Extract SNP IDs from the VCF's ID column (column 3, index 2)
    vcf_ids = set()

    vcf = open(vcf_file, 'r')
    for line in vcf:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 2:
            vcf_ids.add(parts[2])
    vcf.close()

    # Filter predictions based on SNP column
    pred = open(prediction_file, 'r')
    out = open(output_file, 'w')

    header = next(pred)
    out.write(header)

    for line in pred:
        parts = line.strip().split()
        if not parts:
            continue
        snp = parts[0]
        if snp in vcf_ids:
            out.write(line + '\n')

    pred.close()
    out.close()


# Example usage
match_vcf_with_predictions('toy_trimmed2.vcf', '/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/sed_all_chr/chr18.txt', 'matched_snps_output2.vcf')
