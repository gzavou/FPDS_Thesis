def normalize_snp_id(chrom, pos, ref, alt):
    """
    Standardize SNP ID format as chrom_pos_REF_ALT (e.g., '15_31254588_G_A').
    Ensures REF < ALT lexicographically for consistent comparison.
    """
    chrom = str(chrom).replace("chr", "")  # remove 'chr' if present
    if ref > alt:
        ref, alt = alt, ref  # standardize order to minimize REF/ALT issues
    return f"{chrom}_{pos}_{ref}_{alt}"


# Extract normalized SNPs from VCF
vcf_snps = set()
with open("/pool01/projects/abante_lab/ao_prediction_enrollhd_2024/enroll_hd/vcfs/gwa12345.mis4.9064.hg38.chr15.vcf", "r") as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
        snp_id = normalize_snp_id(chrom, pos, ref, alt)
        vcf_snps.add(snp_id)

# Extract normalized SNPs from prediction file
pred_snps = set()
with open("/pool01/projects/abante_lab/genomic_llms/borzoi/test_outputs/test_02.txt", "r") as pred:
    for line in pred:
        if line.startswith("snp"):  # skip header
            continue
        snp_id = line.split("\t")[0]
        # normalize SNP string from prediction if needed
        parts = snp_id.strip().split("_")
        if len(parts) == 4:
            snp_id = normalize_snp_id(*parts)
        pred_snps.add(snp_id)

# Find and print matching SNPs
common_snps = vcf_snps.intersection(pred_snps)
print(f"Matching SNPs: {len(common_snps)}")
with open("matching_snps.txt", "w") as f:
    for snp in sorted(common_snps):
        f.write(snp + "\n")
