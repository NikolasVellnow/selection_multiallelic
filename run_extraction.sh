echo "Running freebayes\n"
freebayes -f /vol/storage/dest_linvilla_bams/dmel-all-chromosome-r6.12.fasta --pooled-continuous -r 3R -L bam_list_RGfixed.txt > variants_3R_pooled.vcf
echo "Compressing VCF-file\n"
bgzip variants_3R_pooled.vcf
echo "Creating index for compressed VCF-file\n"
tabix variants_3R_pooled.vcf.gz
echo "Extracting allele counts with bcftools\n"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' variants_3R_pooled.vcf.gz > allele_counts_snps_3R.tsv

echo "Filtering with awk"
awk 'BEGIN { FS="\t"; OFS="\t" }
{
    # Split ALT field (column 4) on commas.
    n = split($4, altAlleles, ",");
    # Only consider sites with at least 2 alternate alleles (i.e. REF + >=2 ALTs = 3 alleles)
    if(n >= 2) {
        valid = 1;
        # Check that all sample fields (columns 5 to NF) are not missing (i.e. not ".")
        for(i = 5; i <= NF; i++){
            if($i == ".") { valid = 0; break; }
        }
        if(valid) print $0;
    }
}' allele_counts_snps_3R.tsv > filtered_sites.tsv

