process snpCalling {
    input:
        path genome
        path bulk1
        path bulk2
        
    output:
        path "Bulks.vcf", emit: joint_genotyped_vcf
        path "*SNPs.vcf", emit: snp_vcf
        path "*SNPs.table", emit: snp_table

    publishDir "${params.outdir}/variant_calling/", mode: 'copy', pattern: "Bulks.vcf"
    publishDir "${params.outdir}/final_SNPs/", mode: 'copy', pattern: "*SNPs.vcf"
    publishDir "${params.outdir}/final_SNPs/", mode: 'copy', pattern: "*SNPs.table"

    script:
    """
    # indexing genome
    samtools faidx ${genome}

    # genome dict
    gatk CreateSequenceDictionary -R ${genome} -O \$(basename $genome .fasta).dict

    #combining gvcfs
    gatk CombineGVCFs -R ${genome} -V ${bulk1} -V ${bulk2} -O Bulks_combined.g.vcf

    #joint genotyping
    gatk GenotypeGVCFs -R ${genome} -V Bulks_combined.g.vcf -O Bulks.vcf

    #select variants (SNPs)
    gatk SelectVariants -R ${genome} -V Bulks.vcf --select-type-to-include SNP -O Bulks_SNPs.vcf

    #snp table
    gatk VariantsToTable -R ${genome} -V Bulks_SNPs.vcf -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -O Bulks_SNPs.table
    """
}