process gatkHaplotypeCallingDNA {

    input:
        path bam_file
        path genome
        
    output:
        path "*.g.vcf", emit: gvcf

    publishDir "${params.outdir}/variant_calling", mode: 'copy'

    script:
    """
    # indexing genome
    samtools faidx ${genome}

    # genome dict
    gatk CreateSequenceDictionary -R ${genome} -O \$(basename $genome .fasta).dict

    # Add read groups
    gatk AddOrReplaceReadGroups -I ${bam_file} -O \$(basename $bam_file .sorted.bam)_with_rg.bam --RGID LH00000.000 --RGLB Lib1 --RGPL ILLUMINA --RGPU FC00HHH0LT0.0 --RGSM \$(basename $bam_file .sorted.bam)

    # GATK MarkDuplicates
    gatk MarkDuplicates -I \$(basename $bam_file .sorted.bam)_with_rg.bam -O \$(basename $bam_file .sorted.bam).mrkdups.bam -M \$(basename $bam_file .sorted.bam).mrkdups.metrics.txt

    # Samtools CSI Index
    samtools index -c \$(basename $bam_file .sorted.bam).mrkdups.bam

    # GATK HaplotypeCaller
    gatk HaplotypeCaller -R ${genome} -I \$(basename $bam_file .sorted.bam).mrkdups.bam -O \$(basename $bam_file .sorted.bam).g.vcf -ERC GVCF
    """

}