process gatkHaplotypeCallingRNA {

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
    gatk AddOrReplaceReadGroups -I ${bam_file} -O \$(basename $bam_file _Aligned.sortedByCoord.out.bam)_with_rg.bam --RGID LH00000.000 --RGLB Lib1 --RGPL ILLUMINA --RGPU FC00HHH0LT0.0 --RGSM \$(basename $bam_file _Aligned.sortedByCoord.out.bam)

    # GATK MarkDuplicates
    gatk MarkDuplicates -I \$(basename $bam_file _Aligned.sortedByCoord.out.bam)_with_rg.bam -O \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.bam -M \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.metrics.txt

    # GATK SplitNCigarReads
    gatk SplitNCigarReads -R ${genome} -I \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.bam -O \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.spliltNcigarReads.bam --create-output-bam-index false

    # Samtools CSI Index
    samtools index -c \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.spliltNcigarReads.bam

    # GATK HaplotypeCaller
    gatk HaplotypeCaller -R ${genome} -I \$(basename $bam_file _Aligned.sortedByCoord.out.bam).mrkdups.spliltNcigarReads.bam -O \$(basename $bam_file _Aligned.sortedByCoord.out.bam).g.vcf -ERC GVCF
    """

}