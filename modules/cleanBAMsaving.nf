process BAMoutput {

    input:
        path bam_file

    output:
        path "*.bam"
        path "*.bam.csi"
    
    publishDir "${params.outdir}/mapped", mode: 'copy'

    script:
    """
    # Create BAM index files
    mv ${bam_file} \$(basename ${bam_file} _Aligned.sortedByCoord.out.bam).sorted.bam
    samtools index -c \$(basename ${bam_file} _Aligned.sortedByCoord.out.bam).sorted.bam
    """
}