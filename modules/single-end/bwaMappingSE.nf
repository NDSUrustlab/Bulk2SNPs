process BWAmapping_SE {

    input:
        path bwa_index
        val genome
        tuple val(sample_id), path(trimmed_reads)

    output:
        path "*.sorted.bam", emit: bam_file
        path "*.sorted.bam.csi"
    
    publishDir "${params.outdir}/mapped", mode: 'copy'

    script:
    """
    bwa-mem2 mem -t 16 ${genome} ${trimmed_reads} > \$(basename ${trimmed_reads} _trim.fq.gz).sam

    samtools view -Sb \$(basename ${trimmed_reads} _trim.fq.gz).sam > \$(basename ${trimmed_reads} _trim.fq.gz).bam
    samtools sort -@ 16 \$(basename ${trimmed_reads} _trim.fq.gz).bam -o \$(basename ${trimmed_reads} _trim.fq.gz).sorted.bam

    samtools index -c \$(basename ${trimmed_reads} _trim.fq.gz).sorted.bam
    """
        
}