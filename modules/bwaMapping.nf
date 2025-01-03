process BWAmapping {

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
    bwa-mem2 mem -t 16 ${genome} ${trimmed_reads[0]} ${trimmed_reads[1]} > \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).sam

    samtools view -Sb \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).sam > \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).bam
    samtools sort -@ 16 \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).bam -o \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).sorted.bam
    
    samtools index -c \$(basename ${trimmed_reads[0]} _1_trim.fq.gz).sorted.bam
    """
        
}