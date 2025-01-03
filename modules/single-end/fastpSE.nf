process FASTP_SE {

    // tag "Adapter trimming of $sample_id"

    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("*_trim.fq.gz"), emit: trimmed_reads

    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    script:
    """
    fastp -i ${reads[0]} -o ${sample_id}_trim.fq.gz --thread 16 --average_qual 15 --length_required 50
    """
    
}