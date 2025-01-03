process FASTP {
    
    // tag "Adapter trimming of $sample_id"

    input:
        tuple val(sample_id), path(reads)
    
    output:
        tuple val(sample_id), path("*_trim.fq.gz"), emit: trimmed_reads

    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_1_trim.fq.gz -O ${sample_id}_2_trim.fq.gz --thread 16 --detect_adapter_for_pe --average_qual 15 --length_required 50 
    """
    
}
