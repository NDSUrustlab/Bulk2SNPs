process BWA_Index {

    input:
        path genome
    
    output:
        path "*", emit: bwa_index
    
    publishDir "${params.outdir}/index", mode: 'copy'

    script:
    """
    bwa-mem2 index ${genome}
    """
}