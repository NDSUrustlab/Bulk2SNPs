process STAR_Index {

    input:
        path genome
        path gff

    output:
        path "Log.out", emit: index_file
        path "*"
    
    publishDir "${params.outdir}/index", mode: 'copy'

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir . \
    --genomeFastaFiles ${genome} \
    --sjdbGTFfile ${gff} \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFtagExonParentGene ID \
    --sjdbOverhang 100 \
    --limitGenomeGenerateRAM 103020153610
    """

}