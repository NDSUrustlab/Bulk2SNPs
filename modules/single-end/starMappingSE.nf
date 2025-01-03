process STARmapping_SE {

    input:
        path star_index  
        tuple val(sample_id), path(trimmed_reads)

    output:
        path "PASS2_output/*/*.bam", emit: bam_file

    script:
    """
    mkdir -p PASS1_output/\$(basename ${trimmed_reads} _trim.fq.gz) PASS2_output/\$(basename ${trimmed_reads} _trim.fq.gz)

    STAR \
    --runThreadN 16 \
    --genomeDir $star_index \
    --readFilesIn ${trimmed_reads} \
    --outFileNamePrefix PASS1_output/\$(basename ${trimmed_reads} _trim.fq.gz)/\$(basename ${trimmed_reads} _trim.fq.gz)_  \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 20 \
    --outSAMunmapped Within \
    --readFilesCommand zcat

    STAR \
    --runThreadN 16 \
    --genomeDir $star_index \
    --readFilesIn ${trimmed_reads} \
    --outFileNamePrefix PASS2_output/\$(basename ${trimmed_reads} _trim.fq.gz)/\$(basename ${trimmed_reads} _trim.fq.gz)_  \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 20 \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --sjdbFileChrStartEnd PASS1_output/\$(basename ${trimmed_reads} _trim.fq.gz)/\$(basename ${trimmed_reads} _trim.fq.gz)_SJ.out.tab
    """

}