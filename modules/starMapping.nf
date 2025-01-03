process STARmapping {

    input:
        path star_index 
        tuple val(sample_id), path(trimmed_reads)

    output:
        path "PASS2_output/*/*.bam", emit: bam_file

    script:
    """
    mkdir -p PASS1_output/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz) PASS2_output/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)

    STAR \
    --runThreadN 16 \
    --genomeDir $star_index \
    --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \
    --outFileNamePrefix PASS1_output/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)_  \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 20 \
    --outSAMunmapped Within \
    --readFilesCommand zcat

    STAR \
    --runThreadN 16 \
    --genomeDir $star_index \
    --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \
    --outFileNamePrefix PASS2_output/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)_  \
    --outReadsUnmapped Fastx \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 20 \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --sjdbFileChrStartEnd PASS1_output/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)/\$(basename ${trimmed_reads[0]} _1_trim.fq.gz)_SJ.out.tab
    """

}