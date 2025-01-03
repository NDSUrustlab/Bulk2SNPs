process FASTQC_trim {

    // tag "Quality Check on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*"

    script:
    """
    fastqc -t 2 -q ${reads}
    """

}