#!/usr/bin/env nextflow

// Parameters
params.datatype = "DNA" // DNA or RNA
params.seqtype = "PE" //PE for paired-end or SE for single-end sequences
params.reads = "${projectDir}/data"
params.bulk1 = "bulk1"
params.bulk2 = "bulk2"
params.genome = "${projectDir}/genome/genome.fasta"
params.gff = "${projectDir}/genome/annotation.gff"
params.outdir = "${projectDir}/results"
params.index_dir = "${projectDir}/index"
params.help = false

//////////////////////////////////////////////////
// Reading the raw files with different extensions
//////////////////////////////////////////////////

if (params.seqtype == "PE") {

    if (file("${params.reads}/${params.bulk1}").isDirectory()) {


        bulk1_ch = Channel.fromFilePairs("${params.reads}/${params.bulk1}/${params.bulk1}_{1,2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)

        bulk2_ch = Channel.fromFilePairs("${params.reads}/${params.bulk2}/${params.bulk2}_{1,2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)

    } else {

        bulk1_ch = Channel.fromFilePairs("${params.reads}/${params.bulk1}_{1,2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)

        bulk2_ch = Channel.fromFilePairs("${params.reads}/${params.bulk2}_{1,2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)
    }

} else {

    if (file("${params.reads}/${params.bulk1}").isDirectory()) {


        bulk1_ch = Channel.fromPath("${params.reads}/${params.bulk1}/${params.bulk1}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)
        .map { file_path ->
            def sample_id = file_path.baseName
            tuple(sample_id, file_path) // Create a tuple of sample_id and file_path
        }

        bulk2_ch = Channel.fromPath("${params.reads}/${params.bulk2}/${params.bulk2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)
        .map { file_path ->
            def sample_id = file_path.baseName
            tuple(sample_id, file_path) // Create a tuple of sample_id and file_path
        }

    } else {

        bulk1_ch = Channel
        .fromPath("${params.reads}/${params.bulk1}.{fq,fq.gz,fastq,fastq.gz}", checkIfExists: true)
        .map { file_path ->
            def sample_id = file_path.baseName
            tuple(sample_id, file_path) // Create a tuple of sample_id and file_path
        }

        bulk2_ch = Channel.fromPath("${params.reads}/${params.bulk2}.{fq,fq.gz,fastq,fastq.gz}" , checkIfExists: true)
        .map { file_path ->
            def sample_id = file_path.baseName
            tuple(sample_id, file_path) // Create a tuple of sample_id and file_path
        }
    }
}

if( params.help ) {

log.info """
    ==============================================
    ========= B u l k 2 S N P s Pipeline =========
    ==============================================
    Usage:
    nextflow run main.nf --bulk1 <Bulk1_name> --bulk2 <Bulk2_name> 
                         --genome <genome_dir> --gff <annotation.gff3>
    
    Parameters:
      --bulk1        Name of bulk1 [required]
      --bulk2        Name of bulk2 [required]
      --reads        Path to the raw reads [default: data/]
      --genome       Path to the reference genome directory [required]
      --gff          Path to the annotations directory [required]
      --index        Path to the genome index directory
      --datatype     Type of data DNA or RNA [default: DNA]
      --seqtype      Type of sequencing used, paried-end (PE) or single end (SE) [default: PE]
      --outdir       Name of the output directory [default: results/]
      --help         Display this help message
    """
exit 0
}


log.info """\
    B u l k 2 S N P s - P I P E L I N E
    ===================================
    datatype: ${params.datatype}
    sequencing: ${params.seqtype}
    Bulk 1: ${params.bulk1}
    Bulk 2: ${params.bulk2}
    Genome: ${params.genome}
    Annotation file: ${params.gff}
    Output: ${params.outdir}
    """
    .stripIndent(true)

include { FASTQC as QC_BULK1 } from './modules/fastqc.nf'
include { FASTQC as QC_BULK2 } from './modules/fastqc.nf'

include { FASTP_SE as TRIMMING_SE_BULK1 } from './modules/single-end/fastpSE.nf'
include { FASTP_SE as TRIMMING_SE_BULK2 } from './modules/single-end/fastpSE.nf'

include { FASTP as TRIMMING_BULK1 } from './modules/fastp.nf'
include { FASTP as TRIMMING_BULK2 } from './modules/fastp.nf'

include { FASTQC_trim as QC_TRIM_Bulk1 } from './modules/fastqcTrim.nf'
include { FASTQC_trim as QC_TRIM_Bulk2 } from './modules/fastqcTrim.nf'

include { MULTIQC } from './modules/multiQC.nf'

include { STAR_Index as STAR_INDEXDING } from './modules/starIndex.nf'

include { BWA_Index as BWA_INDEXDING } from './modules/bwaIndex.nf'

include { STARmapping_SE as STAR_MAPPING_SE_Bulk1 } from './modules/single-end/starMappingSE.nf'
include { STARmapping_SE as STAR_MAPPING_SE_Bulk2 } from './modules/single-end/starMappingSE.nf'

include { BWAmapping_SE as BWA_MAPPING_SE_Bulk1 } from './modules/single-end/bwaMappingSE.nf'
include { BWAmapping_SE as BWA_MAPPING_SE_Bulk2 } from './modules/single-end/bwaMappingSE.nf'

include { STARmapping as STAR_MAPPING_Bulk1 } from './modules/starMapping.nf'
include { STARmapping as STAR_MAPPING_Bulk2 } from './modules/starMapping.nf'

include { BWAmapping as BWA_MAPPING_Bulk1 } from './modules/bwaMapping.nf'
include { BWAmapping as BWA_MAPPING_Bulk2 } from './modules/bwaMapping.nf'

include { BAMoutput as BAM_OUTPUT } from './modules/cleanBAMsaving.nf'

include { gatkHaplotypeCallingRNA as HAPLOTYPE_CALLING_RNA_BULK1 } from './modules/haplotypeCallingRNA.nf'
include { gatkHaplotypeCallingRNA as HAPLOTYPE_CALLING_RNA_BULK2 } from './modules/haplotypeCallingRNA.nf'

include { gatkHaplotypeCallingDNA as HAPLOTYPE_CALLING_DNA_BULK1 } from './modules/haplotypeCallingDNA.nf'
include { gatkHaplotypeCallingDNA as HAPLOTYPE_CALLING_DNA_BULK2 } from './modules/haplotypeCallingDNA.nf'

include { snpCalling as SNP_CALLING } from './modules/snpCalling.nf'



workflow RNA_SE {

    // Adapter trimming
    bulk1trim_ch = TRIMMING_SE_BULK1(bulk1_ch)
    bulk2trim_ch = TRIMMING_SE_BULK2(bulk2_ch)

    // Quality check (before and after trimming)
    bulk1_QC_ch = QC_BULK1(bulk1_ch)
    bulk2_QC_ch = QC_BULK2(bulk2_ch)
    QC_ch = bulk1_QC_ch.mix(bulk2_QC_ch).collect()

    bulk1_trimQC_ch = QC_TRIM_Bulk1(bulk1trim_ch.trimmed_reads)
    bulk2_trimQC_ch = QC_TRIM_Bulk2(bulk2trim_ch.trimmed_reads)
    trimQC_ch = bulk1_trimQC_ch.mix(bulk2_trimQC_ch).collect()

    // MultiQC reports
    MULTIQC(QC_ch.mix(trimQC_ch).collect())

    genome_ch = Channel.fromPath(params.genome)
    gff_ch = Channel.fromPath(params.gff)

    // ${params.index_dir}.view()
    
    // STAR indexing
    if ( file(params.index_dir).exists() ) {
        // Use the provided index directory
        // println("running from dir")
        star_index_ch = Channel.fromPath(params.index_dir)
        // // STAR mapping
        star_mapping_bulk1_ch = STAR_MAPPING_SE_Bulk1(star_index_ch, bulk1trim_ch.trimmed_reads)
        star_mapping_bulk2_ch = STAR_MAPPING_SE_Bulk2(star_index_ch, bulk2trim_ch.trimmed_reads)
    } else {
        // Generate the index directory
        // println("Generating")
        star_index_ch = STAR_INDEXDING(genome_ch, gff_ch)
        star_mapping_bulk1_ch = STAR_MAPPING_SE_Bulk1(star_index_ch.index_file.map{ file -> file.parent }, bulk1trim_ch.trimmed_reads)
        star_mapping_bulk2_ch = STAR_MAPPING_SE_Bulk2(star_index_ch.index_file.map{ file -> file.parent }, bulk2trim_ch.trimmed_reads)
    }

    // Saving Bam files and their index
    BAM_OUTPUT(star_mapping_bulk1_ch.bam_file.mix(star_mapping_bulk2_ch.bam_file).flatten())

    // Haplotype calling
    bulk1_GHC_ch = HAPLOTYPE_CALLING_RNA_BULK1(star_mapping_bulk1_ch.bam_file, genome_ch)
    bulk2_GHC_ch = HAPLOTYPE_CALLING_RNA_BULK2(star_mapping_bulk2_ch.bam_file, genome_ch)

    // SNP filtering
    SNP_CALLING(genome_ch, bulk1_GHC_ch.gvcf, bulk2_GHC_ch.gvcf)

}

workflow RNA_PE {

    // Adapter trimming
    bulk1trim_ch = TRIMMING_BULK1(bulk1_ch)
    bulk2trim_ch = TRIMMING_BULK2(bulk2_ch)

    // Quality check (before trimming)
    bulk1_QC_ch = QC_BULK1(bulk1_ch)
    bulk2_QC_ch = QC_BULK2(bulk2_ch)
    QC_ch = bulk1_QC_ch.mix(bulk2_QC_ch).collect()

    // Quality check (after trimming)
    bulk1_trimQC_ch = QC_TRIM_Bulk1(bulk1trim_ch)
    bulk2_trimQC_ch = QC_TRIM_Bulk2(bulk2trim_ch)
    trimQC_ch = bulk1_trimQC_ch.mix(bulk2_trimQC_ch).collect()
   
    // MultiQC reports
    MULTIQC(QC_ch.mix(trimQC_ch).collect())

    genome_ch = Channel.fromPath(params.genome)
    gff_ch = Channel.fromPath(params.gff)

    // STAR indexing
    if ( file(params.index_dir).exists() ) {
        // Use the provided index directory
        // println("running from dir")
        star_index_ch = Channel.fromPath(params.index_dir)
        // // STAR mapping
        star_mapping_bulk1_ch = STAR_MAPPING_Bulk1(star_index_ch, bulk1trim_ch.trimmed_reads)
        star_mapping_bulk2_ch = STAR_MAPPING_Bulk2(star_index_ch, bulk2trim_ch.trimmed_reads)
    } else {
        // Generate the index directory
        // println("Generating")
        star_index_ch = STAR_INDEXDING(genome_ch, gff_ch)
        star_mapping_bulk1_ch = STAR_MAPPING_Bulk1(star_index_ch.index_file.map{ file -> file.parent }, bulk1trim_ch.trimmed_reads)
        star_mapping_bulk2_ch = STAR_MAPPING_Bulk2(star_index_ch.index_file.map{ file -> file.parent }, bulk2trim_ch.trimmed_reads)
    }

    // Saving Bam files and their index
    BAM_OUTPUT(star_mapping_bulk1_ch.bam_file.mix(star_mapping_bulk2_ch.bam_file).flatten())

    // Haplotype calling
    bulk1_GHC_ch = HAPLOTYPE_CALLING_RNA_BULK1(star_mapping_bulk1_ch.bam_file, genome_ch)
    bulk2_GHC_ch = HAPLOTYPE_CALLING_RNA_BULK2(star_mapping_bulk2_ch.bam_file, genome_ch)

    // SNP filtering
    SNP_CALLING(genome_ch, bulk1_GHC_ch.gvcf, bulk2_GHC_ch.gvcf)

}

workflow DNA_SE {

    // Adapter trimming
    bulk1trim_ch = TRIMMING_SE_BULK1(bulk1_ch)
    bulk2trim_ch = TRIMMING_SE_BULK2(bulk2_ch)

    // Quality check (before and after trimming)
    bulk1_QC_ch = QC_BULK1(bulk1_ch)
    bulk2_QC_ch = QC_BULK2(bulk2_ch)
    QC_ch = bulk1_QC_ch.mix(bulk2_QC_ch).collect()

    bulk1_trimQC_ch = QC_TRIM_Bulk1(bulk1trim_ch.trimmed_reads)
    bulk2_trimQC_ch = QC_TRIM_Bulk2(bulk2trim_ch.trimmed_reads)
    trimQC_ch = bulk1_trimQC_ch.mix(bulk2_trimQC_ch).collect()

    // MultiQC reports
    MULTIQC(QC_ch.mix(trimQC_ch).collect())

    genome_ch = Channel.fromPath(params.genome)

    // BWA indexing
    if ( file(params.index_dir).exists() ) {
        // Use the provided index directory
        // println("running from dir")
        bwa_index_ch = Channel.fromPath(params.index_dir)

    } else {
        // Generate the index directory
        // println("Generating")
        bwa_index_ch = BWA_INDEXDING(genome_ch)

    }

    // BWA mapping
    bwa_mapping_bulk1_ch = BWA_MAPPING_SE_Bulk1(bwa_index_ch, genome_ch.map{file -> file.Name}, bulk1trim_ch.trimmed_reads)
    bwa_mapping_bulk2_ch = BWA_MAPPING_SE_Bulk2(bwa_index_ch, genome_ch.map{file -> file.Name}, bulk2trim_ch.trimmed_reads)
    
    // Haplotype calling
    bwa_bulk1_GHC_ch = HAPLOTYPE_CALLING_DNA_BULK1(bwa_mapping_bulk1_ch.bam_file, genome_ch)
    bwa_bulk2_GHC_ch = HAPLOTYPE_CALLING_DNA_BULK2(bwa_mapping_bulk2_ch.bam_file, genome_ch)

    // SNP filtering
    SNP_CALLING(genome_ch, bwa_bulk1_GHC_ch.gvcf, bwa_bulk2_GHC_ch.gvcf)

}

workflow DNA_PE {

    // Adapter trimming
    bulk1trim_ch = TRIMMING_BULK1(bulk1_ch)
    bulk2trim_ch = TRIMMING_BULK2(bulk2_ch)

    // Quality check (before trimming)
    bulk1_QC_ch = QC_BULK1(bulk1_ch)
    bulk2_QC_ch = QC_BULK2(bulk2_ch)
    QC_ch = bulk1_QC_ch.mix(bulk2_QC_ch).collect()

    // Quality check (after trimming)
    bulk1_trimQC_ch = QC_TRIM_Bulk1(bulk1trim_ch)
    bulk2_trimQC_ch = QC_TRIM_Bulk2(bulk2trim_ch)
    trimQC_ch = bulk1_trimQC_ch.mix(bulk2_trimQC_ch).collect()
   
    // MultiQC reports
    MULTIQC(QC_ch.mix(trimQC_ch).collect())

    genome_ch = Channel.fromPath(params.genome)

    // BWA indexing
    if ( file(params.index_dir).exists() ) {
        // Use the provided index directory
        // println("running from dir")
        bwa_index_ch = Channel.fromPath(params.index_dir)

    } else {
        // Generate the index directory
        // println("Generating")
        bwa_index_ch = BWA_INDEXDING(genome_ch)

    }

    // BWA mapping
    bwa_mapping_bulk1_ch = BWA_MAPPING_Bulk1(bwa_index_ch, genome_ch.map{file -> file.Name}, bulk1trim_ch.trimmed_reads)
    bwa_mapping_bulk2_ch = BWA_MAPPING_Bulk2(bwa_index_ch, genome_ch.map{file -> file.Name}, bulk2trim_ch.trimmed_reads)
    
    // Haplotype calling
    bwa_bulk1_GHC_ch = HAPLOTYPE_CALLING_DNA_BULK1(bwa_mapping_bulk1_ch.bam_file, genome_ch)
    bwa_bulk2_GHC_ch = HAPLOTYPE_CALLING_DNA_BULK2(bwa_mapping_bulk2_ch.bam_file, genome_ch)

    // SNP filtering
    SNP_CALLING(genome_ch, bwa_bulk1_GHC_ch.gvcf, bwa_bulk2_GHC_ch.gvcf)

}

workflow {
    if (params.datatype == "RNA" && params.seqtype == "SE") {
        RNA_SE()
    } 
    else if (params.datatype == "RNA" && params.seqtype == "PE") {
        RNA_PE()
    } 
    else if (params.datatype == "DNA" && params.seqtype == "SE") {
        DNA_SE()
    } 
    else if (params.datatype == "DNA" && params.seqtype == "PE") {
        DNA_PE()
    } 
    else {
        error "Unsupported datatype (${params.datatype}) or seqtype (${params.seqtype})."
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Please find your SNP table file for downstream QTL anlysis with QTLseqR in --> $params.outdir/final_SNPs\n" : "Oops .. something went wrong" )
}


