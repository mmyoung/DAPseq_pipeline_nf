process MARK_DUPLICATES {
    tag "Mark duplicate"
    label 'Mark_duplication'
    publishDir "${params.output_dir}/deduplicateion_out", mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}_dedup_Q20_sorted.bam"), emit: bam
    tuple val(meta), path("${meta.id}_dedup_Q20_sorted.bam.bai"), emit: bai
    script:

    """
    module load samtools/1.13

    java -jar /project/gzy8899/softwares/picard.jar \\
        MarkDuplicates \\
        --INPUT $bam_file \\
        --OUTPUT ${meta.id}_dedup.bam \\
        --METRICS_FILE ${meta.id}.MarkDuplicates.metrics.txt \\
        --REMOVE_DUPLICATES
        
    samtools view -bq 20  ${meta.id}_dedup.bam | samtools sort - > ${meta.id}_dedup_Q20_sorted.bam

    samtools index ${meta.id}_dedup_Q20_sorted.bam
    
    """
     

}