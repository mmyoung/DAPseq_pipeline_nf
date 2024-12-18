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

    java -Xmx20g -jar /project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env/share/picard-2.18.23-0/picard.jar \\
        MarkDuplicates \\
        I=$bam_file \\
        O=${meta.id}_dedup.bam \\
        M=${meta.id}.MarkDuplicates.metrics.txt \\
        REMOVE_DUPLICATES=true
        
    samtools view -bq 20  ${meta.id}_dedup.bam | samtools sort - > ${meta.id}_dedup_Q20_sorted.bam

    samtools index ${meta.id}_dedup_Q20_sorted.bam
    
    """
     

}
