process BOWTIE2MAP {
    tag "Mapping clean reads to bowtie index"
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*map_sorted.bam"), emit: bam
        tuple val(meta), path("*.log"), emit: log

    script:

    def prefix = "${meta.id}"
    def reads_args = ""
    if (meta.single_end) {
        reads_args = "-U ${reads}"
    } else {
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
    }

    """
    module load samtools/1.13
    /project/zhuzhuzhang/lyang/software/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 \\
    -x ${params.bowtie_idx} \\
    $reads_args \\
    --threads $task.cpus \\
    2> ${prefix}.bowtie2.log \\
    | samtools view -bS | samtools sort - > ${prefix}_map_sorted.bam
    
    """
}