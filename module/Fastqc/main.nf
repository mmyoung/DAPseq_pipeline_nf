process FASTQC {
    tag "FASTQC"

    input:
        tuple val(meta), path(reads)

    output:
        path("${meta.id}/*")

    publishDir "${params.output_dir}/clean_reads_fastQC", mode: 'copy'

    script:
    def prefix = "${meta.id}"
    if (meta.single_end) {
        """
        echo "FastQC for single-end data..."
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        mkdir ${prefix}
        /project/zhuzhuzhang/lyang/software/FastQC/fastqc -o ${prefix} --threads $task.cpus ${prefix}.fastq.gz
        """
    } 
    else {
        """
        echo "FastQC for pair-end data..."
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        mkdir ${prefix}
        /project/zhuzhuzhang/lyang/software/FastQC/fastqc -o ${prefix} --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
        """
    }
}