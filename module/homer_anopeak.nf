process HOMER_ANNOTATEPEAKS {
    tag "$sample_id"
    label 'Peak_annotation'
    publishDir "${params.output_dir}/macs3_output", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    tuple val(sample_id), path(peak)
    path  fasta
    path  gtf

    output:
    tuple val(sample_id), path("*annotatePeaks.txt"), emit: txt

    script:
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${sample_id}.annotatePeaks.txt
    """
}