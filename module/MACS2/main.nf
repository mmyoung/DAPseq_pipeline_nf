nextflow.enable.dsl=2

process MACS2_CALLPEAK {
    tag "MACS2 calling peak for ${meta.id}"
    label 'Peak_calling_with_MACS3'
    publishDir "${params.output_dir}/macs3_output", mode: 'copy'

    conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    tuple val(meta), path(ipbam), val(controlbam)
    val   macs2_gsize

    output:
    tuple val(meta.id), path("*.narrowPeak"), emit: peak
    tuple val(meta.id), path("*.xls")       , emit: xls
    tuple val(meta.id), path("*.gappedPeak"), optional: true, emit: gapped
    tuple val(meta.id), path("*.bed")       , optional: true, emit: bed
    tuple val(meta.id), path("*.bdg")       , optional: true, emit: bdg

    script:
    // Handle case where controlbam is null
    def control = controlbam ? "--control $controlbam" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = meta.single_end ? 'BAM' : 'BAMPE'

    """
    macs3 \\
        callpeak \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --call-summits \\
        --keep-dup 1 \\
        --treatment $ipbam \\
        $control
    """
}