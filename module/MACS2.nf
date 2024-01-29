process MACS2_CALLPEAK {
    tag "$sample_id"
    label 'Peak calling with MACS3'
    publishDir "${params.output_dir}/macs3_output", mode: 'copy'

    conda  "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"


    input:
    tuple val(sample_id), path(ipbam), path(controlbam)
    val   macs2_gsize

    output:
    tuple val(sample_id), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(sample_id), path("*.xls")                   , emit: xls

    tuple val(sample_id), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(sample_id), path("*.bed")       , optional:true, emit: bed
    tuple val(sample_id), path("*.bdg")       , optional:true, emit: bdg

    script:
    def control   = controlbam ? "--control $controlbam" : ''

    """
    macs2 \\
        callpeak \\
        --gsize $macs2_gsize \\
        --format BAMPE \\
        --name $sample_id \\
        --treatment $ipbam \\
        $control
    """
}