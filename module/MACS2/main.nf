nextflow.enable.dsl=2

process MACS2_CALLPEAK {
    tag "MACS2 calling peak for ${meta.id}"
    label 'Peak_calling_with_MACS3'
    publishDir "${params.output_dir}/macs3_output", mode: 'copy'

    // conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    
    tuple val(meta), path(ipbam, stageAs: 'treatment.bam'), val(control_path)
    val   macs2_gsize

    output:
    tuple val(meta.id), path("*.narrowPeak"), emit: peak
    tuple val(meta.id), path("*.xls")       , emit: xls
    tuple val(meta.id), path("*.gappedPeak"), optional: true, emit: gapped
    tuple val(meta.id), path("*.bed")       , optional: true, emit: bed
    tuple val(meta.id), path("*.bdg")       , optional: true, emit: bdg

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = meta.single_end ? 'BAM' : 'BAMPE'

    """
    echo "=== MACS2 Debug ==="
    echo "Sample: ${meta.id}"
    echo "Treatment: $ipbam" 
    echo "Control: $control_path"

    if [ "$control_path" != "null" ] && [ -f "$control_path" ]; then
        echo "Using control: $control_path"
        CONTROL_CMD="--control $control_path"
    else
        echo "No control used"
        CONTROL_CMD=""
    fi

    macs3 \\
        callpeak \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --call-summits \\
        --keep-dup 1 \\
        --treatment $ipbam \\
        \$CONTROL_CMD
    """
}    