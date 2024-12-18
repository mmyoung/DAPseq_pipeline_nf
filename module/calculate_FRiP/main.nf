process CAL_FRIP {
    conda '/project/gzy8899/lyang/DAPseq_pipeline_nf/env_requirment.yml'

    tag "FRiP calculation"
    label 'FRiP'
    publishDir "${params.output_dir}/FRiP_score", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    input:
    tuple val(meta), path(peak)
    tuple val(meta), path(bam)

    output:
    path("*.txt"), emit: txt

    shell:

    """
    module load samtools

    # Count reads in the peaks
    peak_read=\$(samtools view -c -L ${peak} ${bam})

    # Count total reads in the BAM file
    total_read=\$(samtools view -c ${bam})

    # Calculate the FRiP score
    score=\$(awk "BEGIN {print \$peak_read / \$total_read}")

    # Output results to a text file
    printf "%s\\t%d\\t%d\\t%.4f\\n" "${meta.id}" \$peak_read \$total_read \$score > ${meta.id}_FRiP_score.txt
    """
}