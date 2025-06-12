process CAL_FRIP {

    // conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    tag "FRiP calculation"
    label 'FRiP'
    publishDir "${params.output_dir}/FRiP_score", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    input:
    tuple val(meta), path(peak), path(bam)

    output:
    path("*.txt"), optional: true, emit: txt

    shell:

    """
    module load samtools

    if [ -s ${peak} ]; then

        # Count reads in the peaks
        peak_read=\$(samtools view -c -L <(cat ${peak}|awk '\$7>5{print \$0}') ${bam})

        # Count total reads in the BAM file
        total_read=\$(samtools view -c ${bam})

        # Calculate the FRiP score
        score=\$(awk "BEGIN {print \$peak_read / \$total_read}")

        # Output results to a text file
        printf "%s\\t%d\\t%d\\t%.4f\\n" "${meta}" \$peak_read \$total_read \$score > ${meta}_FRiP_score.txt

    fi
    """
}