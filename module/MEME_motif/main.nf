process MEME_MOTIF {
    tag "$sample_id"
    label 'Peak_annotation'
    publishDir "${params.output_dir}/meme_output", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    tuple val(sample_id), path(peak)
    path  fasta

    output:
    path("*"), emit: txt

    script:
    """
    /project/zhuzhuzhang/lyang/software/bedtools2/bin/bedtools getfasta -fi ${fasta} -bed ${peak} -fo ${sample_id}.peak.fasta
    meme ${sample_id}.peak.fasta -dna -o ${sample_id}_meme
    """
}