process MEME_MOTIF {
    conda '/project/gzy8899/lyang/DAPseq_pipeline_nf/env_requirment.yml'

    tag "MEME motif analysis"
    label 'Peak_annotation'
    publishDir "${params.output_dir}/meme_output", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

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