process MEME_MOTIF {
    
    // conda '/project/gzy8899/lyang/DAPseq_pipeline_nf/env_export.yml'
    conda  "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"
    
    tag "MEME motif analysis"
    label 'Peak_annotation'
    publishDir "${params.output_dir}/meme_output", mode: 'copy'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    // Using the top 100 peaks for motif analysis

    input:
    tuple val(sample_id), path(peak)
    path  fasta

    output:
    path("*"), optional: true,emit: txt

    shell:
    """
    if [ -s ${peak} ]; then

    /project/zhuzhuzhang/lyang/software/bedtools2/bin/bedtools getfasta -fi ${fasta} -bed <(sort -k 7 -n -r ${peak}|head -100|awk 'BEGIN{OFS="\t"}{print \$1,(\$2+\$10-30 >0)?\$2+\$10-30:1,\$2+\$10+30}') -fo ${sample_id}.peak.fasta
    meme ${sample_id}.peak.fasta -dna -revcomp -mod anr -nmotifs 2 -minw 8 -maxw 32 -maxsize 10000000 -o ${sample_id}_meme

    fi
    """
}