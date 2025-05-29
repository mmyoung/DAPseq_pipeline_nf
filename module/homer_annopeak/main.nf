process HOMER_ANNOTATEPEAKS {
    tag "Homer peak annotation for ${sample_id}"
    label 'Peak_annotation'
    publishDir "${params.output_dir}/homer_annotation", mode: 'copy'

    conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    tuple val(sample_id), path(peak, stageAs: 'peaks_homer.narrowPeak')
    path fasta, stageAs: 'genome_homer.fa' 
    path gtf, stageAs: 'annotation.gtf'

    output:
    tuple val(sample_id), path("*annotatePeaks.txt"), optional: true, emit: txt
    tuple val(sample_id), path("*.annStats.txt"), optional: true, emit: stats

    script:
    """
    if [ -s peaks_homer.narrowPeak ]; then
        echo "Processing peak file: peaks_homer.narrowPeak"
        echo "Peak file size: \$(wc -l < peaks_homer.narrowPeak) lines"
        
        if [[ annotation.gtf == *.gtf ]]; then
            echo "Using GTF format annotation"
            annotatePeaks.pl \\
                peaks_homer.narrowPeak \\
                genome_homer.fa \\
                -gtf annotation.gtf \\
                -gid \\
                -cpu ${task.cpus ?: 1} \\
                -annStats ${sample_id}.annStats.txt \\
                > ${sample_id}.annotatePeaks.txt
        else
            echo "Using GFF format annotation" 
            annotatePeaks.pl \\
                peaks_homer.narrowPeak \\
                genome_homer.fa \\
                -gff annotation.gtf \\
                -keepAll \\
                -cpu ${task.cpus ?: 1} \\
                -annStats ${sample_id}.annStats.txt \\
                > ${sample_id}.annotatePeaks.txt
        fi
        
        echo "Annotation completed. Output file size: \$(wc -l < ${sample_id}.annotatePeaks.txt) lines"
    else
        echo "Peak file peaks_homer.narrowPeak is empty or does not exist"
        touch ${sample_id}.annotatePeaks.txt
        touch ${sample_id}.annStats.txt
    fi
    """
}