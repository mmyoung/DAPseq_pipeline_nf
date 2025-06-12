process COVERAGE {

    // conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"
    tag "Calculate the reads coverage in each base ... "
    publishDir "${params.output_dir}/coverage_out", mode: 'copy'

    input:
        tuple val(meta), path(bam_file), path(bai)

    output:
        tuple path("${meta.id}_base.depth"), path("${meta.id}_depth.pdf")

    script: 
        """
        module load samtools/1.13
    
        samtools depth -a ${bam_file} | cut -f 3 | sort | uniq -c > ${meta.id}_base.depth
        Rscript ${projectDir}/bin/line_plot.r ${meta.id}_base.depth ${meta.id} ${meta.id}_depth.pdf
        """
}