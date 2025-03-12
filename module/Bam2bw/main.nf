process BAM2BW {
    tag "Making .bw files from .bam files"

    conda  "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"
    
    publishDir "${params.output_dir}/bw_output", mode: 'copy'

    input:
        tuple val(meta), path(bam_file), path(bai)

    output:
        path("${meta.id}_sorted_bam.bw")

    script: 
        """
        bamCoverage -b ${bam_file} -o ${meta.id}_sorted_bam.bw
        """
}