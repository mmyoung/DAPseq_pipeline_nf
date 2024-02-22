process TRIMGALORE {
    tag "Adapter and low-quality based trimming"
    publishDir "${params.output_dir}/trimm", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*{trimmed,val}*.fq.gz"), emit: reads
        tuple val(meta), path("*report.txt")          , emit: log
        tuple val(meta), path("*.html")               , emit: html    , optional: true
        tuple val(meta), path("*.zip")                , emit: zip     , optional: true

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Clipping presets have to be evaluated in the context of SE/PE
    def c_r1   = params.prime5_trim_len > 0             ? "--clip_r1 ${params.prime5_trim_len}"                         : ''
    def c_r2   = params.prime5_trim_len > 0             ? "--clip_r2 ${params.prime5_trim_len}"                         : ''
    def tpc_r1 = params.prime3_trim_len > 0 ? "--three_prime_clip_r1 ${params.prime3_trim_len}"             : ''
    def tpc_r2 = params.prime3_trim_len > 0 ? "--three_prime_clip_r2 ${params.prime3_trim_len}"             : ''

    def prefix = "${meta.id}"
    print"${meta.single_end}"
    if (meta.single_end) {
        """
        /project/zhuzhuzhang/lyang/software/TrimGalore-0.6.10/trim_galore \\
        --path_to_cutadapt /project/zhuzhuzhang/lyang/software/miniconda3/envs/cutadaptenv/bin/cutadapt \\
        --cores $cores \\
        --gzip \\
        --basename $prefix \\
        $c_r1 \\
        $tpc_r1 \\
        $reads
        """
    } else {
        """
        /project/zhuzhuzhang/lyang/software/TrimGalore-0.6.10/trim_galore \\
        --path_to_cutadapt /project/zhuzhuzhang/lyang/software/miniconda3/envs/cutadaptenv/bin/cutadapt \\
            --cores $cores \\
            --paired \\
            --gzip \\
            --basename $prefix \\
            $c_r1 \\
            $c_r2 \\
            $tpc_r1 \\
            $tpc_r2 \\
            ${reads[0]} \\
            ${reads[1]}
        """
    }
}