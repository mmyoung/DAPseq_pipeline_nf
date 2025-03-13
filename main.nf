nextflow.enable.dsl=2

//params.help = false
params.threads = 1
params.output_dir = './results'

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to call DAP-seq peaks with MACS3 and peak annotation and motif analyses'
    log.info '--------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run ./ -params-file params.yml'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--prime5_trim_len    Trim the first ? number of bases in 5-end of each read.'
    log.info '	--prime3_trim_len    Trim the first ? number of bases in 3-end of each read.'
    log.info '  --bowtie_idx bowtie_idx    bowtie index prefix for analysis.'
    log.info '  --data_dir raw_data    The folder where the raw .fq files are.'
    log.info '	--fq_sheet fq_sheet    A tab-delimited file storing the samples information, with three columns: sample_id,fq1,fq2,single_end,control.'
    log.info '	--fasta    Genome fasta file for the analyzing species.'
    log.info '	--gtf    Genome gtf file for the analyzing species.'
    log.info '  --output_dir OUTDIR   Name for directory for saving the results. Default: results/'
    log.info '  --fq_dir raw_data    The folder where the raw .fq files are.'
    log.info '  --gsize The size of analyzing genome.'
    exit 1
}


include {INPUT_CHECK} from "./subworkflow/input_check"
include {TRIMGALORE} from "./module/Trimgalore"
include {FASTQC} from "./module/Fastqc"
include {BAM2BW} from "./module/Bam2bw"
include {BOWTIE2MAP} from "./module/Bowtie2Mapping"
include {COVERAGE} from "./module/CoverageCal"
include {MARK_DUPLICATES} from './module/markduplicates'
include {MACS2_CALLPEAK} from "./module/MACS2"

include {HOMER_ANNOTATEPEAKS} from "./module/homer_annopeak"
include {MEME_MOTIF} from "./module/MEME_motif"
include {CAL_FRIP} from "./module/calculate_FRiP"

// Define a process to collect completion signals
process COMPLETION_CHECK {
    input:
        val(count)
        val(output_dir)
    
    output:
        val(true), emit: ready
    
    exec:
        println "All $count samples have completed processing"
}

process SAMPLE_REPORTING{

    conda  "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"
    tag "report"
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
        path(sample_sheet)
        val(output_dir)
        val(ready_signal)

    output:
        path("read_peak.num.summary")
        path("report.html")

    script:

        def reportPath = "${projectDir}/bin/report.Rmd"

        """
        sed "1d" ${sample_sheet} | while IFS=',' read ID fq1 _ _ _ || [[ -n "\${ID}" ]] ;do
            raw_num=`cat ${output_dir}/trimm/\${fq1}_trimming_report.txt |grep "Total reads processed:"|sed s/" "//g|cut -d ":" -f 2`
            peak_num=`cat ${output_dir}/macs3_output/\${ID}_peaks.narrowPeak | wc -l `

            mapped_reads=`cat ${output_dir}/FRiP_score/\${ID}_FRiP_score.txt | cut -f 3`
            FRiP_score=`cat ${output_dir}/FRiP_score/\${ID}_FRiP_score.txt | cut -f 4`

            printf "\${ID}\t\${raw_num}\t\${mapped_reads}\t\${peak_num}\t\${FRiP_score}\n"
        done >read_peak.num.summary

        mkdir -p ${params.output_dir}/report/
        cp read_peak.num.summary ${params.output_dir}/report/

        Rscript -e "rmarkdown::render(input = '${reportPath}', 
                                    output_file='report.html', 
                                    params=list(summary_table='${params.output_dir}/report/read_peak.num.summary',parent_path='${output_dir}'))"

        cp report.html ${params.output_dir}/report/
        """

}


workflow {

    INPUT_CHECK(params.fq_sheet)
    
    TRIMGALORE(INPUT_CHECK.out.reads)
    BOWTIE2MAP(TRIMGALORE.out.reads)
//  MARK_DUPLICATES(BOWTIE2MAP.out.bam)
    
    FASTQC(TRIMGALORE.out.reads)
//  MARK_DUPLICATES
    BOWTIE2MAP
        .out
        .bam
//      .join(MARK_DUPLICATES.out.bai, by: [0])
        .join(BOWTIE2MAP.out.bai, by: [0])
        .set {ch_genome_bam_bai}

//    ch_genome_bam_bai|(BAM2BW & COVERAGE)

    BAM2BW(ch_genome_bam_bai)
    COVERAGE(ch_genome_bam_bai)

//    ch_genome_bam_bai
//       .combine(ch_genome_bam_bai)
//        .map { 
//            meta1, bam1, bai1, meta2, bam2, bai2 ->
//                meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ], [ bai1, bai2 ] ] : [ meta1, [ bam1, null ], [ bai1, null ] ] 
//        }
//        .set { ch_ip_control_bam_bai }

//    ch_genome_bam_bai.view()

    ch_genome_bam_bai
        .map { 
            meta, bam, bai -> 
                [ meta , bam, bai ] 
        }
        .set { ch_ip_control_bam }

//    ch_ip_control_bam.view()

// call peaks with macs
    MACS2_CALLPEAK(ch_ip_control_bam, params.gsize)
// calculate FRiP score based on the peak coverage 
    CAL_FRIP(MACS2_CALLPEAK.out.peak,BOWTIE2MAP.out.bam)
// annotate peak distribution using HOMER
    HOMER_ANNOTATEPEAKS(MACS2_CALLPEAK.out.peak, params.fasta, params.gtf)
// analyze motif of peaks using MEME suite
    MEME_MOTIF(MACS2_CALLPEAK.out.peak, params.fasta)
// Count the number of samples and wait for all processes to complete
    CAL_FRIP.out.txt
        .count()
        .set { sample_count }
    
    // Create a completion check that waits for all processes
    COMPLETION_CHECK(sample_count, params.output_dir)
    
    // Run the reporting at the end of the main workflow after completion check
    SAMPLE_REPORTING(
        params.fq_sheet,
        params.output_dir,
        COMPLETION_CHECK.out.ready
    )
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

