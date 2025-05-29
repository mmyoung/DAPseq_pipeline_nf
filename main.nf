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
            min5fold_peak_num=`cat ${output_dir}/macs3_output/\${ID}_peaks.narrowPeak | awk '\$7>5{print \$0}'| wc -l `
            max_peak_foldch=`cat ${output_dir}/macs3_output/\${ID}_peaks.narrowPeak | awk '\$7>max{max=\$7}END{print max}'`
            mapping_ratio=`cat ${output_dir}/trimm/\${fq1}_trimming_report.txt | grep "overall alignment rate"| sed s/"overall alignment rate"//g`

            mapped_reads=`cat ${output_dir}/FRiP_score/\${ID}_FRiP_score.txt | cut -f 3`
            peak_reads=`cat ${output_dir}/FRiP_score/\${ID}_FRiP_score.txt | cut -f 2`
            FRiP_score=`cat ${output_dir}/FRiP_score/\${ID}_FRiP_score.txt | cut -f 4`

            printf "\${ID}\t\${raw_num}\t\${mapped_reads}\t\${mapping_ratio}\t\${peak_num}\t\${min5fold_peak_num}\t\${max_peak_foldch}\t\${peak_reads}\t\${FRiP_score}\n"
        done >read_peak.num.summary

        mkdir -p ${params.output_dir}/report/
        cp read_peak.num.summary ${params.output_dir}/report/

        Rscript -e "rmarkdown::render(input = '${reportPath}', 
                                    output_file='report.html', 
                                    params=list(summary_table='${params.output_dir}/report/read_peak.num.summary',parent_path='${output_dir}'))"

        mv ${projectDir}/bin/report.html ${params.output_dir}/report/
        """
}

workflow {
    INPUT_CHECK(params.fq_sheet)
    
    TRIMGALORE(INPUT_CHECK.out.reads)
    BOWTIE2MAP(TRIMGALORE.out.reads)
    FASTQC(TRIMGALORE.out.reads)
    
    BOWTIE2MAP
        .out
        .bam
        .join(BOWTIE2MAP.out.bai, by: [0])
        .set {ch_genome_bam_bai}

    BAM2BW(ch_genome_bam_bai)
    COVERAGE(ch_genome_bam_bai)

    // CORRECTED LOGIC FOR YOUR SAMPLE STRUCTURE
    
    // Identify control samples (these are the samples that serve as controls for others)
    // Based on your CSV: the last 3 samples that have empty control field but are referenced by others
    def control_sample_names = [
        'multiDAP_At_seedling_shoot_Halo_beads',
        'multiDAP_At_seedling_root_Halo_beads', 
        'multi_ampDAP_At_Halo_beads'
    ]
    
    // Split all BAM samples into controls and treatments
    ch_genome_bam_bai
        .branch { meta, bam, bai ->
            // Control samples: those that are referenced as controls by other samples
            controls: control_sample_names.contains(meta.id)
            // Treatment samples: all others (both with and without controls)
            treatments: !control_sample_names.contains(meta.id)
        }
        .set { sample_branches }

    // Debug: Check the split
    sample_branches.controls.count().view { "Control samples: $it" }
    sample_branches.treatments.count().view { "Treatment samples: $it" }

    // Prepare available controls (keyed by sample ID)
    sample_branches.controls
        .map { meta, bam, bai -> [ meta.id, bam, bai ] }
        .set { available_controls }

    available_controls.view { "Available control: $it" }

    // Process treatment samples
    sample_branches.treatments
        .branch { meta, bam, bai ->
            // Samples that need a control (have non-empty control field)
            with_control: meta.control && meta.control != ''
            // Samples that don't need a control (empty control field)
            without_control: !meta.control || meta.control == ''
        }
        .set { treatment_branches }

    // Debug: Check treatment split
    treatment_branches.with_control.count().view { "Treatments with controls: $it" }
    treatment_branches.without_control.count().view { "Treatments without controls: $it" }

    // For treatments WITH controls: join with their specified control
    treatment_branches.with_control
        .map { meta, bam, bai -> [ meta.control, meta, bam, bai ] }  // Key by required control ID
        .join(available_controls, by: 0, remainder: true)  // Join with available controls
        .map { control_id, meta, treatment_bam, treatment_bai, control_bam, control_bai ->
            [ meta, treatment_bam, control_bam ?: null ]
        }
        .set { treatments_with_controls }

    // For treatments WITHOUT controls: use null as control
    treatment_branches.without_control
        .map { meta, bam, bai -> [ meta, bam, null ] }
        .set { treatments_without_controls }

    // Combine both types of treatments
    treatments_with_controls.mix(treatments_without_controls)
        .set { ch_ip_control_bam }

    // Debug: Check final MACS2 input
    ch_ip_control_bam.count().view { "Total samples for MACS2: $it" }
    ch_ip_control_bam.view { "MACS2 input: ${it[0].id} (control: ${it[2] ? 'YES' : 'NO'})" }

    // Call peaks with MACS - should now process ~61 samples (22 without controls + 39 with controls)
    MACS2_CALLPEAK(ch_ip_control_bam, params.gsize)
    
    // Remove any duplicates from peak channel (precautionary)
    MACS2_CALLPEAK.out.peak
        .unique { sample_id, peak_file -> sample_id }
        .set { unique_peak_channel }
    
    // Debug: Check final peak channel
    unique_peak_channel.count().view { "Total unique peaks: $it" }
    
    // Use the cleaned channel for all downstream processes
    CAL_FRIP(unique_peak_channel, BOWTIE2MAP.out.bam)
    HOMER_ANNOTATEPEAKS(unique_peak_channel, params.fasta, params.gtf)
    MEME_MOTIF(unique_peak_channel, params.fasta)
    
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