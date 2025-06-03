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
    log.info '  --control_samples    Comma-separated list of control sample names (optional, auto-detected if not provided)'
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
            echo -e "Sample_ID\\tRaw_reads\\tMapped_reads\\tMapping_ratio\\tPeak_number\\tMin5fold_peaks\\tMax_peak_foldchange\\tPeak_reads\\tFRiP_score" > read_peak.num.summary
            
            sed "1d" ${sample_sheet} | while IFS=',' read ID fq1 _ _ _ || [ -n "\$ID" ]; do
                [ -z "\$ID" ] && continue
                
                # Check if peak file exists
                peak_file="${output_dir}/macs3_output/\$ID"_peaks.narrowPeak""
                frip_file="${output_dir}/FRiP_score/\$ID"_FRiP_score.txt""
                
                if [ -f "\$peak_file" ] && [ -f "\$frip_file" ]; then
                    
                    raw_num=\$(cat ${output_dir}/trimm/\$fq1"_trimming_report.txt" | grep "Total reads processed:" | sed 's/ //g' | cut -d ":" -f 2 || echo "NA")
                    peak_num=\$(cat "\$peak_file" | wc -l)
                    min5fold_peak_num=\$(awk '\$7>5' "\$peak_file" | wc -l)
                    max_peak_foldch=\$(awk '\$7>max{max=\$7}END{print max}' "\$peak_file")
                    mapping_ratio=\$(cat ${output_dir}/alignment/\$ID".bowtie2.log" | grep "overall alignment rate" | sed 's/overall alignment rate//g' || echo "NA")
                    
                    mapped_reads=\$(cut -f 3 "\$frip_file")
                    peak_reads=\$(cut -f 2 "\$frip_file")
                    FRiP_score=\$(cut -f 4 "\$frip_file")
                    
                    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "\$ID" "\$raw_num" "\$mapped_reads" "\$mapping_ratio" "\$peak_num" "\$min5fold_peak_num" "\$max_peak_foldch" "\$peak_reads" "\$FRiP_score"
                fi
            done >> read_peak.num.summary

            mkdir -p ${params.output_dir}/report/
            cp read_peak.num.summary ${params.output_dir}/report/

            Rscript -e "rmarkdown::render(input = '${reportPath}', output_file='report.html', params=list(summary_table='${params.output_dir}/report/read_peak.num.summary',parent_path='${output_dir}'))"

            mv ${projectDir}/bin/report.html ${params.output_dir}/report/
        """
}

// SIMPLIFIED AND CORRECTED APPROACH:

workflow {
    INPUT_CHECK(params.fq_sheet)
    
    INPUT_CHECK.out.reads.count().view { "=== TOTAL INPUT SAMPLES: $it ===" }
    
    TRIMGALORE(INPUT_CHECK.out.reads)
    BOWTIE2MAP(TRIMGALORE.out.reads)
    FASTQC(TRIMGALORE.out.reads)
    
    BOWTIE2MAP
        .out
        .bam
        .join(BOWTIE2MAP.out.bai, by: [0])
        .set {ch_genome_bam_bai}

    ch_genome_bam_bai.count().view { "=== TOTAL BAM SAMPLES: $it ===" }

    BAM2BW(ch_genome_bam_bai)
    COVERAGE(ch_genome_bam_bai)

        // FIXED CONTROL DETECTION - Simple and reliable

    def known_controls = params.control_samples ? 
    params.control_samples.split(',').collect { it.trim() } : 
    []

    log.info "Using known control samples: ${known_controls}"

    // Simple branching - no complex channel operations
    ch_genome_bam_bai
        .branch { meta, bam, bai ->
            controls: known_controls.contains(meta.id)
            treatments: !known_controls.contains(meta.id)
        }
        .set { sample_branches }

        // Debug the basic split
    sample_branches.controls.count().view { "=== CONTROL SAMPLES: $it ===" }
    sample_branches.treatments.count().view { "=== TREATMENT SAMPLES: $it ===" }

    // SIMPLE CONTROL LOOKUP - No complex operations
    sample_branches.controls
        .map { meta, bam, bai -> 
            println "CONTROL AVAILABLE: ${meta.id} -> ${bam}"
            return [ meta.id, [bam], [bai] ]
        }
        .set { control_lookup }

    // Process treatments
    sample_branches.treatments
        .branch { meta, bam, bai ->
            needs_control: meta.control && meta.control.trim() != ''
            no_control: !meta.control || meta.control.trim() == ''
        }
        .set { treatment_branches }

    treatment_branches.needs_control.count().view { "=== TREATMENTS NEEDING CONTROLS: $it ===" }
    treatment_branches.no_control.count().view { "=== TREATMENTS WITHOUT CONTROLS: $it ===" }

    // Perform the join
    treatment_branches.needs_control
        .map { meta, bam, bai ->
                meta.control ? [ meta.control, meta, [ bam ], [ bai ] ] : null}
        .combine(control_lookup, by: 0)
        .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
        .set { treatments_with_controls }

    // For treatments without controls
    treatment_branches.no_control
        .map { meta, bam, bai -> 
            return [meta, [bam, null], [bai, null]] 
        }
        .set { treatments_without_controls }

    // Combine all treatments for MACS2
    treatments_with_controls
        .mix(treatments_without_controls)
        .map {
            meta, bams, bais ->
                [ meta , bams[0], bams[1] ]
        }
        .set { macs_input_ready }

    // Final debug
    macs_input_ready.count().view { "=== TOTAL FOR MACS2: $it ===" }


    // MACS2 and downstream
    //MACS2_CALLPEAK(all_treatments_for_macs, params.gsize)
    MACS2_CALLPEAK(macs_input_ready, params.gsize)

    // Extract the ID from the meta map to match MACS2's structure
    bowtie_bam_for_frip = BOWTIE2MAP.out.bam
        .map { meta, bam -> [meta.id, bam] }

    MACS2_CALLPEAK.out.peak.count().view { "=== TOTAL PEAK FILES: $it ===" }
    bowtie_bam_for_frip.count().view { "=== TOTAL BAM FILES: $it ===" }

    // Now they both have [sample_id, file] structure
    cal_frip_input = MACS2_CALLPEAK.out.peak.join(bowtie_bam_for_frip, by: 0)
    cal_frip_input.out.txt.count().view { "=== TOTAL FRIP SCORE FILES: $it ===" }
    CAL_FRIP( cal_frip_input )

    CAL_FRIP.out.txt.count().set { sample_count }

    // annotate peak distribution using HOMER
    HOMER_ANNOTATEPEAKS(MACS2_CALLPEAK.out.peak, params.fasta, params.gtf)
    // analyze motif of peaks using MEME suite
    MEME_MOTIF(MACS2_CALLPEAK.out.peak, params.fasta)


    COMPLETION_CHECK(sample_count, params.output_dir)
    SAMPLE_REPORTING(params.fq_sheet, params.output_dir, COMPLETION_CHECK.out.ready)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}