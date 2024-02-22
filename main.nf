nextflow.enable.dsl=2

params.help = false
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



workflow {

    INPUT_CHECK(params.fq_sheet)
    
    TRIMGALORE(INPUT_CHECK.out.reads)
    BOWTIE2MAP(TRIMGALORE.out.reads)
    MARK_DUPLICATES(BOWTIE2MAP.out.bam)
    
    FASTQC(TRIMGALORE.out.reads)
    MARK_DUPLICATES
        .out
        .bam
        .join(MARK_DUPLICATES.out.bai, by: [0])
        .set {ch_genome_bam_bai}

    ch_genome_bam_bai|(BAM2BW & COVERAGE)

    ch_genome_bam_bai
        .combine(ch_genome_bam_bai)
        .map { 
            meta1, bam1, bai1, meta2, bam2, bai2 ->
                meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ], [ bai1, bai2 ] ] : null
        }
        .set { ch_ip_control_bam_bai }

    ch_ip_control_bam_bai
    .map { 
        meta, bams, bais -> 
            [ meta , bams[0], bams[1] ] 
    }
    .set { ch_ip_control_bam }

    MACS2_CALLPEAK(ch_ip_control_bam, params.gsize)

    HOMER_ANNOTATEPEAKS(MACS2_CALLPEAK.out.peak, params.fasta, params.gtf)
    MEME_MOTIF(MACS2_CALLPEAK.out.peak, params.fasta)
    
}