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
    log.info 'nextflow run ./ -params-file *'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '	--sample_sheet sample_sheet    A tab-delimited file storing the samples information, with three columns: sample, ipbam, ctrlbam.'
    log.info '	--fasta    Genome fasta file for the analyzing species.'
    log.info '	--gtf    Genome gtf file for the analyzing species.'
    log.info '  --output_dir OUTDIR   Name for directory for saving the results. Default: results/'
    log.info '  --data_dir raw_data    The folder where the raw .bam files are.'
    log.info '  --gsize The size of analyzing genome.'
    exit 1
}

include {MACS2_CALLPEAK} from "./module/MACS2"
include {HOMER_ANNOTATEPEAKS} from "./module/homer_annopeak"
include {MEME_MOTIF} from "./module/MEME_motif"


sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: ',') 
workflow {

    sample_sheet.map { row -> row.ipbam = params.raw_bam_dir + "/" + row.ipbam; row }
                .map { row -> row.ctrlbam = params.raw_bam_dir + "/" + row.ctrlbam; row }
                .map { row -> [row.sample, file(row.ipbam, checkIfExists: true), file(row.ctrlbam, checkIfExists: true)] }.set{ all_in }
    all_in.view()
    MACS2_CALLPEAK(all_in, params.gsize)
    HOMER_ANNOTATEPEAKS(MACS2_CALLPEAK.out.peak, params.fasta, params.gtf)
    MEME_MOTIF(MACS2_CALLPEAK.out.peak, params.fasta)
}