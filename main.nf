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
    log.info '    nextflow run ./ -params-file *'
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

include {CALL_MACS} from "./module/MACS2"
include {ANNO_PEAK} from "./module/homer_anopeak"
include {CALL_MEME} from "./module/MEME_motif"


sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: ',') 

workflow {

    aln_in = sample_sheet.map { row -> row.ipbam = params.data_dir + "/" + row.ipbam; row }
                .map { row -> row.ctrlbam = params.data_dir + "/" + row.ctrlbam; row }
                .map { row -> [row.sample, file(row.ipbam), file(row.ctrlbam)] }
    CALL_MACS(aln_in, params.gsize)
    ANNO_PEAK(CALL_MACS.out.peak, params.fasta, params.gtf)
    CALL_MEME(CALL_PEAK.out.peak)
}


