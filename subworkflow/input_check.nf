workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel.fromPath(params.fq_sheet, checkIfExists: true)
            .ifEmpty { exit 1, "bam sheet not found" }
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { reads }

    emit:
    reads      // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.control = row.control


    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(params.fq_dir + "/" + row.fq1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fq1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(params.fq_dir + "/" + row.fq1) ] ]
    } else {
        if (!file(params.fq_dir + "/" + row.fq2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fq2}"
        }
        fastq_meta = [ meta, [ file(params.fq_dir + "/" + row.fq1), file(params.fq_dir + "/" + row.fq2) ] ]
    }
    return fastq_meta
}