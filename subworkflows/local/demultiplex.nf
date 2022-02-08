//
// Check input samplesheet and get read channels
//

include { QIIME_IMPORT } from '../../modules/local/qiime_import'
include { QIIME_DEMUX  } from '../../modules/local/qiime_demux'
include { QIIME_EXPORT } from '../../modules/local/qiime_export'

workflow QIIME_DEMULTIPLEX {
    take:
    reads // [meta, R1, R2, I1]
    barcodes
    single_end

    main:
    ch_versions = QIIME_IMPORT ( reads, single_end ).versions
    QIIME_DEMUX (QIIME_IMPORT.out.reads, barcodes, single_end)
    QIIME_EXPORT (QIIME_DEMUX.out.reads)

    emit:
    reads    = QIIME_EXPORT.out.reads         // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}
