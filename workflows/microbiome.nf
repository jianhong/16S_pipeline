/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMicrobiome.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.samplesheet, params.metadata ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input bcl folder not specified!' }
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Input samplesheet not specified!' }
if (!params.skip_demultiplex) {
    if (params.barcodes) { ch_barcodes = file(params.barcodes, checkIfExists: true) } else { exit 1, 'Barcodes not specified!' }
}
if (params.metadata) { ch_metadata = file(params.metadata, checkIfExists: true) } else { exit 1, 'Metadata not specified!' }
if (params.silva_nr99) { ch_silva_nr99 = file(params.silva_nr99, checkIfExists: true) } else { exit 1, 'Training set not specified!' }
if (params.silva_tax) { ch_silva_tax = file(params.silva_tax, checkIfExists: true) } else { exit 1, 'Species assignment data not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { QIIME_DEMULTIPLEX } from '../subworkflows/local/demultiplex'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { BCL2FASTQ                   } from '../modules/local/bcl2fastq'
include { REMOVE_PRIMERS              } from '../modules/local/remove_primers'
include { SYNC_BARCODES               } from '../modules/local/sync_barcodes'
include { FILTERING                   } from '../modules/local/filtering'
include { DADA2                       } from '../modules/local/dada2'
include { PHYLOSEQ                    } from '../modules/local/phyloseq'
include { KRONA                       } from '../modules/local/krona'
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MICROBIOME {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))

    //
    // MODULE: Convert .bcl (basecall) files to .fastq files.
    //
    if(!params.skip_bcl2fastq){
        BCL2FASTQ(ch_input, ch_samplesheet, params.bcl2fastq)
        ch_versions = ch_versions.mix(BCL2FASTQ.out.versions)
        bcl_undetermined = BCL2FASTQ.out.reads
    }else{
        bcl_undetermined = Channel.fromPath("${params.input}/**.{fastq,fq}.gz")
    }
    if(params.verbose){
        bcl_undetermined.view{'bcl_undetermined:'+it.baseName}
    }

    if(!params.skip_demultiplex){
        //
        // convert reads files to format [meta, [reads1, reads2]]
        //
        if(params.single_end){
            ch_reads = bcl_undetermined.flatten()
                                .filter{ it.baseName =~ params.R1_pattern }
                                .map{[[id: it.baseName.replaceAll(params.basename_pattern, ''), single_end: true], it]}
        }else{
            ch_reads = bcl_undetermined.flatten()
                                .filter{ it.baseName =~ params.R12_pattern }
                                .map{[[id: it.baseName.replaceAll(params.basename_pattern, ''), single_end: false], it ]}
                                .groupTuple(sort: true)
        }
        if(params.verbose) {
            ch_reads.view{'ch_reads:'+it}
        }

        //
        // MODULE: Run FastQC
        //
        if(!params.skip_fastqc){
            FASTQC (
                ch_reads
            )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: trimmomatic to remove primers from raw reads
        //
        REMOVE_PRIMERS(ch_reads, params.single_end)
        ch_versions = ch_versions.mix(REMOVE_PRIMERS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(REMOVE_PRIMERS.out.log.collect().ifEmpty([]))

        //
        // MODULE: Make sure forward, reverse, and index read files are all aligned
        //
        if(!params.single_end){
            ch_reads2 = bcl_undetermined.flatten()
                                .filter{ it.baseName =~ params.I1_pattern }
                                .map{[[id: it.baseName.replaceAll(params.basename_pattern, ''), single_end: false], it]}
                                .join(REMOVE_PRIMERS.out.paired)
            //ch_reads2.view()
            SYNC_BARCODES(ch_reads2)
            ch_versions = ch_versions.mix(SYNC_BARCODES.out.versions)
            ch_reads3 = SYNC_BARCODES.out.reads
        }else{
            ch_reads3 = REMOVE_PRIMERS.out.paired.map{[[id:it[0].id], it[1], []]}
                        .join(bcl_undetermined.flatten()
                            .filter{ it.baseName =~ params.I1_pattern }
                            .map{[[id: it.baseName.replaceAll(params.basename_pattern, '')], it]})
        }

        if(params.verbose) {
            ch_reads3.view{'ch_reads3:'+it}
        }

        //
        // MODULE: demultiplex
        //
        QIIME_DEMULTIPLEX(ch_reads3, ch_barcodes, params.single_end)
        ch_versions = ch_versions.mix(QIIME_DEMULTIPLEX.out.versions)
        ch_reads4 = QIIME_DEMULTIPLEX.out.reads
    }else{
        ch_reads = bcl_undetermined.flatten()
                            .filter{ it.baseName =~ params.R12_pattern }
                            .map{[[id: it.baseName.replaceAll(params.basename_pattern, ''), single_end: params.single_end], it]}
                            .groupTuple()

        if(params.verbose){
            ch_reads.view{'ch_reads:'+it}
        }
        //
        // MODULE: Run FastQC
        //
        if(!params.skip_fastqc){
            FASTQC (
                ch_reads
            )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: trimmomatic to remove primers
        //
        REMOVE_PRIMERS(ch_reads, params.single_end)
        ch_versions = ch_versions.mix(REMOVE_PRIMERS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(REMOVE_PRIMERS.out.log.collect().ifEmpty([]))

        ch_reads2 = REMOVE_PRIMERS.out.paired
        if(params.single_end){
            ch_reads3 = ch_reads2.map{[[id:'demultiplexed'], it[1]]}
        }else{
            ch_reads3 = ch_reads2.map{[it[1], it[2]]}.flatten().map{[[id:'demultiplexed'], it]}
        }
        ch_reads4 = ch_reads3.groupTuple()
    }

    //
    // MODULE: Filter
    //
    FILTERING(ch_reads4, params.single_end)
    ch_versions = ch_versions.mix(FILTERING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FILTERING.out.qc.collect().ifEmpty([]))

    //
    // MODULE: Run dada2
    //
    DADA2(FILTERING.out.reads, ch_silva_nr99, ch_silva_tax)
    ch_versions = ch_versions.mix(DADA2.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DADA2.out.qc.collect().ifEmpty([]))

    //
    // MODULE: Run phyloseq
    //
    PHYLOSEQ(DADA2.out.robj, ch_metadata)
    ch_versions = ch_versions.mix(PHYLOSEQ.out.versions)

    //
    // MODULE: Run Krona
    //
    KRONA(PHYLOSEQ.out.krona)
    ch_versions = ch_versions.mix(KRONA.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMicrobiome.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
