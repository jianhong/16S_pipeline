process QIIME_IMPORT {
    tag "$meta.id"

    conda (params.enable_conda ? "qiime2::qiime2=2021.11 qiime2::q2cli=2021.11 qiime2::q2-demux=2021.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://martinjf/default/qiime2:2021.8' :
        'quay.io/qiime2/core:2021.11' }"

    input:
    tuple val(meta), path(reads1, stageAs: "sync/forward.fastq.gz"), path(reads2, stageAs: "sync/reverse.fastq.gz"), path(index, stageAs: "sync/barcodes.fastq.gz")

    output:
    tuple val(meta), path("${prefix}_emp-paired-end-sequences.qza")       , emit: reads
    path "versions.yml", emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    qiime tools import \\
        --input-path sync \\
        --output-path ${prefix}_emp-paired-end-sequences.qza \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(qiime info | sed -n '/: [0-9.]/p' | sed 's/.*/    &/g')
    END_VERSIONS
    """
}
