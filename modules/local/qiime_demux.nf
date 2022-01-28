process QIIME_DEMUX {
    tag "$meta.id"

    conda (params.enable_conda ? "qiime2::qiime2=2021.11 qiime2::q2cli=2021.11 qiime2::q2-demux=2021.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qiime:1.9.1--py_3' :
        'quay.io/qiime2/core:2021.11' }"

    input:
    tuple val(meta), path(qza)
    path barcodes

    output:
    tuple val(meta), path("${prefix}_demux-full.qza")   , emit: reads
    tuple val(meta), path("${prefix}_demux-details.qza"), emit: details
    tuple val(meta), path("${prefix}_demux-summary.qzv"), emit: summary
    path "versions.yml", emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    qiime demux emp-paired \\
        --m-barcodes-file $barcodes \\
        --i-seqs $qza \\
        --o-per-sample-sequences ${prefix}_demux-full.qza \\
        --o-error-correction-details ${prefix}_demux-details.qza \\
        $args

    qiime demux summarize \\
        --i-data ${prefix}_demux-full.qza \\
        --o-visualization ${prefix}_demux-summary.qzv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
