process QIIME_DEMUX {
    tag "$meta.id"

    conda (params.enable_conda ? "qiime2::qiime2=2021.11 qiime2::q2cli=2021.11 qiime2::q2-demux=2021.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://martinjf/default/qiime2:2021.8' :
        'quay.io/qiime2/core:2021.11' }"

    input:
    tuple val(meta), path(qza)
    path barcodes
    val single_end

    output:
    tuple val(meta), path("${prefix}_demux-full.qza")   , emit: reads
    tuple val(meta), path("${prefix}_demux-details.qza"), emit: details
    tuple val(meta), path("${prefix}_demux-summary.qzv"), emit: summary
    path "versions.yml", emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp
    export TMPDIR=./tmp
    export TMP=./tmp
    export TEMP=./tmp

    if [ "${single_end}" == "true" ]; then
        qiime demux emp-single \\
            --m-barcodes-file $barcodes \\
            --i-seqs $qza \\
            --o-per-sample-sequences ${prefix}_demux-full.qza \\
            --o-error-correction-details ${prefix}_demux-details.qza \\
            $args
    else
        qiime demux emp-paired \\
            --m-barcodes-file $barcodes \\
            --i-seqs $qza \\
            --o-per-sample-sequences ${prefix}_demux-full.qza \\
            --o-error-correction-details ${prefix}_demux-details.qza \\
            $args
    fi

    qiime demux summarize \\
        --i-data ${prefix}_demux-full.qza \\
        --o-visualization ${prefix}_demux-summary.qzv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(qiime info | sed -n '/: [0-9.]/p' | sed 's/.*/    &/g')
    END_VERSIONS
    """
}
