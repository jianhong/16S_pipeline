process SYNC_BARCODES {
    tag "$meta.id"
    tag "process_high"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(barcodes), path(reads1, stageAs: 'R1.fastq.gz'), path(reads2)

    output:
    tuple val(meta), path("${prefix}_R1.paired.fastq.gz"), path(reads2), path("${prefix}_I1.synced.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sync_paired_end_reads.py \\
        $barcodes $barcodes \\
        $reads1 \\
        ${prefix}_I1.synced.fastq.gz \\
        ${prefix}_R1.paired.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
