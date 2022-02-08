process REMOVE_PRIMERS {
    tag "$meta.id"
    tag 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2' :
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    val single_end

    output:
    tuple val(meta), path("${prefix}_R1.paired.fastq.gz"), path("${prefix}_R2.paired.fastq.gz"), emit: paired
    tuple val(meta), path("${prefix}_R1.unpaired.fastq.gz"), path("${prefix}_R2.unpaired.fastq.gz"), optional:true, emit: unpaired
    path '*.trim_out.log', emit: log
    path "versions.yml", emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "${single_end}" == "true" ]; then
        touch ${prefix}_R2.paired.fastq.gz
        trimmomatic SE \\
                    -threads $task.cpus \\
                    ${reads[0]} \\
                    ${prefix}_R1.paired.fastq.gz \\
                    $args \\
                    > ${prefix}.trim_out.log 2>&1
    else
        trimmomatic PE \\
                    -threads $task.cpus \\
                    ${reads[0]} ${reads[1]} \\
                    ${prefix}_R1.paired.fastq.gz \\
                    ${prefix}_R1.unpaired.fastq.gz \\
                    ${prefix}_R2.paired.fastq.gz \\
                    ${prefix}_R2.unpaired.fastq.gz \\
                    $args \\
                    > ${prefix}.trim_out.log 2>&1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
