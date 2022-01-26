process BCL2FASTQ {
    tag "$raw"
    tag 'process_high'

    input:
    path raw
    path samplesheet
    val bcl2fastq

    output:
    path 'fastq/Undetermined_*.fastq.gz' , emit: reads
    path "versions.yml"                  , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    mkdir -p fastq
    ${bcl2fastq} \\
            --runfolder-dir $raw \\
            --output-dir fastq/ \\
            --sample-sheet $samplesheet \\
            --processing-threads $task.cpus \\
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl2fastq: \$(echo \$(bcl2fastq --version 2>&1) | sed 's/^.*bcl2fastq v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
