process FILTERING {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("$prefix")     , emit: reads
    path '*.png'                         , emit: qc
    path "versions.yml"                  , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript

    # Applying dada2 pipeline to bioreactor time-series
    ## Following tutorial http://benjjneb.github.io/dada2_pipeline_MV/tutorial.html
    ## and here http://benjjneb.github.io/dada2_pipeline_MV/bigdata.html

    pkgs <- c("dada2", "ShortRead", "ggplot2")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    DEMUXD_READS = "$reads"
    TRIM = grepl("trimming_reads", "$args")
    OUTFOLDER = "$prefix"

    # Filtering and Trimming --------------------------------------------------
    #note s1 - is run number
    #note r1 - is forward, r2 - is reverse
    # Forward and Reverse Filenames
    files <- list.files(DEMUXD_READS)
    fnFs.s1 <- files[grepl("_R1_", files)]
    fnRs.s1 <- files[grepl("_R2_", files)]

    # Sort to ensure filenames are in the same order
    fnFs.s1 <- sort(fnFs.s1)
    fnRs.s1 <- sort(fnRs.s1)

    sample.names.1 <- sapply(strsplit(fnFs.s1,".fastq", fixed=TRUE), `[`, 1)

    # Fully Specify the path for the fnFs and fnRs
    fnFs.s1 <- file.path(DEMUXD_READS, fnFs.s1)
    fnRs.s1 <- file.path(DEMUXD_READS, fnRs.s1)

    # Examine qulaity profiles of the forward and reverse reads
    p <- plotQualityProfile(fnFs.s1[[1]])
    ggsave('Forward_quality_profile_s1.png', plot=p)
    p <- plotQualityProfile(fnRs.s1[[1]])
    ggsave('Reverse_quality_profile_s1.png', plot=p)

    # Perform filtering and trimming

    #update trimLeft based on plot output
    # For the first sequencing run
    dir.create(OUTFOLDER)
    filtFs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_F_filt.fastq.gz"))
    filtRs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_R_filt.fastq.gz"))
    if(TRIM){ # Reads look really good quality don't filter here
        trimLeft = c(10, 0)
        truncLen = c(150, 150)
    }else{
        trimLeft = c(0, 0)
        truncLen = c(0, 0)
    }
    for (i in seq_along(fnFs.s1)){
        fastqPairedFilter(c(fnFs.s1[i], fnRs.s1[i]), c(filtFs.s1[i], filtRs.s1[i]),
                          trimLeft=trimLeft, truncLen=truncLen,
                          maxN=0, maxEE=2, truncQ=2,
                          compress=TRUE, verbose=TRUE,
                          rm.phix=TRUE)
    }
    """
}
