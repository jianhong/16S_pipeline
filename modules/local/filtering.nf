process FILTERING {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("$prefix"), path("*.rds") , emit: reads
    path '*.png'                                    , emit: qc
    path "versions.yml"                             , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
    FILTER_STATS = "filter_stats.rds"
    NCORE <- ifelse(.Platform[["OS.type"]]!="windows", as.numeric("$task.cpus"), FALSE)

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
    p_F <- plotQualityProfile(fnFs.s1, aggregate=TRUE)
    p_R <- plotQualityProfile(fnRs.s1, aggregate=TRUE)
    ggsave('Forward_quality_profile_aggregate.png', plot=p_F)
    ggsave('Reverse_quality_profile_aggregate.png', plot=p_R)
    if(TRIM){ # Reads look really good quality don't filter here
        getTrimRange <- function(x){
            l <- lapply(x[["layers"]], function(.ele) .ele[["data"]])
            m <- lapply(x[["layers"]], function(.ele) .ele[["mapping"]])
            id <- lapply(m, function(.ele) any(grepl("Mean", as.character(.ele[["y"]]))))
            id <- which(unlist(id))
            if(length(id)>0){
                d <- l[[id[1]]]
                pos <- which(d[, "Mean"] < 30)
                if(length(pos)>0){
                    pos <- which(d[, "Mean"] >= 30)
                    ## start pos
                    pos_l <- pos[1]
                    ## end pos
                    pos_r <- pos[length(pos)]
                    c(pos_l, pos_r-pos_l+1)
                }else{
                    c(0, 0)
                }
            }else{
                c(0, 0)
            }
        }
        trim_range_L <- getTrimRange(p_F)
        trim_range_R <- getTrimRange(p_R)
        trimLeft = c(trim_range_L[1], trim_range_R[1])
        truncLen = c(trim_range_L[2], trim_range_R[2])
    }else{
        trimLeft = c(0, 0)
        truncLen = c(0, 0)
    }

    # Perform filtering and trimming

    #update trimLeft based on plot output
    # For the first sequencing run
    dir.create(OUTFOLDER)
    filtFs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_F_filt.fastq.gz"))
    filtRs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_R_filt.fastq.gz"))

    out <- filterAndTrim(fwd=fnFs.s1, filt=filtFs.s1,
                        rev=fnRs.s1, filt.rev=filtRs.s1,
                        trimLeft=trimLeft, truncLen=truncLen,
                        maxN=0, maxEE=2, truncQ=2,
                        compress=TRUE, verbose=TRUE,
                        rm.phix=TRUE, multithread=NCORE)
    saveRDS(out, FILTER_STATS)
    """
}
