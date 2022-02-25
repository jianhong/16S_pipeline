process FILTERING {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)
    val single_end

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
    TRIM = FALSE
    OUTFOLDER = "$prefix"
    FILTER_STATS = "filter_stats.rds"
    NCORE <- ifelse(.Platform[["OS.type"]]!="windows", as.numeric("$task.cpus"), FALSE)
    PAIRED_END <- "$single_end" != "true"

    args <- strsplit("${args}", "\\\\s+")[[1]]
    parse_args <- function(options, args){
        out <- lapply(options, function(.ele){
            if(any(.ele[-3] %in% args)){
                if(.ele[3]=="logical"){
                    TRUE
                }else{
                    id <- which(args %in% .ele[-3])[1]
                    x <- args[id+1]
                    mode(x) <- .ele[3]
                    x
                }
            }
        })
    }
    grepl("trimming_reads", "$args")
    option_list <- list("trim"=c("--trimming_reads", "-t", "logical"),
                        "trimLeft"=c("--trim_left", "-a", "integer"),
                        "trimRight"=c("--trim_right", "-b", "integer"),
                        "truncLenLeft"=c("--trunc_length_left", "-m", "integer"),
                        "truncLenRight"=c("--trunc_length_right", "-n", "integer"))
    opt <- parse_args(option_list, args)
    if(!is.null(opt[["trim"]])){
        TRIM <- TRUE
    }
    if(!is.null(opt[["trimLeft"]])){
        trimLeft <- trimLeft1 <- opt[["trimLeft"]]
    }else{
        trimLeft1 <- 0
        trimLeft <- NULL
    }
    if(!is.null(opt[["trimRight"]])){
        trimLeft <- c(trimLeft1, opt[["trimRight"]])
    }
    if(!is.null(opt[["truncLenLeft"]])){
        truncLen <- truncLen1 <- opt[["truncLenLeft"]]
    }else{
        truncLen1 <- 0
        truncLen <- NULL
    }
    if(!is.null(opt[["truncLenRight"]])){
        truncLen <- c(truncLen1, opt[["truncLenRight"]])
    }

    # Filtering and Trimming --------------------------------------------------
    #note s1 - is run number
    #note r1 - is forward, r2 - is reverse
    # Forward and Reverse Filenames
    if(dir.exists(DEMUXD_READS)){
        files <- list.files(DEMUXD_READS)
    }else{
        files <- strsplit(DEMUXD_READS, "\\\\s+")[[1]]
    }

    fnFs.s1 <- files[grepl("_R1[_.]", files)]
    fnRs.s1 <- files[grepl("_R2[_.]", files)]

    # Sort to ensure filenames are in the same order
    fnFs.s1 <- sort(fnFs.s1)
    fnRs.s1 <- sort(fnRs.s1)

    sample.names.1 <- sapply(strsplit(fnFs.s1, "_R1[_.].*?(fastq|fq)", fixed=FALSE), `[`, 1)
    if(PAIRED_END){
        ## match pairs
        sample.names.2 <- sapply(strsplit(fnRs.s1,"_R2[_.].*?(fastq|fq)", fixed=FALSE), `[`, 1)
        sample.names.shared <- intersect(sample.names.1, sample.names.2)
        fnFs.s1 <- fnFs.s1[match(sample.names.shared, sample.names.1)]
        fnRs.s1 <- fnRs.s1[match(sample.names.shared, sample.names.2)]
    }

    # Fully Specify the path for the fnFs and fnRs
    if(dir.exists(DEMUXD_READS)){
        fnFs.s1 <- file.path(DEMUXD_READS, fnFs.s1)
        fnRs.s1 <- file.path(DEMUXD_READS, fnRs.s1)
    }

    # Examine qulaity profiles of the forward and reverse reads
    p <- plotQualityProfile(fnFs.s1[[1]])
    ggsave('Forward_quality_profile_s1.png', plot=p)
    p_F <- plotQualityProfile(fnFs.s1, aggregate=TRUE)
    ggsave('Forward_quality_profile_aggregate.png', plot=p_F)
    if(PAIRED_END){
        p <- plotQualityProfile(fnRs.s1[[1]])
        ggsave('Reverse_quality_profile_s1.png', plot=p)
        p_R <- plotQualityProfile(fnRs.s1, aggregate=TRUE)
        ggsave('Reverse_quality_profile_aggregate.png', plot=p_R)
    }
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
        if(is.null(trimLeft[1]) || is.null(truncLen[1])){
            trim_range_L <- getTrimRange(p_F)
            if(PAIRED_END){
                trim_range_R <- getTrimRange(p_R)
                trimLeft = c(trim_range_L[1], trim_range_R[1])
                truncLen = c(trim_range_L[2], trim_range_R[2])
            }else{
                trimLeft = c(trim_range_L[1])
                truncLen = c(trim_range_L[2])
            }
        }
    }else{
        if(PAIRED_END){
            trimLeft = c(0, 0)
            truncLen = c(0, 0)
        }else{
            trimLeft = c(0)
            truncLen = c(0)
        }
    }

    # Perform filtering and trimming

    #update trimLeft based on plot output
    # For the first sequencing run
    dir.create(OUTFOLDER)
    filtFs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_F_filt.fastq.gz"))
    filtRs.s1 <- file.path(OUTFOLDER, paste0(sample.names.1,"_R_filt.fastq.gz"))
    if(PAIRED_END){
        out <- filterAndTrim(fwd=fnFs.s1, filt=filtFs.s1,
                            rev=fnRs.s1, filt.rev=filtRs.s1,
                            trimLeft=trimLeft, truncLen=truncLen,
                            maxN=0, maxEE=2, truncQ=2,
                            compress=TRUE, verbose=TRUE,
                            rm.phix=TRUE, multithread=NCORE)
    }else{
        out <- filterAndTrim(fwd=fnFs.s1, filt=filtFs.s1,
                            trimLeft=trimLeft, truncLen=truncLen,
                            maxN=0, maxEE=2, truncQ=2,
                            compress=TRUE, verbose=TRUE,
                            rm.phix=TRUE, multithread=NCORE)
    }
    saveRDS(out, FILTER_STATS)
    """
}
