process DADA2 {
    tag "$meta.id"
    tag "process_high"

    conda (params.enable_conda ? "bioconda::bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads), path(stats)
    path train_set
    path species_assignment

    output:
    tuple val(meta), path("*.rds")       , emit: robj
    path '*.{png,csv}'                   , emit: qc
    path 'dada2.out.txt'                 , emit: log
    path "versions.yml"                  , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    # Applying dada2 pipeline to bioreactor time-series
    ## Following tutorial http://benjjneb.github.io/dada2_pipeline_MV/tutorial.html
    ## and here http://benjjneb.github.io/dada2_pipeline_MV/tutorial.html

    pkgs <- c("dada2", "ggplot2")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    set.seed(4)

    FILTPATH <- "$reads"
    NCORE <- ifelse(.Platform[["OS.type"]]!="windows", as.numeric("$task.cpus"), FALSE)
    TRAIN_SET <- "$train_set"
    SPECIES_ASSIGNMENT <- "$species_assignment"
    STATS <- "$stats"
    SEQL1 <- 0
    SEQL2 <- 0
    SAMPLENAMES <- "samplenames.1.rds"
    SEQTAB_S1 <- "seqtab.s1.rds"
    SEQTAB <- "seqtab.nochim.rds"
    TAXTAB <- "taxtab.rds"

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
    option_list <- list("seqlen1"=c("--seq1", "-a", "integer"),
                        "seqlen2"=c("--seq2", "-b", "integer"))
    opt <- parse_args(option_list, args)
    if(!is.null(opt[["seqlen1"]])){
        SEQL1 <- opt[["seqlen1"]]
    }
    if(!is.null(opt[["seqlen2"]])){
        SEQL2 <- opt[["seqlen2"]]
    }

    # stats
    getN <- function(x) sum(getUniques(x))
    track <- readRDS(STATS)

    # Find filenames ----------------------------------------------------------

    # Forward and reverse filenames
    filts.s1 <- list.files(FILTPATH, full.names=TRUE)

    # Sort to ensure fileneames are in the same order
    filts.s1 <- sort(filts.s1)
    sample.names.1 <- sapply(strsplit(basename(filts.s1),"_"), `[`, 1)
    names(filts.s1) <- sample.names.1


    # Separate forward and reverse samples
    filtFs.s1 <- filts.s1[grepl("_F_filt",filts.s1)]
    filtRs.s1 <- filts.s1[grepl("_R_filt",filts.s1)]

    sample.names.1 <- sapply(strsplit(basename(filtFs.s1), "_"), `[`, 1)
    saveRDS(sample.names.1, SAMPLENAMES)

    # Dereplication -----------------------------------------------------------

    # Learn Error Rates
    ## aim to learn from about 1M total reads - so just need subset of samples
    ## source: http://benjjneb.github.io/dada2_pipeline_MV/bigdata.html
    filts.learn.s1 <- sample(sample.names.1, 36)

    derepFs.s1.learn <- derepFastq(filtFs.s1[filts.learn.s1], verbose=TRUE)
    derepRs.s1.learn <- derepFastq(filtRs.s1[filts.learn.s1], verbose=TRUE)

    # Sample Inference --------------------------------------------------------

    dadaFs.s1.learn <- dada(derepFs.s1.learn, err=NULL, selfConsist=TRUE, multithread=NCORE)
    dadaRs.s1.learn <- dada(derepRs.s1.learn, err=NULL, selfConsist=TRUE, multithread=NCORE)
    rm(derepFs.s1.learn, derepRs.s1.learn)

    # # Visualize estimated error rates
    p<- plotErrors(dadaFs.s1.learn[[1]], nominalQ=TRUE)
    ggsave("dada_errors_F_s1.png", plot=p)
    p<- plotErrors(dadaRs.s1.learn[[1]], nominalQ=TRUE)
    ggsave("dada_errors_R_s1.png", plot=p)


    # Just keep the error profiles
    errFs.s1 <- dadaFs.s1.learn[[1]][["err_out"]]
    errRs.s1 <- dadaRs.s1.learn[[1]][["err_out"]]
    rm(dadaFs.s1.learn, dadaRs.s1.learn)

    # Now sample inference for entire dataset
    # Run 1
    derepFs.s1 <- vector("list", length(sample.names.1))
    derepRs.s1 <- vector("list", length(sample.names.1))
    dadaFs.s1 <- vector("list", length(sample.names.1))
    dadaRs.s1 <- vector("list", length(sample.names.1))
    names(dadaFs.s1) <- sample.names.1
    names(dadaRs.s1) <- sample.names.1
    names(derepFs.s1) <- sample.names.1
    names(derepRs.s1) <- sample.names.1
    for (sam in sample.names.1){
        message("Processing:", sam, "\n")
        derepFs.s1[[sam]] <- derepFastq(filtFs.s1[[sam]])
        derepRs.s1[[sam]] <- derepFastq(filtRs.s1[[sam]])
        dadaFs.s1[[sam]] <- dada(derepFs.s1[[sam]], err=errFs.s1, multithread=NCORE)
        dadaRs.s1[[sam]] <- dada(derepRs.s1[[sam]], err=errRs.s1, multithread=NCORE)
    }

    # Run 1: Merge Paired Reads
    mergers.s1 <- mergePairs(dadaFs.s1, derepFs.s1, dadaRs.s1, derepRs.s1, verbose=TRUE)
    head(mergers.s1)
    track <- cbind(track, sapply(dadaFs.s1, getN), sapply(dadaRs.s1, getN),
                    sapply(mergers.s1, getN))
    # Run 1: Clear up space
    rm(derepFs.s1, derepRs.s1, dadaFs.s1, dadaRs.s1)


    # Construct Sequence Table ------------------------------------------------

    #To use only the forward reads
    #follow https://github.com/benjjneb/dada2/issues/134

    seqtab.s1 <- makeSequenceTable(mergers.s1)
    saveRDS(seqtab.s1, SEQTAB_S1)
    dim(seqtab.s1)
    # Inspect the distributioh of sequence lengths
    seqlenTab <- table(nchar(colnames(seqtab.s1)))
    write.csv(seqlenTab, "seqlenTab.csv", row.names=FALSE)
    if(SEQL1==0 && SEQL2==0){## auto detect cutoff range
        maxV <- which.max(seqlenTab)
        maxV <- as.numeric(names(seqlenTab)[maxV])
        if(length(maxV)>1){
            if(any(abs(diff(maxV))>2)){
                stop("can not determine cutoff SEQ1 and SEQ2")
            }
            maxV <- maxV[median(seq_along(maxV))]
        }
        SEQL1 <- maxV - 2
        SEQL2 <- maxV + 2
    }
    #Trim sequences of interest
    seqtab.s1 <- seqtab.s1[,nchar(colnames(seqtab.s1)) %in% seq(SEQL1,SEQL2)]
    # Inspect the distributioh of sequence lengths
    seqlenTab <- table(nchar(colnames(seqtab.s1)))
    write.csv(seqlenTab, "seqlenTab.filt.csv", row.names=FALSE)

    # Remove Chimeras ---------------------------------------------------------

    seqtab.s1.nochim <- removeBimeraDenovo(seqtab.s1, method='consensus', multithread=NCORE, verbose=TRUE)
    freq_chimeric <- c(number_row=nrow(seqtab.s1.nochim), number_col=ncol(seqtab.s1.nochim), frequency=sum(seqtab.s1.nochim)/sum(seqtab.s1))
    write.csv(t(freq_chimeric), "freq_chimeric.csv", row.names=FALSE)
    saveRDS(seqtab.s1.nochim, "seqtab.s1.nochim.rds")

    track <- cbind(track, rowSums(seqtab.s1.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names.1
    write.csv(track, "processing_tracking.csv")

    # Merge Sequence Tables Together ------------------------------------------

    seqtab.nochim <- seqtab.s1.nochim
    saveRDS(seqtab.nochim, SEQTAB)


    # Simplify naming ---------------------------------------------------------

    seqtab <- seqtab.nochim

    # Assign Taxonomy ---------------------------------------------------------
    # Following: http://benjjneb.github.io/dada2_pipeline_MV/species.html

    # Assign using Naive Bayes RDP
    taxtab <- assignTaxonomy(colnames(seqtab), TRAIN_SET, multithread=NCORE)

    # improve with exact genus-species matches
    # this step is pretty slow, should improve in later releases
    # - note: Not allowing multiple species matches in default setting
    taxtab <- addSpecies(taxtab, SPECIES_ASSIGNMENT, verbose=TRUE)
    saveRDS(taxtab, TAXTAB)

    # How many sequences are classified at different levels? (percent)
    classify_levels <- colSums(!is.na(taxtab))/nrow(taxtab)
    write.csv(t(classify_levels), "classify_levels.csv", row.names=FALSE)

    # copy the log file
    file.copy(".command.log", "dada2.out.txt")
    """
}
