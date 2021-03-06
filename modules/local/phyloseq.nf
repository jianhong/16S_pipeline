process PHYLOSEQ {
    tag "$meta.id"
    tag "process_high"
    tag "error_ignore"

    conda (params.enable_conda ? "bioconda::bioconductor-phyloseq=1.38.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.38.0--r41hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-phyloseq:1.38.0--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(rds)
    path metadata

    output:
    tuple val(meta), path("$prefix")         , emit: phyloseq
    tuple val(meta), path("$krona_folder/*") , emit: krona, optional: true
    path 'phyloseq.out.txt'                  , emit: log
    path "versions.yml"                      , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    krona_folder = "4Krona"
    """
    #!/usr/bin/env Rscript
    # Applying dada2 pipeline to bioreactor time-series
    ## Following tutorial http://benjjneb.github.io/dada2_pipeline_MV/tutorial.html
    ## and here http://benjjneb.github.io/dada2_pipeline_MV/tutorial.html

    pkgs <- c("phyloseq")
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

    MAPPING <- "$metadata"
    NCORE <- as.numeric("$task.cpus")
    SEQTAB_S1 <- "seqtab.s1.rds"
    SEQTAB <- "seqtab.nochim.rds"
    TAXTAB <- "taxtab.rds"
    SAMPLENAMES <- "samplenames.1.rds"
    OUTFOLDER <- "$prefix"
    SAMPLEID_COL <- "SampleID"
    KRONA_FOLDER <- "$krona_folder"

    # Make phyloseq object ----------------------------------------------------
    # import data from dada2 output
    sample.names.1 <- readRDS(SAMPLENAMES)
    seqtab.s1 <- readRDS(SEQTAB_S1)
    seqtab <- readRDS(SEQTAB)
    taxtab <- readRDS(TAXTAB)

    # Import mapping
    map1 <- read.csv(MAPPING, stringsAsFactors = FALSE)
    map1 <- map1[map1[, SAMPLEID_COL] %in% sample.names.1,]
    map <- as.data.frame(map1) # without this line get sam_data slot empty error from phyloseq
    rownames(map) <- map[, SAMPLEID_COL]

    # Make refseq object and extract sequences from tables
    refseq <- colnames(seqtab)
    names(refseq) <- paste0('seq_', seq_along(refseq))
    colnames(seqtab) <- names(refseq[match(colnames(seqtab), refseq)])
    rownames(taxtab) <- names(refseq[match(rownames(taxtab), refseq)])

    # Write the taxtable, seqtable, and refseq to ascii ------------------------
    dir.create(OUTFOLDER, recursive=TRUE)
    write.table(seqtab, file=file.path(OUTFOLDER, 'seqtab.nochim.tsv'), quote=FALSE, sep='\\t')
    write.table(taxtab, file=file.path(OUTFOLDER, 'taxtab.nochim.tsv'), quote=FALSE, sep='\\t')
    write.table(refseq, file=file.path(OUTFOLDER, 'refseqs.nochim.tsv'), quote=FALSE, sep='\\t', col.names = FALSE)

    # Combine into phyloseq object
    ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(map), tax_table(taxtab))
    saveRDS(ps, file.path(OUTFOLDER, 'phyloseq.rds'))

    # copy the log file
    file.copy(".command.log", "phyloseq.out.txt")

    # output file for krona
    df <- psmelt(ps)
    cn <- c("Abundance", SAMPLEID_COL, rank_names(ps))
    cn <- intersect(cn, colnames(df))
    if(SAMPLEID_COL %in% cn){
        for(i in rank_names(ps)){
            if(!i %in% cn){
                df[, i] <- NA
            }
        }
        df <- df[, c("Abundance", SAMPLEID_COL, rank_names(ps)), drop=FALSE]
        df[, SAMPLEID_COL] <- as.factor(gsub("[^0-9a-zA-Z]", "_", df[, SAMPLEID_COL])) ## trim the sample names
        dir.create(KRONA_FOLDER, recursive=TRUE)
        for(lvl in levels(df[, SAMPLEID_COL])){
            write.table(
                df[which(df[, SAMPLEID_COL] == lvl & df[, "Abundance"] != 0), -2, drop=FALSE],
                file = file.path(KRONA_FOLDER, paste0(lvl, ".txt")),
                sep = "\\t", row.names = FALSE, col.names = FALSE,
                na = "", quote = FALSE)
        }
    }
    """
}
