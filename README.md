# [![Cite with Zenodo](https://zenodo.org/badge/451935888.svg)](https://zenodo.org/badge/latestdoi/451935888)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**16S_pipeline** is a bioinformatics best-practice analysis pipeline for A pipeline for 16S rRNA gene sequencing and shotgun metagenomic sequencing.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

0. Prepare fastq files ([`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html))
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Remove primers ([`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic))
3. Sync barcodes ([`fastq_pair_filter.py`](https://gist.github.com/588841/))
4. Demultiplex ([`qiime2::demux`](https://docs.qiime2.org/2021.11/plugins/available/demux/))
5. Filter reads ([`DATA2`](http://benjjneb.github.io/dada2/))
6. Run dada2 ([`DATA2`](http://benjjneb.github.io/dada2/))
7. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run jianhong/16S_pipeline -profile test,YOURPROFILE
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run jianhong/16S_pipeline -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '[path to raw reads files]' --barcodes '[path to barcodes tsv file]' --metadata '[path to metadata csv file]'
    ```

    Run it on cluster.

    First prepare a profile config file named as [profile.config](docs/usage.md).

    ```console
    // submit by slurm
    process.executor = "slurm"
    process.clusterOptions = "-J microbiome"
    params {
        // Input data
        input  = 'path/to/your/fastqfiles'
        skip_bcl2fastq = true
        barcodes = 'path/to/your/barcodes.tsv'
        metadata = 'path/to/your/metadata.csv'

        // Genome references
        silva_nr99 = 'https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1'
        silva_tax = 'https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1'
    }
    ```

    Then run:

    ```console
    nextflow run jianhong/16S_pipeline -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> -c profile.config
    ```

## Create `conda` container for `nextflow`

1. Install [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

    ```console
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

2. Create [`nextflow`](https://www.nextflow.io/) environment.

    ```console
    conda create -y --name microbiome bioconda::nextflow=21.10.6
    ```

3. Create profile config file named as [profile.config](docs/usage.md).

    ```console
    // submit by slurm
    process.executor = "slurm"
    process.clusterOptions = "-J JL21"
    params {
        // Input data
        input  = 'Spinach_done' // replace Spinach_done by your own file
        barcodes = '16S_pipeline_JL21/0_mapping/barcodes.tsv'
        metadata = '16S_pipeline_JL21/0_mapping/metadata.csv'

        // report email
        email = 'your@email.addr'
    }
    ```

4. Activate [`nextflow`](https://www.nextflow.io/) environment and Run the pipeline.

    ```console
    conda activate microbiome
    module load bcl2fastq/2.20
    nextflow run jianhong/16S_pipeline -r main -profile conda -c profile.config
    ```

## Documentation

The 16S_pipeline pipeline comes with documentation about the pipeline [usage](docs/usage.md), and [output](docs/output.md).

## Credits

16S_pipeline was originally written by Jianhong Ou, and Jeff Letourneau.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use  16S_pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.451935888](https://doi.org/10.5281/zenodo.451935888)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
