# 16S_pipeline: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

16S_pipeline is a bioinformatics analysis pipelline used for 16S microbiome analysis. Supported is paired-end illumina data.
The pipeline is based on [DADA2](http://benjjneb.github.io/dada2/).

## Raw reads, barcodes and metadata

You will need to create a barcodes and metadata tables with information about the samples you would like to analyze before running the pipeline. Use this parameter to specify its location.

```console
--input '[path to raw reads files]' --barcodes '[path to barcodes tsv file]' --metadata '[path to metadata csv file]'
```

### Raw reads

There are two choices for the raw reads input.

1. The fastq.gz files exported by [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html). The fastq.gz with R1, R2 and I1 files should be put into one folder and the folder path will be the parameter of `input`. And `--skip_bcl2fastq` should be set.
2. The Illumina intensity files ready to be handled by `bcl2fastq`. The folder will be the parameter of `input`. Please make sure the `bcl2fastq` is installed. For `module` package available system (check by command `module avail`), the `bcl2fastq` could be load into the `PATH` by following sample code:

```console
module load bcl2fastq/2.20
```

### barcodes

It has to be a tab-separated file with 2 columns, and a header row as shown in the examples below.
It will be used to do demultiplex by [qiime2::demux](https://docs.qiime2.org/2021.11/plugins/available/demux/).

```console
sample-id barcode-sequence
#q2:types categorical
Spinach1  TGTGCGATAACA
Spinach2  GATTATCGACGA
Spinach3  GCCTAGCCCAAT
....
```

An [example samplesheet](../assets/barcodes.tsv) has been provided with the pipeline.

### Metadata

The `metadata` contain the information about the samples you wold like to analyze.
It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```console
SampleID,Character1,Character2,Character3
SAMPLE1,A,treatment,1
SAMPLE2,B,treatment,1
SAMPLE3,C,control,1
```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `SampleID`            | Custom sample name. This entry will be unique for each sample and should not contain any special characters. |
| `CharacterX`      | metadata for each sample. The column name can be anything related with the samples.  |

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run jianhong/16S_pipeline --input raw --barcodes barcodes.tsv --metadata metadata.csv -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull jianhong/16S_pipeline
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [16S_pipeline releases page](https://github.com/16S_pipeline/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customize the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the pipeline is failing after multiple re-submissions of the `BCL2FASTQ` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[9d/172ca5] NOTE: Process `NFCORE_MICROBIOME:MICROBIOME:BCL2FASTQ` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_MICROBIOME:MICROBIOME:BCL2FASTQ'

Caused by:
    Process `NFCORE_MICROBIOME:MICROBIOME:BCL2FASTQ` terminated with an error exit status (137)

Command executed:
    bcl2fastq \
            --runfolder-dir raw \
            --output-dir fastq/ \
            --sample-sheet samplesheet.csv \
            --processing-threads 2

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    bcl2fastq --runfolder-dir raw --output-dir fastq/ <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `BCL2FASTQ` process. The quickest way is to search for `process BCL2FASTQ` in the [jianhong/16S_pipeline Github repo](https://github.com/jianhong/16S_pipeline/search?q=process+BCL2FASTQ). We have standardized the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/local/bcl2fastq.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to `label process_high`. The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/jianhong/16S_pipeline/blob/c7f5684d49151ac7974e7eabc7915f5e4a0fd3aa/conf/base.config#L39-L43) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `BCL2FASTQ` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: BCL2FASTQ {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Module specific parameters

Module specific parameters can be passed by config file in following format by setting the `ext.args`.

```nextflow
process {
    withName: FILTERING {
        ext.args   = '--trimming_reads --trim_left 10 --trim_right 0 --trunc_length_left 150 trunc_length_right 150'
        ext.prefix = 's1'
        publishDir = [
            path: { "${params.outdir}/4_filter" },
            mode: 'copy',
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

For `FILTERING`, the options are

```{R}
--trimming_reads, -t        "logical",  Trim reads or not.
--trim_left, -a,            "integer",  Default 0. The number of nucleotides to remove
                                        from the start of the R1 read. If both
                                        trunc_length_left and trim_left are provided,
                                        filtered reads will have length
                                        trunc_length_left-trim_left.
--trim_right, -b,           "integer",  Default 0. The number of nucleotides to remove
                                        from the start of the R2 reads. If both
                                        trunc_length_right and trim_right are provided,
                                        filtered reads will have length
                                        trunc_length_right-trim_right.
--trunc_length_left, -m,    "integer",  Default 0 (no truncation). Truncate R1 reads
                                        after trunc_length_left bases. Reads shorter
                                        than this are discarded.
--trunc_length_right, -n,   "integer",  Default 0 (no truncation). Truncate R2 reads
                                        after trunc_length_right bases. Reads shorter
                                        than this are discarded.
```

For `DADA2`, the options are

```{R}
--seq1, -a,           "integer",  The number of minimal sequence length should be kept.
--seq2, -b,           "integer",  The number of maximum sequence length should be kept.
```

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
