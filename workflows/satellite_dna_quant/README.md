# Satellite DNA quantification with kallisto

## Rationale
This workflow quantifies a set of small targets (ie. centromeric satellite sequences) from whole-genome sequencing data using `kallisto`.  We assume that the starting point for each sample is a bam file containing both aligned and unaligned reads.  First the unaligned reads and their mates are stripped from the bam file, sorted and reverted to fastq format.  Then tarets are quantified from the (paired-end) fastq files using `kallisto`, from an index which is assumed to be built before the workflow is started.

## Parameters
All adjustable parameters are in the file `config.yaml`.  They are listed below.  All file paths are relative to the working directory unless given as full paths.

* `index`: an index of quantification targets built with a command like `kallisto index -i {index_path} {fasta file of targets}`.
* `bams`: text file listing bam files on which to perform quantification, one file per sample
* `threads`: default number of threads to use if not specified in individual tasks
* `memory`: default memory request for LSF (in GB) if not specificed in individual tasks
* `out`: root directory for `kallisto` output; will have one subdirectory for each sample
* `working`: directory for keeping intermediate fastq files (recommended: use `/netscr` or `/lustre` space)

## Prerequisites

* `python` >= 3.3
* `snakemake` >= 3.8.1
* `bwa`
* `samblaster`
* `samtools` >= 1.0
* `picard`
* `kallisto` 

Note that the version of Snakemake required is newer than the one installed by the Killdevil system admins.  The newer version allows the use of YAML config files instead of JSON, which I like.

## Invocation
Invoke the workflow like
```
snakemake --snakefile quant_sats_from_bam.snake --configfile config.yaml --directory /working/dir --cluster launch.py
```

The script `launch.py` gets some job metadata from the Snakemake engine and uses it to assemble a call to the LSF scheduler (`bsub`).  When per-task resource requirements (processors, memory etc.) are not specified, defaults will be used that are in my experience appropriate for the Killdevil environment.  See workflow in `../wgs`, in this same repository, for an example `launch.py` script.