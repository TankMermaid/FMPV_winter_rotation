# WGS alignment with Snakemake

## Rationale
This workflow performs alignment using `bwa mem` with optical duplicates marked by `samblaster`.  Alignment happens in parallel over "read groups" (generally taken to be lane-by-sample combinations) and the results are merged and sorted into a single BAM file per sample.  Read groups are assigned to libraries and libraries to samples in the `run_map.txt` file.  It is assumed that each read group corresponds to a pair of FASTQ files, and that the FASTQ files for each sample are organized in a single directory whose name is the sample name.  If this is not the case, either move files or better yet create some symlinks to achieve the right (virtual) directory structure.

## Parameters
All adjustable parameters are in the file `config.yaml`.  They are listed below.  All file paths are relative to the working directory unless given as full paths.

* `reference`: reference genome for alignment, assumed to be already indexed for `bwa`
* `logs`: destination for LSF job reports
* `tmpdir`: temporary files directory; it's recommended not to rely on the system-wide `/tmp` as it has a tendency to fill up
* `readbuffer`: how many reads to keep in memory during sorting and merging of BAM files; tweak up or down to control balance between memory usage and number of temporary files created
* `threads`: default number of threads to use if not specified in individual tasks
* `memory`: default memory request for LSF (in GB) if not specificed in individual tasks
* `runs`: tab-delimited file mapping read groups to libraries and libraries to samples
* `aligned`: root directory for alignment output
* `fastq`: root directory for raw FASTQ files, assumed to have one subdirectory per sample

## Prerequisites

* `python` >= 3.3
* `snakemake` >= 3.8.1
*  `bwa`
*  `samblaster`
*  `samtools` >= 1.0
*  `picard`

Note that the version of Snakemake required is newer than the one installed by the Killdevil system admins.  The newer version allows the use of YAML config files instead of JSON, which I like.

## Invocation
Invoke the workflow like
```
snakemake --snakefile align.snake --configfile config.yaml --directory /working/dir --cluster launch.py
```

The script `launch.py` gets some job metadata from the Snakemake engine and uses it to assemble a call to the LSF scheduler (`bsub`).  When per-task resource requirements (processors, memory etc.) are not specified, defaults will be used that are in my experience appropriate for the Killdevil environment.

I used a wrapper script like `launch.sh` to start the master Snakemake job.