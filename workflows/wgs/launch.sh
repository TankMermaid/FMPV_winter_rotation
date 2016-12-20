#! /bin/bash

WD=/netscr/apm/do_wgs
NODES=64

~/bin/snakemake --snakefile align.snake \
	--configfile config.yaml \
	--directory $WD --printshellcmds \
	--cluster /proj/pmdvlab/do_wgs/launch.py --jobs $NODES
