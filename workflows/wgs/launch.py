#!/usr/bin/env python3
import os
import sys

sys.path.append( os.path.expanduser("~/lib/python3") )

from snakemake.utils import read_job_properties

DEFAULT_PROPERTIES = {	"queue": ("day", "-q {queue}"), "threads": ("1", "-n {threads}"),
			"memory": (8, "-M {memory}"),
			"exclusive": ("", "-x"), "log_dir": ("jobs/job%J.out", "-o {log_dir}") }

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)["params"]
#os.system("echo " + str(read_job_properties(jobscript)))

# os.system("echo {} >>test.out".format(str(job_properties)))

def get_job_property(prop):
	
	if prop in job_properties:
		return (job_properties[prop], DEFAULT_PROPERTIES[prop][1])
	elif prop in DEFAULT_PROPERTIES:
		return DEFAULT_PROPERTIES[prop]
	else:
		return (None, None)

this_job = ""
for prop in DEFAULT_PROPERTIES.keys():
	# os.system("echo \"{}\" >>test.out".format(str(get_job_property(prop))))
	(val, fmt) = get_job_property(prop)
	if val is not None:
		# print(prop)
		# print(val)
		if val:
			this_job += " " + str(fmt).format(**{ prop: str(val) })

cmd = "bsub -R \"span[hosts=1]\" {opts} {script}".format(opts = this_job, script = jobscript)
os.system(cmd)
#os.system("echo \"{}\" >>test.out".format(cmd))
