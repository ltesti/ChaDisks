# @ shell=/bin/bash
#
# Template script for Lupus
# Marco Tazzari (ESO)
#
# @ error = job_err_<starname>_$(jobid)
# @ output = job_out_<starname>_$(jobid)
# @ job_name = <starname>_$(jobid)
# @ initialdir = <fit_dir>
# @ job_type = parallel
# @ node_usage = shared
# @ node = <node>
# @ tasks_per_node = <tasks_per_node>
#### @ first_node_tasks = 6
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit = 48:00:00
# @ notification = always
# @ notify_user = mtazzari@rzg.mpg.de
# @ restart = no
# @ queue

export MP_TASK_AFFINITY=cpu:1

# run the program
poe python fit_main.py -f fit_<starname>.dat > output_<starname>.out


# Notes for Hydra
# @ blocking = unlimited	# cannot be used on Hydra
# @ node = 5  				# Hydra does not allow to use <11 CPUs/node
on Hydra node must *always* be specified according to:
	node && (total_tasks || tasks_per_node || tasks_per_host) || task_geometry