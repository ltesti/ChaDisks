# @ shell=/bin/bash
#
# Template script for Hydra
# Marco Tazzari (ESO)
# 24 Mar 2016
#
# @ error = job_err_<starname>_$(jobid)
# @ output = job_out_<starname>_$(jobid)
# @ job_name = <starname>_$(jobid)
# @ initialdir = <fit_dir>
# @ job_type = parallel
# @ node_usage = shared
# @ node = <node>  
# @ total_tasks = <total_tasks>
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit = 24:00:00
# @ notification = always
# @ notify_user = mtazzari@rzg.mpg.de
# @ restart = no
# @ queue

export MP_TASK_AFFINITY=cpu:1

# run the program
poe python fit_main.py -f fit_<starname>.dat > output_<starname>.out
