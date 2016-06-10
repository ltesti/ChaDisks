# @ shell=/bin/bash
#
# Template script for c2pap
# Marco Tazzari (ESO)
# 16 May 2016
#
# @ group =  pr84wo
# @ error = job_err_Z1_$(jobid)
# @ output = job_out_Z1_$(jobid)
# @ job_name = Z1_$(jobid)
# @ initialdir = /gpfs/work/pr84wo/ru49pep2/lupus/g/Z1
# @ job_type = parallel
# @ class = parallel
# @ node_usage = shared
# @ blocking = unlimited
# @ total_tasks = 55
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit = 48:00:00
# @ notification = always
# @ notify_user = mtazzari@eso.org
# @ restart = no
# @ queue

export MP_TASK_AFFINITY=cpu:1

# run the program
poe python fit_main.py -f fit_Z1.dat > output_Z1.out

