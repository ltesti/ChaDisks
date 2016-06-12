# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# from __future__ import (division, print_function, absolute_import,
#                         unicode_literals)

import pyvfit
from pyvfit import wizard

import sys
import os
import json
import numpy as np
from subprocess import call
from astropy.io import ascii
import optparse

p = optparse.OptionParser()
p.add_option("-d", "--dry-run", action="store_true", dest="dry_run")
p.add_option("-a", "--analysis", action="store_true", dest="analysis_only")
(options,args) = p.parse_args()
dry_run = options.dry_run
analysis_only = options.analysis_only

WORKING_DIR = "./"
TEMPLATES_DIR = os.path.join(WORKING_DIR, "templates")
ROOT_FIT_DIR = os.path.join(WORKING_DIR, "")  # 'general_test/'
UVDATA_DATABASE_DIR = '../uvdata/'
PARAMETERS_TABLE = 'table_wizard.txt'
TOTAL_TASKS = '55'
NODE = '5'
HYDRA_ROOT_PREFIX = "/ptmp/mtazzari/lupus/g/"

# *** CREATE parameters input files (e.g. fit_xx.dat) *** #
wle = 0.890  # mm  Observing wavelength (to be used for both wle_out_mm and wle_mm)

# MODEL PARAMETERS
# radial grid
nrad = 500
rmin = 0.1  # AU
rmax = 500.  # AU

# dust
nsize = 200
amin = 5.e-7  # [cm]
amax = 10.  # [cm]
a0min = {'mid': 1.e-6,  # cm
         'sur': 1.e-6}
a0max = {'sur': 1.e-4}  # cm
q = {'mid': 3.0, 'sur': 3.5}

bmax = {'sur': 0.}
ncomp = 4
opt_const_filenames = {"Si": 'in_silicate.dat', "Ca": 'in_aC_ACH2_Zubko.dat', "WaterIce": 'in_H2Oice_Warren.dat'}
fractions = {"Si": 5.4e-2, "Ca": 20.6e-2, "WaterIce": 43.999e-2}

# wle grid
nwle = 200  # if > 400 the model crashes
wlemin = 1.e-5  # [cm]
wlemax = 1.  # [cm]

wle_out_mm = [wle]
sigma_prescription = 'g'

par_model = {'nrad': nrad, 'rmin': rmin, 'rmax': rmax, 'nsize': nsize, 'amin': amin, 'amax': amax, 'a0min': a0min,
             'a0max': a0max, 'q': q, 'bmax': bmax, 'ncomp': ncomp, 'opt_const_filenames': opt_const_filenames,
             'fractions': fractions, 'nwle': nwle, 'wlemin': wlemin, 'wlemax': wlemax, 'wle_out_mm': wle_out_mm,
             'sigma_prescription': sigma_prescription,}

# OBSERVATION PARAMETERS
# data_filename: read from input table
wle_mm = [wle]
disk_ffcontr = [0.]
weight_corr = [1.]

par_obs = {'wle_mm': wle_mm, 'disk_ffcontr': disk_ffcontr, 'weight_corr': weight_corr}

# FITS PARAMETERS
nwalkers = 80
nburnin = 200
nstep_dump = 100
initial_step = 0
ntotsteps = 20000
fit_mode = 'singlewle'
par_names = [r'$\gamma$', r'$\Sigma_{0}$', r'$R_{c}$', r'$i$', r'$PA$', r'$\Delta \alpha$ 0',
             r'$\Delta \delta$ 0']
iwle_chosen = 0

fixed_pars = {'a0max_fixed': 1.0213, 'bmax_fixed': 0.}
# a_amax = 1.0213cm  chosen to reproduce kappa_890mu=3.37cm2/g (like in Ansdell+2016)

# to be changed to [-2., 2.] for Exp-taper models.
range_gamma = [-2., 2.]  # gamma>2 induces crashes in the disk model (zbrent exceeds_
range_SS = np.log10(np.array([0.05, 300.]))
range_RR = [2., 300.]
range_inc = [0., 90.]
range_PA = [0., 180.]
range_x0 = [-2., 2.]
range_y0 = [-2., 2.]
par_range = [range_gamma, range_SS.tolist(), range_RR, range_inc, range_PA, range_x0, range_y0]

# initialization
# par_min = np.array([0.5, np.log10(0.3), 15., 1., 1., -0.1, -0.1])
# par_max = np.array([1.5, np.log10(180.), 180., 89., 179., 0.1, 0.1])

par_fit = {'nwalkers': nwalkers, 'nburnin': nburnin, 'nstep_dump': nstep_dump, 'initial_step': initial_step,
           'ntotsteps': ntotsteps, 'fit_mode': fit_mode, 'par_names': par_names, 'iwle_chosen': iwle_chosen,
           'fixed_pars': fixed_pars, 'par_range': par_range, 'p0': None, 'restarting': False,
           'sigma_prescription': sigma_prescription,  # can be discarded?
           'log_level': 'INFO',  # can be discarded?
}

# *** PREPARE FIT ENVIRONMENT *** #
with open(os.path.join(TEMPLATES_DIR, '_source.py'), 'r') as f:
    source_file_tmpl = f.read()

with open(os.path.join(TEMPLATES_DIR, 'hydra_script.sh'), 'r') as g:
    ll_script_tmpl = g.read()
fit_main_filename = os.path.join(TEMPLATES_DIR, 'fit_main.py')

# Import stellar parameters and other parameters from input table
# read table
ll_submit = []
table_parameters = ascii.read(PARAMETERS_TABLE, comment='#')

for source in table_parameters:

    # read source parameters
    starname = source['ID']
    source_label = source['Name']
    dist = np.float(source['Dist'])
    temp = np.float(source['Teff'])
    lum = np.float(source['Lstar'])
    mass = np.float(source['Mstar'])
    iflar = np.float(source['iflar'])

    print(starname, dist, temp, lum, mass, iflar, source_label)
    # if dry_run: continue

    # create directory
    fit_dir = os.path.join(ROOT_FIT_DIR, starname)
    if not dry_run:
        print("Creating directory {0}".format(fit_dir))
        call("mkdir -p {0}".format(fit_dir), shell=True)

    # STAR PARAMETERS
    par_star = {'starname': starname, 'lum': lum, 'temp': temp, 'mass': mass, 'dist': dist,}

    # MODEL PARAMETERS - UPDATE
    par_model.update(iflar=iflar)

    # OBSERVATION PARAMETERS - UPDATE
    obs_filename = "alma_b7_" + starname + ".txt"
    par_obs.update(data_filenames=[obs_filename])

    # create total parameters
    par = dict(model=par_model, star=par_star, obs=par_obs, fit=par_fit)

    input_filename = os.path.join(fit_dir, "fit_" + starname + '.dat')
    print("Creating {0}".format(input_filename))
    if not dry_run:
        with open(input_filename, 'w') as j:
            json.dump(par, j, sort_keys=True, indent=1)

    # fill scripts with proper values
    source_file = source_file_tmpl.replace("<starname>", starname)
    source_file = source_file.replace("<bestfit_pars>", 'None')
    source_file = source_file.replace("<source_label>", source_label.replace('<space>', ' '))
    analysis_filename = os.path.join(fit_dir, "_source.py")
    if not dry_run or analysis_only:
        print("Creating {0}".format(analysis_filename))
        with open(analysis_filename, 'w') as h:
            h.write(source_file)

    ll_script = ll_script_tmpl.replace("<starname>", starname)
    ll_script = ll_script.replace("<fit_dir>", os.path.join(HYDRA_ROOT_PREFIX, starname))
    ll_script = ll_script.replace("<total_tasks>", TOTAL_TASKS)
    ll_script = ll_script.replace("<node>", NODE)
    # ll_script = ll_script.replace("<tasks_per_node>", '16')
    # ll_script = ll_script.replace("<node>", '3')
    # ll_script = ll_script.replace("<first_node_tasks>", '11')
    script_filename = os.path.join(fit_dir, "script_" + starname + ".sh")
    if not dry_run:
        print("Creating {0}".format(script_filename))
        with open(script_filename, 'w') as h:
            h.write(ll_script)

    fit_main_dest = os.path.join(fit_dir, "fit_main.py")
    if not dry_run:
        print("Creating {0}".format(fit_main_dest))
        call("cp {0} {1}".format(fit_main_filename, fit_main_dest), shell=True)

    uvdata_orig = os.path.join(UVDATA_DATABASE_DIR, obs_filename)
    uvdata_dest = os.path.join(os.path.join(fit_dir, "uvdata"), obs_filename)
    if not dry_run:
        print("Copying {0} {1}".format(uvdata_orig, uvdata_dest))
        call("mkdir -p {0}".format(os.path.join(fit_dir, "uvdata")), shell=True)
        call("cp {0} {1}".format(uvdata_orig, uvdata_dest), shell=True)

    ll_submit.append("llsubmit {0}".format(os.path.join(fit_dir, "script_" + starname + ".sh")))

if not dry_run:
    script_filename = os.path.join(WORKING_DIR, "ll_submit.sh")
    print("Creating {0}".format(script_filename))
    with open(script_filename, 'w') as h:
    	h.write("\n".join(ll_submit))
