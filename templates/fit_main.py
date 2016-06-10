#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import pyvfit
from pyvfit.probability import Singlewle_Prob
import py2layer
from py2layer import TwoLayer_g_7K as TwoLayer_model

import os
import optparse
import sys
import numpy as np
import logging

LOG_FORMAT = '%(asctime)-15s %(name)-25s %(levelname)-8s %(message)s'
LOG_LEVEL = 'INFO'
LOG_FILENAME = 'local.log'
DIR_UVDATA = './uvdata/'

def setup(par):
    # read observations
    obs = pyvfit.ObsData(par['star']['starname'], par_star=par['star'], par_obs=par['obs'], dir_uvdata=DIR_UVDATA)

    # create model
    model = TwoLayer_model(obs.star, par=par['model'])
    model.compute_grids()

    if par['fit']['fit_mode'] == 'singlewle':
        obs.choose_wle(par['fit']['iwle_chosen'])
        model.set_wle_out(obs.wle*1.e3)  # obs.wle is [m]. and I should pass value in [mm].

    # create probability
    prob = Singlewle_Prob(np.array(par['fit']['par_range']), par['fit']['fixed_pars'])

    # create fitter
    fitter = pyvfit.Fitter(par['fit']['nwalkers'], par['fit']['nburnin'], par['fit']['ntotsteps'],
                           par['fit']['nstep_dump'], par['fit']['initial_step'],
                           par['fit']['p0'], np.array(par['fit'].get('initial_par_range', par['fit']['par_range'])))

    # create imagers
    imagers = []
    for i in xrange(prob.nwle):
        imagers.append(pyvfit.Imager(obs.uvtables[i].u, obs.uvtables[i].v, obs.uvtables[i].wle, obs.star.dist.cm,
                            obs.disk_ffcontr[i], maxuv_factor = 20.))

    return obs, model, prob, fitter, imagers


if __name__=='__main__':
    p=optparse.OptionParser()
    p.add_option("-f", "--file", action="store",dest="filename")
    p.add_option("-t", "--test", action="store_true",dest="is_test", help="Test if everything works.", default=False)

    (options,args) = p.parse_args()
    # Initialize the MPI-based pool used for parallelization.
    if options.is_test is False:
        # provare ad implementare il plan here
        from emcee.utils import MPIPool
        pool = MPIPool(loadbalance=True)
        if not pool.is_master():
            # Wait for instructions from the master process.
            pool.wait()
            sys.exit(0)

    log = logging.getLogger()
    log.setLevel(LOG_LEVEL)
    logfile = logging.FileHandler(filename=LOG_FILENAME, mode='w')
    logfile.setFormatter(logging.Formatter(LOG_FORMAT))
    log.addHandler(logfile)

    logging.info("")
    logging.info("****************************************************")
    logging.info("*****************      PyVFit      *****************")
    logging.info("****************************************************")

    # read fit parameters
    fit_filename = options.filename
    par = pyvfit.read_headers(fit_filename)  # in this way, it is possible to override some par values, e.g. the par_range

    # register software version
    par.update(software=dict(pyvfit=pyvfit.__version__, py2layer=py2layer.__version__))

    obs, model, prob, fitter, imagers = setup(par)

    if options.is_test == False:
        #pyvfit.run_emcee(obs, model, prob, fitter, imagers, pyvfit.lnprob, pool, par)
        pyvfit.run_emcee(obs, model, prob, fitter, imagers, pyvfit.lnprob_cuda, pool, par)
    else:
        print("Executing TEST")
        pyvfit.run_emcee(obs, model, prob, fitter, imagers, pyvfit.lnprob, None, par, is_test=True, test_p=[1., 1., 100., 10., 50., 0.01, 0.01])



