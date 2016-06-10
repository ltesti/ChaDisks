#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from pyvfit import analysis, get_last_snapshot
from pyvfit.utils import util_plot, util, compute_percentiles, clean_chain, \
    get_chain

import numpy as np
from subprocess import call
import optparse
import logging
import textwrap
import time
import glob
import os

# uncomment to enable logging
# logger = logging.getLogger("analysis_main")
# console = logging.StreamHandler()
# console.setLevel(logging.INFO)
# logger.addHandler(console)

tasks = ["Plot MCMC",
         "Print Fit Results",
         "Compute Bestfit model",
         "Print Bestfit model structure",
         "Do Bestfit model uv plot",
         "Compute set of models",
         "Compute visibilities for set of models",
         "Do uv plot for set of models",
         "Compute physical quantities",
         "Prepare scripts for maps"]

p = optparse.OptionParser()
p.add_option("-s", action="store", dest="flags", default='0',
             help="Select commands to execute:\n"+
                  "\n".join(["{0}. {1}".format(i, task) for i, task in enumerate(tasks)]))
p.add_option("-d", "--diag", action="store_true", dest="chain_diagnostics",
             default=False, help="Produce chain diagnostics plots.")
p.add_option("--cpu", action="store", dest="ncpus", default=4,
             help="Number of CPUs to use for models and visibilities computation")
# p.add_option("--dry", "--dry-run", action="store_true", dest="dry_run",
#              help="Execute a dry run.")

(options, args) = p.parse_args()
ncpus = options.ncpus
chain_diagnostics = options.chain_diagnostics
flags_str = options.flags
# dry_run = options.dry_run

do_2d_contours = False
do_plot_surface_density = False
do_plot_temp_mid = False


def logtoreal(flatchain, bestfit_pars, logpars):
    flatchain_real = flatchain.copy()
    bestfit_pars_real = bestfit_pars.copy()

    for logpar in logpars:
        flatchain_real[:, logpar] = 10. ** flatchain[:, logpar]
        bestfit_pars_real[logpar] = 10. ** bestfit_pars[logpar]

    return flatchain_real, bestfit_pars_real


def do_plot_mcmc(flatchain_real, bestfit_pars_real, mcmc_fig_filename,
                 **kwargs):
    import corner

    # make corner plot
    print("Making corner plot...")
    fig = corner.corner(flatchain_real, truths=bestfit_pars_real, **kwargs)
    fig.savefig(mcmc_fig_filename)
    print("done")


def print_fit_results(flatchain_real, bestfit_pars_real, starname, labels,
                      mcmc_results_filename):
    from astropy.table import Table

    # write table of results
    percentiles = compute_percentiles(flatchain_real)
    t = Table(names=('field', 'info', 'bestfit', '16th', '50th', '84th'),
              dtype=('S15', 'S15', 'float64', 'float64', 'float64', 'float64'))
    t.add_row(('starname', starname, np.nan, np.nan, np.nan, np.nan))
    for i in xrange(len(bestfit_pars_real)):
        t.add_row(('par{0}'.format(i), labels[i], bestfit_pars_real[i],
                   percentiles[i, 0], percentiles[i, 1], percentiles[i, 2]))

    print(repr(t))
    t.write(mcmc_results_filename, format='ascii.fixed_width_two_line')


def do_bestfit_model(bestfit_pars, headers, bestfit_dir,
                     bestfit_model_filename, bestfit_uvmaps_filename, **kwargs):
    from fit_main import setup

    print("Bestfit parameters:\n{0} \n".format(repr(bestfit_pars)))

    obs, model, prob, fitter, imagers, = setup(headers)
    ref_coords = prob.get_parameters(bestfit_pars)

    # to compute the bestfit model (and residuals) we must shift the model in the disk center
    ref_coords['delta_delta'] = [0.]
    ref_coords['delta_alpha'] = [0.]

    analysis.compute_models(None, headers, model, prob,
                                bestfit_model_filename, ncpus=1,
                                single_model_pars=bestfit_pars)
    analysis.compute_visibilities(obs, imagers, ref_coords, bestfit_model_filename,
                                      bestfit_uvmaps_filename, working_dir=bestfit_dir,
                                      ncpus=1, **kwargs)


def do_bestfit_disk_structure(bestfit_pars, headers,
                              bestfit_structure_filename):
    from fit_main import setup
    from py2layer.utils import print_model_structure

    obs, model, prob, fitter, imagers = setup(headers)
    res = model.compute(*prob.get_parameters(bestfit_pars)['model'],
                        return_opacity=True)
    print_model_structure(res, model, bestfit_structure_filename)


def do_plot_uv_bestfit(bestfit_pars, headers, uvbinsize, uvlim, Jylims,
                       uvdata_dir, bestfit_dir, plot_filename, **kwargs):
    # TODO: generalize for nwle>1
    from fit_main import setup

    obs, model, prob, fitter, imagers, = setup(headers)

    ref_coords = prob.get_parameters(bestfit_pars)
    inc = np.array([ref_coords['inc'], ref_coords['inc']])
    PA = np.array([ref_coords['PA'], ref_coords['PA']])

    uvtables = [os.path.join(uvdata_dir, headers['obs']['data_filenames'][0]),
                os.path.join(bestfit_dir, headers['obs']['data_filenames'][0][:-4] + "_mod.txt")]

    unit_multiplier = 1.
    wavelength = [headers['obs']['wle_mm'], headers['obs']['wle_mm']]
    ismodel = [False, True]
    colors = ['k', 'b']
    linestyle = ['.', '-']

    apply_shift = dict(delta_alpha=[ref_coords['delta_alpha'][0],
                                    ref_coords['delta_alpha'][0]],
                       delta_delta=[ref_coords['delta_delta'][0],
                                    ref_coords['delta_delta'][0]])

    from pyvfit.utils.util_plot import plot_uvplot
    fig = plot_uvplot(uvtables, [uvbinsize[0], uvbinsize[0]], wavelength, inc,
                      PA, np.array(Jylims[0]) * unit_multiplier,
                      np.array(Jylims[0][1]) * unit_multiplier, linestyle,
                      colors, uvlim=uvlim[0], ismodel=ismodel, fontsize=30,
                      fontsize_ticklabels=20, unit_multiplier=unit_multiplier,
                      apply_shift=apply_shift)

    ax = fig.get_axes()
    plot_label = kwargs.get("inset_caption", "")
    ax[0].text(uvlim[0][1] * 0.85, ax[0].get_ylim()[1] * unit_multiplier * 0.85,
               plot_label, fontsize=30, fontweight='bold', ha='right')
    for tick in ax[1].xaxis.get_major_ticks():
        tick.label.set_fontsize(25)
    for tick in ax[1].yaxis.get_major_ticks():
        tick.label.set_fontsize(25)
    for tick in ax[0].yaxis.get_major_ticks():
        tick.label.set_fontsize(25)

    fig.tight_layout()
    fig.savefig(plot_filename)


def do_models(flatchain, headers, models_filename, **kwargs):
    from fit_main import setup

    # how to select the models to be computed?  E.g. solved using NUM_STEPS/10
    obs, model, prob, fitter, imagers = setup(headers)
    analysis.compute_models(flatchain, headers, model, prob,
                                models_filename, **kwargs)


def do_visibilities(bestfit_pars, headers, models_filename, uvmaps_filename, **kwargs):
    from fit_main import setup

    obs, model, prob, fitter, imagers = setup(headers)
    ref_coords = prob.get_parameters(bestfit_pars)
    analysis.compute_visibilities(obs, imagers, ref_coords, models_filename,
                                      uvmaps_filename, **kwargs)


def do_plot_uv_density(bestfit_pars, uvmaps_filename, uvlim, plots_dir, plot_uv_density_basename, **kwargs):
    from fit_main import setup

    loaded = util.load_models(uvmaps_filename)
    uvmaps = loaded['models']
    headers = loaded['headers']

    obs, model, prob, fitter, imagers, = setup(headers)

    ref_coords = prob.get_parameters(bestfit_pars)
    util_plot.plot_visibilities(uvmaps, headers, obs,
                                headers['analysis']['uvbinsize'],
                                headers['analysis']['uvdist_bin'],
                                plots_dir, plot_uv_density_basename,
                                ref_coords, uvlim,
                                export_ascii=False, **kwargs)


def do_compute_physical_quantities(models_filename, bestfit_structure_filename, results_filename,
                                   gas_to_dust=100., cumul_thrs=0.95):
    """
    Adds some physical quantities of the models to the results_#.txt file.

    :return:
    """
    # TODO: compute the average temperature, e.g. <T_mid> = (T_mid(R) 2piRcos(i) dR)/(int 2piR cos(i) dR)
    from astropy.io import ascii
    from scipy.integrate import cumtrapz
    from pyvfit.utils import util
    from pyvfit.constants import AU, pc

    mJy = 1.e-26  # erg/s/cm2/Hz

    # open the models file
    loaded = util.load_models(models_filename)
    headers = loaded['headers']
    models = loaded['models']
    nmodels = len(models)

    # distance to the star
    dist2 = (headers['star']['dist']*pc)**2.

    m_dust = np.zeros(nmodels)
    flux_int = np.zeros(nmodels)
    r_out = np.zeros(nmodels)
    temp_mid = np.zeros(nmodels)

    for i, res in enumerate(models):
        cos_inc = np.cos(res['inc'])

        # cumulative mass
        sigma_cumul = 2. * np.pi * cumtrapz(res['gridrad'] * res['surf_dens'],
                                            res['gridrad'], initial=0.)

        temp_avg = 2. * np.pi * cumtrapz(res['gridrad'] * res['surf_dens']*res['temp_mid'],
                                            res['gridrad'], initial=0.)/sigma_cumul
    
        if not np.isnan(res['intensity']).any():
            # cumulative integrated flux
            flux_cumul = 2. * np.pi * cos_inc / dist2 / mJy * cumtrapz(
                res['gridrad'] * res['intensity'][0], res['gridrad'], initial=0.)

        else:
            print("WARNING: NaN found in model {0}, parameters: {1}".format(i, res['parameters']))
            continue

        flux_cumul_norm = flux_cumul / flux_cumul[-1]

        # outer radius: radius that contains cumul_thrs% of the integrated flux
        flux_iout = np.where(flux_cumul_norm >= cumul_thrs)[0][0]

        m_dust[i] = sigma_cumul[flux_iout] / gas_to_dust
        flux_int[i] = flux_cumul[flux_iout]
        r_out[i] = res['gridrad'][flux_iout]/AU
        temp_mid[i] = temp_avg[flux_iout]

    # compute percentiles
    q = [16., 50., 84.]
    m_dust_perc = np.percentile(m_dust, q=q)
    flux_int_perc = np.percentile(flux_int, q=q)
    r_out_perc = np.percentile(r_out, q=q)
    temp_mid_perc = np.percentile(temp_mid, q=q)

    try:
        # compute bestfit model quantities
        w = ascii.read(bestfit_structure_filename)

        # _bf : _bestfit
        flux_cumul_bf = np.array(w['I0_cumul'])
        flux_cumul_norm_bf = flux_cumul_bf / flux_cumul_bf[-1]
        r_bf = np.array(w['R'])
        sigma_cumul_bf = np.array(w['Sigma_cumul'])

        temp_avg_bf = 2. * np.pi * cumtrapz(w['R']*AU * w['Sigma'] * w['T_mid'],
                                            w['R']*AU, initial=0.)/sigma_cumul_bf
        
        flux_iout_bf = np.where(flux_cumul_norm_bf >= cumul_thrs)[0][0]

        m_dust_bf = sigma_cumul_bf[flux_iout_bf] / gas_to_dust
        flux_int_bf = flux_cumul_bf[flux_iout_bf]
        r_out_bf = r_bf[flux_iout_bf]
        temp_mid_bf = temp_avg_bf[flux_iout_bf]

    except IOError:
        print("{0} not found. Bestfit model not available.".format(bestfit_structure_filename))
        print("Execute command 3: bestfit_model structure.")
        r_out_bf = np.nan
        m_dust_bf = np.nan
        flux_int_bf = np.nan

    t = ascii.read(results_filename)
    t.add_row(['r_out', 'AU', r_out_bf, r_out_perc[0], r_out_perc[1], r_out_perc[2]])
    t.add_row(['m_dust', 'g', m_dust_bf, m_dust_perc[0], m_dust_perc[1], m_dust_perc[2]])
    t.add_row(['flux_int', 'mJy', flux_int_bf, flux_int_perc[0], flux_int_perc[1], flux_int_perc[2]])
    t.add_row(['temp_mid', 'K', temp_mid_bf, temp_mid_perc[0], temp_mid_perc[1], temp_mid_perc[2]])

    print("Updated table in {0}:".format(results_filename))
    print(repr(t))
    t.write(results_filename, format='ascii.ipac')

    # t['T_mid'] = res['temp_mid']
    # t['T_sur'] = res['temp_sur']
    # t['H_sur'] = res['hsurf'] * res['gridrad'] / AU
    # t['kappa_0'] = res['opacity'][0]


def search_files(directory, pattern):
    import glob

    files = glob.glob(os.path.join(directory, pattern))
    if len(files) == 0:
        raise IOError("No {0} found in {1}".format(pattern, directory))

    print("{0} found in {1}".format(files, directory))

    return files


def do_prepare_maps(original_ms_filename, maps_dir, bestfit_dir):
    """
    Prepare the scripts to be executed in CASA.

    """
    mod_files = search_files(bestfit_dir, '*_mod.txt')
    for f in mod_files:
        print("cp {0} {1}".format(f, maps_dir))
        call("cp {0} {1}".format(f, maps_dir), shell=True)

    mod_basename = os.path.basename(mod_files[0]).split('.')[0]
    # TODO: update for multi wavelength fits

    res_files = search_files(bestfit_dir, '*_res.txt')
    for f in res_files:
        print("cp {0} {1}".format(f, maps_dir))
        call("cp {0} {1}".format(f, maps_dir), shell=True)

    res_basename = os.path.basename(res_files[0]).split('.')[0]

    try:
        msfiles = search_files(maps_dir, original_ms_filename)
    except IOError:
        print("WARNING: MS table {0} not found in {1}. "
              "Proceeding anyway.".format(original_ms_filename, maps_dir))
        msfiles = [original_ms_filename]
        # ans = raw_input("Do you wish to proceed anyway?")
    ms_basename = os.path.basename(msfiles[0]).split('.')[0]

    fill_script_filename = os.path.join(maps_dir, "do_fill_ms.py")
    with open(fill_script_filename, 'w') as f:
        f.write(textwrap.dedent(
        """
        # Automatically generated script to fill the ms tables.
        # Marco Tazzari, ESO.

        # Execute from the CASA shell with:
        #   execfile("do_fill_ms.py")

        import os
        import time
        import numpy as np

        #  Function that reads the original ms table and writes output tables for the models and residual
        #  NB. No checks done
        def write_ascii_to_ms(origdata_ms,mod_txt,mod_ms,res_txt,res_ms,tb,dualpol=True,dores=True):
            print " "
            print "Create CASA ms with model ascii uvtable"
            print " "


            # Read ASCII file
            start=time.clock()

            f_in = open(mod_txt, "r")

            # Initialize arrays
            um = []; vm = []; Rem = []; Imm = []; Wem = [];

            # Go over all lines in ascii file and store data
            for line in f_in:

              # Store data in arrays
              variable = line.strip().split()

              um.append(float(variable[0]))
              vm.append(float(variable[1]))
              Rem.append(float(variable[2]))
              Imm.append(float(variable[3]))
              Wem.append(float(variable[4]))

            # Close ascii file
            f_in.close()

            if dores:
              f_in = open(res_txt, "r")

              # Initialize arrays
              ur = []; vr = []; Rer = []; Imr= []; Wer = [];

              # Go over all lines in ascii file and store data
              for line in f_in:

                # Store data in arrays
                variable = line.strip().split()

                ur.append(float(variable[0]))
                vr.append(float(variable[1]))
                Rer.append(float(variable[2]))
                Imr.append(float(variable[3]))
                Wer.append(float(variable[4]))

              # Close ascii file
              f_in.close()

            print "Finished reading model visibilities..."
            end = time.clock()
            print "It took", end-start, "seconds"

            # Read ms table and replace data column
            start=time.clock()


            # Create new tables
            syscommand = 'rm -rf '+ mod_ms
            os.system(syscommand)
            syscommand = 'cp -r '+ origdata_ms + ' ' + mod_ms
            os.system(syscommand)

            if dores:
              syscommand = 'rm -rf '+ res_ms
              os.system(syscommand)
              syscommand = 'cp -r '+ origdata_ms + ' ' + res_ms
              os.system(syscommand)

            # Open ms table
            tb.open(mod_ms)

            # Check if CORRECTED_DATA column is present:
            # (MODEL_DATA column will be modified when cleaning)
            all_columns = tb.colnames()
            correct     = False
            if 'CORRECTED_DATA' in all_columns:
              correct_data_m = tb.getcol("CORRECTED_DATA")
              correct = True

            # Get column UVW, DATA, WEIGHT and DATA_DESC_ID (spw info!)
            data_m   = tb.getcol("DATA")
            uvw_m    = tb.getcol("UVW")
            weight_m = tb.getcol("WEIGHT")
            spw_m    = tb.getcol("DATA_DESC_ID")
            tb.close()

            if dores:
              tb.open(res_ms)

              # Check if CORRECTED_DATA column is present:
              # (MODEL_DATA column will be modified when cleaning)
              all_columns = tb.colnames()
              correct     = False
              if 'CORRECTED_DATA' in all_columns:
                correct_data_r = tb.getcol("CORRECTED_DATA")
                correct = True

              # Get column UVW, DATA, WEIGHT and DATA_DESC_ID (spw info!)
              data_r   = tb.getcol("DATA")
              uvw_r    = tb.getcol("UVW")
              weight_r = tb.getcol("WEIGHT")
              spw_r    = tb.getcol("DATA_DESC_ID")
              tb.close()


            # Go over msdata array and input new model values into
            # the XX (data[0,0,:]) and YY (data[1,0,:]) components:
            print "Size of current ms file: ", len(spw_m)
            for i in range(len(spw_m)):

              data_m[0,0,i] = complex(Rem[i],Imm[i])
              weight_m[0,i] = Wem[i]
              if dualpol:
                data_m[1,0,i] = complex(Rem[i],Imm[i])
                weight_m[1,i] = Wem[i]
              if correct is True:
                correct_data_m[0,0,i] = complex(Rem[i],Imm[i])
                if dualpol:
                  correct_data_m[1,0,i] = complex(Rem[i],Imm[i])
              if dores:
                data_r[0,0,i] = complex(Rer[i],Imr[i])
                weight_r[0,i] = Wer[i]
                if dualpol:
                  data_r[1,0,i] = complex(Rer[i],Imr[i])
                  weight_r[1,i] = Wer[i]
                if correct is True:
                  correct_data_r[0,0,i] = complex(Rer[i],Imr[i])
                  if dualpol:
                    correct_data_r[1,0,i] = complex(Rer[i],Imr[i])

            # Put back modified weight column on original ms

            # Put back modified column
            tb.open(mod_ms,nomodify=False)
            tb.putcol("DATA",data_m)
            tb.putcol("WEIGHT",weight_m)
            if correct is True:
              tb.putcol("CORRECTED_DATA",correct_data_m)
            tb.flush()
            # Close table tools
            tb.close()

            if dores:
              tb.open(res_ms,nomodify=False)
              tb.putcol("DATA",data_r)
              tb.putcol("WEIGHT",weight_r)
              if correct is True:
                tb.putcol("CORRECTED_DATA",correct_data_r)
              tb.flush()

            # Close table tools
            tb.close()

            print "Finished with model     ms: ", mod_ms
            if dores:
              print "     and with residuals ms: ", res_ms
            end = time.clock()
            print "It took", end-start, "seconds"
            print ""

        write_ascii_to_ms("{0}.ms",
                              "{1}.txt", "{1}.ms",
                              "{2}.txt", "{2}.ms", tb)
        """.format(ms_basename, mod_basename, res_basename)))
    print("{0} has been written correctly".format(fill_script_filename))

    clean_script_filename = os.path.join(maps_dir, "do_clean.py")
    with open(clean_script_filename, 'w') as f:
        f.write(textwrap.dedent(
        """
        # Automatically generated script to perform the cleaning.
        # Marco Tazzari, ESO.

        # Execute from the CASA shell with:
        #   execfile("do_clean.py")

        mythres = "20e-6mJy"
        niter = 1000
        mycell = '0.025arcsec'
        mysize = [300, 300]
        interactive = True

        # obs
        clean(vis="{0}.ms", imagename="{0}", selectdata=False, interactive=interactive,
              niter=niter, threshold=mythres, cell=mycell, imsize=mysize)
        # model
        clean(vis="{1}.ms", imagename="{1}", selectdata=False, interactive=interactive,
              niter=niter, threshold=mythres, cell=mycell, imsize=mysize)
        # residuals
        clean(vis="{2}.ms", imagename="{2}", selectdata=False, interactive=interactive,
              niter=niter, threshold=mythres, cell=mycell, imsize=mysize)
        """.format(ms_basename, mod_basename, res_basename)))
    print("{0} has been written correctly".format(clean_script_filename))

    export_script_filename = os.path.join(maps_dir, "do_exportfits.py")
    with open(export_script_filename, 'w') as f:
            f.write(textwrap.dedent(
        """
        # Automatically generated script to export the cleaned images to fits.
        # Marco Tazzari, ESO.

        # Execute from the CASA shell with:
        #   execfile("do_exportfits.py")

        # export images to fits
        exportfits(imagename="{0}.image", fitsimage="{0}.fits")
        exportfits(imagename="{1}.image", fitsimage="{1}.fits")
        exportfits(imagename="{2}.image", fitsimage="{2}.fits")
        """.format(ms_basename, mod_basename, res_basename)))
    print("{0} has been written correctly".format(export_script_filename))


# if do_plot_surface_density:
#     call("mkdir -p " + PLOTS_DIR, shell=True)
#
#     loaded = util.load_models(models_filename)
#     models = loaded['models']
#
#     figure = util_plot.plot_physical_quantity(models, 'surf_dens', True,
#                                               [-2., 2.],
#                                               ylabel=r'$\Sigma(R)$ (g/cm$^{2}$)')
#     figure.savefig(PLOTS_DIR + "_surface_density.pdf")
#
# if do_plot_temp_mid:
#     call("mkdir -p " + PLOTS_DIR, shell=True)
#
#     loaded = util.load_models(models_filename)
#     models = loaded['models']
#     figure = util_plot.plot_physical_quantity(models, 'temp_mid', True,
#                                               [0.93, 2], loglin=True,
#                                               ylabel=r'$T_{\mathrm{mid}}(R)$ (K)',
#                                               yticks=[1., 1.301, 1.477, 1.6989,
#                                                       1.845, 2.],
#                                               yticks_labels=[10, 20, 30, 50, 70,
#                                                              100],
#                                               nbins=(200, 300))
#
#     figure.savefig(PLOTS_DIR + "_temp_mid.pdf")
#
#
# if do_2d_contours:
#
#     # flatchain = util_chain.chain_thin(chain.copy(), start_mcmc, stop_mcmc, step_mcmc, logpars)
#     from corner import hist2d
#     import matplotlib as mpl
#
#     import matplotlib.pyplot as plt
#
#     # print(mpl.matplotlib_fname())  # matplotlibrc
#     mpl.rcParams['axes.linewidth'] = 1.7
#     mpl.rcParams['xtick.major.size'] = 14
#     mpl.rcParams['xtick.minor.size'] = 7
#     mpl.rcParams['xtick.major.width'] = 1.7
#     mpl.rcParams['xtick.minor.width'] = 1.2
#     mpl.rcParams['xtick.major.pad'] = 8
#     mpl.rcParams['xtick.labelsize'] = 'large'
#     mpl.rcParams['ytick.major.size'] = 14
#     mpl.rcParams['ytick.minor.size'] = 7
#     mpl.rcParams['ytick.major.width'] = 1.7
#     mpl.rcParams['ytick.minor.width'] = 1.2
#     mpl.rcParams['ytick.major.pad'] = 8
#     mpl.rcParams['ytick.labelsize'] = 'large'
#     mpl.rcParams['font.family'] = 'sans-serif'
#     mpl.rcParams['font.weight'] = 500
#     mpl.rcParams['font.size'] = 16
#     # if we want to plot gamma-Rc and p-Rout in two panels
#     # from matplotlib import gridspec
#     # fig = plt.figure(figsize=(10, 5))
#     # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
#     # ax = plt.subplot(gs[0]), plt.subplot(gs[1])
#
#     fig, ax = plt.subplots(1, 1, figsize=(5, 5))
#     levels = 1. - np.exp(-0.5 * np.arange(1., 3.1, 1) ** 2.)
#     # correspond to 1, 2 and 3 sigma for a 2D hystogram
#
#     hist2d(flatchain[:, 0], flatchain[:, 2], ax=ax, bins=30, levels=levels,
#            plot_datapoints=False, smooth=1.5, plot_contours=True,
#            fill_contours=True, plot_density=False, color='#2D89E5')
#
#     # add point for best-fit model
#     ax.plot(BESTFIT_PARS[0], BESTFIT_PARS[2], marker='*',
#             markerfacecolor='#FFFF66', markeredgecolor='#FFFF66', markersize=16)
#
#     if headers['fit']['sigma_prescription'] == 'g':
#         ax.set_xlabel(r'$\gamma$', fontsize=20)
#         ax.set_ylabel(r'$R_c$', fontsize=20)
#     elif headers['fit']['sigma_prescription'] == 'pl':
#         ax.set_xlabel(r'$p$', fontsize=20)
#         ax.set_ylabel(r'$R_{\mathrm{out}}$', fontsize=20)
#
#     # set axes limits
#     ax.set_ylim(0, 300.)
#     yticks = np.arange(0, 301, 100, dtype='int')
#     ax.set_yticks(yticks)
#     ax.set_yticklabels(yticks)
#     if headers['fit']['sigma_prescription'] == 'g':
#         xticks = np.arange(-2, 2.1, 1, dtype='int')
#         ax.set_xlim(-2., 2.)
#         ax.set_xticks(xticks)
#         ax.set_xticklabels(xticks)
#     elif headers['fit']['sigma_prescription'] == 'pl':
#         xticks = np.arange(-2, 4.1, 1, dtype='int')
#         ax.set_xlim(-2., 4.)
#         ax.set_xticks(xticks)
#         ax.set_xticklabels(xticks)
#
#     ax.minorticks_on()
#     sourcelabel_xcoord = (xticks[-1] - xticks[0]) * 0.15 + xticks[0]
#     sourcelabel_ycoord = yticks[-1] * 0.85
#     ax.text(sourcelabel_xcoord, sourcelabel_ycoord, SOURCE_LABEL, fontsize=20,
#             fontweight='bold')
#
#     fig.tight_layout()
#
#     fig.savefig(PLOTS_DIR + "dist_2D_" + STARNAME + ".pdf")


def main():
    from _source import STARNAME, SOURCE_LABEL, ACC_RATIO_THRS, BESTFIT_PARS, \
        NUM_STEPS, NUM_MODELS, get_bestfit_model, mcmc_range, Jylims, Jyticks, \
        Jyticklabels, Jyunit, UVLIM, UVBINSIZE, ORIGINAL_MS_FILENAME

    try:
        from _source import refine_chain
        print("INFO: refine_chain imported from _source.py.")
    except ImportError:
        print("WARNING: refine_chain not found in _source.py. Continuing anyway.")

        def refine_chain(flatchain):
            return flatchain

    # log_format = '%(asctime)-20s %(levelname)-8s %(name)-8s %(message)s'
    # logfilename = 'analysis_' + time.strftime("%Y_%m_%d_%H%M%S",
    #                                        time.localtime()) + '.log'  # for testing purposes
    # logging.basicConfig(filename=logfilename, filemode='w', format=log_format,
    #                     level=logging.INFO)
    # console.setFormatter(logging.Formatter(log_format))

    print("Starting analysis of {0}".format(STARNAME))

    # directories
    WORKING_DIR = './'
    SNAPSHOTS_DIR = os.path.join(WORKING_DIR, "snapshots/")
    UVDATA_DIR = os.path.join(WORKING_DIR, "uvdata")
    ANALYSIS_DIR = os.path.join(WORKING_DIR, 'analysis/')

    PLOTS_DIR = os.path.join(ANALYSIS_DIR, "plots/")
    BESTFIT_DIR = os.path.join(ANALYSIS_DIR, "bestfit/")
    CHAIN_DIAG_DIR = os.path.join(ANALYSIS_DIR, "chain_diag/")
    MAPS_DIR = os.path.join(ANALYSIS_DIR, "maps/")

    call("mkdir -p " + ANALYSIS_DIR, shell=True)
    call("mkdir -p " + PLOTS_DIR, shell=True)
    call("mkdir -p " + BESTFIT_DIR, shell=True)
    call("mkdir -p " + CHAIN_DIAG_DIR, shell=True)
    call("mkdir -p " + MAPS_DIR, shell=True)

    # filenames
    mcmc_fig_filename = os.path.join(PLOTS_DIR, "triangle_{0}.pdf".format(STARNAME))
    mcmc_results_filename = os.path.join(ANALYSIS_DIR, "results_{0}.txt".format(STARNAME))
    bestfit_model_filename = os.path.join(BESTFIT_DIR, 'bestfit_model_{0}.dat'.format(STARNAME))
    bestfit_uvmaps_filename = os.path.join(BESTFIT_DIR, 'bestfit_uvmaps.dat')
    bestfit_structure_filename = os.path.join(BESTFIT_DIR, "bestfit_disk_structure_{0}.txt".format(STARNAME))
    plot_uv_filename = os.path.join(PLOTS_DIR, "bestfit_uvplot_{0}.pdf".format(STARNAME))
    models_filename = os.path.join(ANALYSIS_DIR, 'models.dat')
    uvmaps_filename = os.path.join(ANALYSIS_DIR, 'uvmaps.dat')
    plot_uv_density_basename = "uvplot_{0}".format(STARNAME)

    # analysis - chain parameters
    filename_snaphots = get_last_snapshot(SNAPSHOTS_DIR)
    snapshot = util.load_fit(filename_snaphots)
    headers = snapshot['headers']

    assert STARNAME == headers['star']['starname'], \
        "Starname given in _source is {0}, while in the headers is {1}".format(
        STARNAME, headers['star']['starname'])

    logpars = headers['fit'].get('logpar', [1])
    last_step = int(headers['results']['current_step'])
    start_mcmc, stop_mcmc, step_mcmc = max(last_step - NUM_STEPS,
                                           0), last_step, 1
    print("Considering MCMC steps from {0} to {1}, step {2}".format(start_mcmc,
                                                                    stop_mcmc,
                                                                    step_mcmc))

    # select steps in the chain
    raw_chain, headers, raw_blobs = get_chain(filename_snaphots)
    print("Chain and blobs shape (raw): {0}, {1}".format(raw_chain.shape, raw_blobs.shape))
    raw_chain = raw_chain[:, start_mcmc:stop_mcmc, :]
    raw_blobs = raw_blobs[:, start_mcmc:stop_mcmc]
    print("Chain and blobs shape (raw) after step selection: {0}, {1}".format(raw_chain.shape,
          raw_blobs.shape))
    print("Autocorrelation time:\n{0}\n\n".format(
        repr(headers['results']['autocorrelation_time'])))

    # clean the chain
    chain, blobs = clean_chain(raw_chain, headers, raw_blobs,
                               acc_ratio_thrs=ACC_RATIO_THRS,
                               chain_diag_dir=CHAIN_DIAG_DIR,
                               savefigs=chain_diagnostics)
    print("Chain and blobs shape (clean): {0}, {1}".format(chain.shape, blobs.shape))

    # extract the flatchain
    nwalkers, nsteps, npars = chain.shape
    flatchain = chain.reshape(nwalkers * nsteps, npars).copy()
    print("Flatchain shape after cleaning: {0}".format(flatchain.shape))

    if BESTFIT_PARS is None:
        BESTFIT_PARS = get_bestfit_model(flatchain.copy(), blobs)
    print("Bestfit model parameters (log): {0}\n".format(repr(BESTFIT_PARS)))

    flatchain = refine_chain(flatchain.copy())
    print("Flatchain shape after refining: {0}".format(flatchain.shape))

    flags = list(flags_str)

    while len(flags) > 0:
        flag = int(flags.pop(0))
        print(textwrap.dedent(
        """
        ****************************************************
        ****
        ****    {0}. {1}
        ****
        ****************************************************
        """.format(flag, tasks[flag].upper())))

        # convert to the real space
        if flag == 0 or flag == 1 or flag == 8:
            if logpars:
                flatchain_real, bestfit_pars_real = logtoreal(flatchain,
                                                              BESTFIT_PARS,
                                                              logpars)
            print("Bestfit model parameters (real): {0}\n".format(repr(bestfit_pars_real)))

        if flag == 0:
            do_plot_mcmc(flatchain_real.copy(), bestfit_pars_real.copy(), mcmc_fig_filename,
                         range=mcmc_range, labels=headers['fit']['par_names'],
                         labels_fontsize=20, verbose=False,
                         quantiles=[0.16, 0.50, 0.84], plot_datapoints=False,
                         plot_contours=True, show_titles=False,
                         title_args={"fontsize": 16}, smooth=0.6,
                         truth_color='r', bins=25)

        if flag == 1:
            print_fit_results(flatchain_real.copy(), bestfit_pars_real.copy(), STARNAME,
                              headers['fit']['par_names'],
                              mcmc_results_filename)

        if flag == 2:
            do_bestfit_model(BESTFIT_PARS, headers, BESTFIT_DIR,
                             bestfit_model_filename, bestfit_uvmaps_filename, uvbinsize=UVBINSIZE)

        if flag == 3:
            do_bestfit_disk_structure(BESTFIT_PARS, headers,
                                      bestfit_structure_filename)

        if flag == 4:
            do_plot_uv_bestfit(BESTFIT_PARS, headers, UVBINSIZE, UVLIM, Jylims,
                               UVDATA_DIR, BESTFIT_DIR, plot_uv_filename,
                               inset_caption=SOURCE_LABEL)
            # TODO: use Jyticks, Jyticklabels

        if flag == 5:
            flatchain_analysis = flatchain[np.random.choice(flatchain.shape[0], NUM_MODELS), :]

            do_models(flatchain_analysis.copy(), headers, models_filename, ncpus=ncpus)

        if flag == 6:
            do_visibilities(BESTFIT_PARS, headers, models_filename,
                            uvmaps_filename, uvbinsize=UVBINSIZE, ncpus=ncpus)

        if flag == 7:
            inset_caption = [SOURCE_LABEL for wl in UVBINSIZE]
            do_plot_uv_density(BESTFIT_PARS, uvmaps_filename, UVLIM, PLOTS_DIR,
                               plot_uv_density_basename, Jylims=Jylims,
                                Jyunit=Jyunit, Jyticks=Jyticks,
                                Jyticklabels=Jyticklabels, inset_caption=inset_caption)

        if flag == 8:
            print_fit_results(flatchain_real.copy(), bestfit_pars_real.copy(),
                              STARNAME, headers['fit']['par_names'],
                              mcmc_results_filename)

            do_compute_physical_quantities(models_filename, bestfit_structure_filename,
                                           mcmc_results_filename)

        if flag == 9:
            do_prepare_maps(ORIGINAL_MS_FILENAME, MAPS_DIR, BESTFIT_DIR)
            # TODO: copy _mod.txt _res.txt to the maps/directory

            # TODO: execute the CASA scripts


if __name__ == '__main__':
    main()
