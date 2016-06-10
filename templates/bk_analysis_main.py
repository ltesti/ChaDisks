#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import pyvfit
from pyvfit.utils import util_plot, util, util_chain
from pyvfit import analysis, get_last_snapshot, read_headers
from pyvfit.constants import AU

import numpy as np
from subprocess import call
import optparse
import sys
import os

import corner
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 16


def get_peaks(flatchain):
    # use peaks
    npars = flatchain.shape[1]
    bestfit_pars = np.zeros(npars)
    for i_par in xrange(npars):
        n, bins = np.histogram(chain[:, :, i_par].flatten(), bins=200)
        bestfit_pars[i_par] = bins[np.argmax(n)]

    return bestfit_pars

def get_min_chisquare(flatchain, blobs):

    return flatchain[np.argmin(blobs), :]


def clean_chain(chain, h, blobs, acc_ratio_thrs=None, savefigs=False, **kwargs):
    """
    Flags the walkers with acceptance ratio lower than acc_ratio_thrs.

    """
    figsize = kwargs.get('figsize', (11., 7.))
    working_dir = kwargs.get('chain_diag_dir', './')
    nwalkers = h['fit']['nwalkers']
    par_names = h['fit']['par_names']

    #### 2. Trace plots for some walkers
    steps = np.arange(chain.shape[1])
    npars = chain.shape[2]
    acc_ratio = h['results']['acceptance_ratio']

    if not acc_ratio_thrs:
        acc_ratio_thrs = np.percentile(acc_ratio, [50.])
        print(acc_ratio_thrs)

    flagged = acc_ratio < acc_ratio_thrs
    not_flagged = np.logical_not(flagged)
    print("Walkers left after cleaning: {0}".format(len(acc_ratio[not_flagged])))

    if savefigs:
        # ** Trace plots of all walkers **
        print("Plotting raw traces for parameters: ", end='')
        for i_par in xrange(npars):
            plt.figure(figsize=figsize)
            print(" {0}".format(i_par), end='')
            for walker in xrange(nwalkers):
                plt.plot(steps, chain[walker, :, i_par], alpha=0.5)
                plt.ylabel(par_names[i_par])
            plt.tight_layout()
            plt.savefig(os.path.join(working_dir, "trace_raw_{0}.png".format(i_par)))

        print("", end="\n")

        plt.figure(figsize=figsize)
        for iw in xrange(nwalkers):
            plt.plot(steps, blobs[iw, :], alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(working_dir, "trace_chi_square.png"))

        # histogram of acceptance_ratio of all walkers
        plt.figure(figsize=figsize)
        plt.hist(acc_ratio, bins=100)
        plt.savefig(os.path.join(working_dir, "acc_ratio_raw.png"))

        # histogram of acceptance_ratio of cleaned walkers
        plt.figure(figsize=figsize)
        plt.hist(acc_ratio[not_flagged])
        plt.savefig(os.path.join(working_dir, "acc_ratio_cleaned.png"))

        # ** Trace plots of not_flagged walkers **
        walkers_not_flagged = np.arange(0, nwalkers, 1)[not_flagged]
        print("Plotting cleaned traces for parameters: ", end='')
        for i_par in xrange(npars):
            plt.figure(figsize=figsize)
            print(" {0}".format(i_par), end='')
            for walker in walkers_not_flagged:
                plt.plot(steps, chain[walker, :, i_par], alpha=0.5)
                plt.ylabel(par_names[i_par])
            plt.tight_layout()
            plt.savefig(os.path.join(working_dir, "trace_clean_{0}.png".format(i_par)))

        print("", end="\n")

    return chain[not_flagged, :, :], blobs[not_flagged, :]


def compute_cumulative(x, y):
    nx = x.shape[0]
    assert nx == y.shape[0]

    dx = np.zeros(nx)
    for i in xrange(nx-1):
        dx[i] = x[i+1]-x[i]
    dx[-1] = dx[-2]

    return np.cumsum(x*y*dx)


from _object import STARNAME, SOURCE_LABEL, ACC_RATIO_THRS, bestfit_pars

print("Starting analysis of {0}".format(STARNAME))

# directories
working_dir = './'
snapshots_dir = os.path.join(working_dir, "snapshots/")
uvdata_dir = os.path.join(working_dir, "uvdata")

analysis_dir = os.path.join(working_dir, 'analysis/')
plots_dir = os.path.join(analysis_dir, "plots/")
bestfit_dir = os.path.join(analysis_dir, "bestfit/")
chain_diag_dir = os.path.join(analysis_dir, "chain_diag/")

call("mkdir -p " + analysis_dir, shell=True)
call("mkdir -p " + plots_dir, shell=True)
call("mkdir -p " + bestfit_dir, shell=True)
call("mkdir -p " + chain_diag_dir, shell=True)

# analysis - chain parameters
fit_filename = "fit_" + STARNAME + ".dat"
filename_snaphots = get_last_snapshot(snapshots_dir)

snapshot =  util.load_fit(filename_snaphots)
headers = snapshot['headers']
chain = snapshot['chain']
blobs = snapshot['blobs']
print("Shape of read-in chain: ({0}, {1}, {2})".format(*chain.shape))

logpars = headers['fit'].get('logpar', [1])
last_step = int(headers['results']['current_step'])
assert last_step == chain.shape[1], "Last step:{0} while chain.shape[1]:{1}".format(last_step, chain.shape[1])
start_mcmc, stop_mcmc, step_mcmc = max(last_step-1500 , 0), last_step, 1
print(start_mcmc, stop_mcmc)
start_analysis, stop_analysis, step_analysis = max(last_step-400, 0), last_step, 100

# analysis - computation parameters
ncpus = 20
models_filename = os.path.join(analysis_dir, 'models.dat')
uvmaps_filename = os.path.join(analysis_dir, 'uvmaps.dat')
bestfit_uvmaps_filename = os.path.join(bestfit_dir, 'bestfit_uvmaps.dat')

# plots
mcmc_chain_extents = None
mcmc_fig_filename = "triangle_" + STARNAME + ".pdf"
nwle = 1
Jylims = [[(-0.1, 0.5), (-0.03, 0.03)]]
Jyticks = [[[0., 0.1, 0.2, 0.3, 0.4, 0.5], [-0.02, 0.02]]]
Jyticklabels = Jyticks
Jyunit = [['(Jy)'], ['(Jy)']]
uvupperlim = [800.]
uvbinsize = [30.e3]


p = optparse.OptionParser()
p.add_option("-s", "--select", action="store",dest="selected_flags", default='10000000')
p.add_option("-d", "--diag", action="store_true",dest="savefigs", default=False)
(options,args) = p.parse_args()

# flags
do_plot_mcmc = False if options.selected_flags[0]=='0' else True
do_models = False if options.selected_flags[1]=='0' else True
do_visibilities = False if options.selected_flags[2]=='0' else True
do_bestfit = False if options.selected_flags[3]=='0' else True
do_plot_uv = False if options.selected_flags[4]=='0' else True
do_2d_contours = False if options.selected_flags[5]=='0' else True
do_bestfit_disk_structure = False if options.selected_flags[6]=='0' else True
do_plot_uv_bestfit = False if options.selected_flags[7]=='0' else True

# imaging: TO DO
do_plot_surface_density = False
do_plot_temp_mid = False

raw_chain, headers, raw_blobs = util_chain.get_chain(filename_snaphots)
print("Chain shape (raw): ", raw_chain.shape, raw_blobs.shape)
raw_chain = raw_chain[:, start_mcmc:stop_mcmc, :]
raw_blobs = raw_blobs[:, start_mcmc:stop_mcmc]
print("Chain shape (raw) after step selection: ", raw_chain.shape, raw_blobs.shape)
# clean chain
chain, blobs = clean_chain(raw_chain, headers, raw_blobs, acc_ratio_thrs=ACC_RATIO_THRS,
                                                          chain_diag_dir=chain_diag_dir,
                                                          savefigs=options.savefigs)
print("Chain shape after cleaning ", chain.shape, blobs.shape)

nwalkers, nsteps, npars = chain.shape
flatchain = chain.reshape(nwalkers*nsteps, npars)
print("Flatchain shape: ", flatchain.shape)

if bestfit_pars is None:
    # bestfit_pars = get_peaks(flatchain)         # choose the distribution peaks
    bestfit_pars = get_min_chisquare(flatchain, blobs)   # choose the min(Chisquare)


if do_plot_mcmc:
    call("mkdir -p "+plots_dir, shell=True)

    # convert to the real space
    bestfit_pars_real = bestfit_pars.copy()
    flatchain_real = flatchain.copy()
    if logpars is not None:
        for logpar in logpars:
            bestfit_pars_real[logpar] = 10.**bestfit_pars[logpar]
            flatchain_real[:, logpar] = 10.**flatchain[:, logpar]

    labels = headers['fit']['par_names']
    print("Bestfit model parameters:\n", repr(bestfit_pars_real))

    # make corner plot
    print("Making corner plot...", end='')
    fig = corner.corner(flatchain_real, labels=labels, range=mcmc_chain_extents,
                        truths=bestfit_pars_real, truth_color='r', bins=25,
                        show_titles=False,
                        title_args={"fontsize": 16},
                        quantiles=[0.16, 0.50, 0.84], labels_fontsize=20,
                        plot_contours=True, verbose=False,
                        smooth=0.6)

    fig.savefig(os.path.join(plots_dir, mcmc_fig_filename))
    print("done")

    print("Autocorrelation time:\n{0}\n".format(repr(headers['results']['autocorrelation_time'])))
    print()

    # write table of results
    percentiles = util_chain.compute_percentiles(flatchain_real)
    from astropy.table import Table
    t = Table(names=('field', 'info', 'bestfit', '16th', '50th', '84th'), dtype=('S15', 'S15', 'float64', 'float64', 'float64', 'float64'))
    t.add_row(('starname', STARNAME, np.nan, np.nan, np.nan, np.nan))
    for i in xrange(npars):
        t.add_row(('par{0}'.format(i), labels[i], bestfit_pars_real[i], percentiles[i, 0], percentiles[i, 1], percentiles[i, 2]))

    print(t)
    t.write(os.path.join(analysis_dir ,"results_{0}.txt".format(STARNAME)), format='ascii.fixed_width_two_line')


if do_bestfit_disk_structure:

    from scipy.integrate import cumtrapz
    from fit_main import setup
    obs, model, prob, fitter, imagers = setup(headers)

    res = model.compute(*prob.get_parameters(bestfit_pars)['model'], return_opacity=True)

    # determine location of 0.88, 1.3, 3.0 mm in gridwle
    iw_0_88mm = np.where(res['gridwle']>=0.880/10.)[0][-1]
    iw_1_3mm = np.where(res['gridwle']>=1.3/10.)[0][-1]
    iw_3mm = np.where(res['gridwle']>=3./10.)[0][-1]

    cos_inc = np.cos(res['inc'])
    dist2 = model.star.dist.cm**2.
    mJy = 1.e-26  # erg/s/cm2/Hz/sr

    from astropy.table import Table
    col_names = ('R', 'Sigma', 'T_mid', 'T_sur', 'H_sur', 'I')
    # t = Table(names=col_names, dtype=('float64', 'float64', 'float64', 'float64', 'float64', 'float64'))
    t = Table()
    t['R'] = res['gridrad']/AU
    t['Sigma'] = res['surf_dens']
    t['Sigma_cumul'] = 2*np.pi*cumtrapz(res['gridrad']*res['surf_dens'], res['gridrad'], initial=0.)
    t['T_mid'] = res['temp_mid']
    t['T_sur'] = res['temp_sur']
    t['H_sur'] = res['hsurf']*res['gridrad']/AU
    t['kappa_0'] = res['opacity'][0]
    t['opt_depth_0_88mm'] = res['surf_dens']*res['opacity_mid'][iw_0_88mm]/cos_inc
    t['opt_depth_1_3mm'] = res['surf_dens']*res['opacity_mid'][iw_1_3mm]/cos_inc
    t['opt_depth_3mm'] = res['surf_dens']*res['opacity_mid'][iw_3mm]/cos_inc
    t['R'].unit = 'AU'
    t['Sigma'].unit = 'g/cm^2'
    t['Sigma_cumul'].unit = 'g'
    t['T_mid'].unit = 'K'
    t['T_sur'].unit = 'K'
    t['H_sur'].unit =  'AU'
    t['opt_depth_0_88mm'].unit = ''
    t['opt_depth_1_3mm'].unit = ''
    t['opt_depth_3mm'].unit = ''

    for iw in xrange(res['intensity'].shape[0]):
        t['I{0}'.format(iw)] = res['intensity'][iw]
        t['I{0}'.format(iw)].unit = 'erg/s/cm2/Hz/sr'

    for iw in xrange(res['intensity'].shape[0]):
        t['I{0}_cumul'.format(iw)] =  2.*np.pi*cos_inc/dist2/mJy*cumtrapz(res['gridrad']*res['intensity'][0], res['gridrad'], initial=0.)
        t['I{0}_cumul'.format(iw)].unit = 'mJy'

    print(t)
    t.write(os.path.join(bestfit_dir ,"bestfit_disk_structure_{0}.txt".format(STARNAME)), format='ascii.ipac') #, delimiter='\t') # 'ascii.fixed_width_two_line',


if do_models:

    from fit_main import setup
    obs, model, prob, fitter, imagers = setup(headers)
    x= chain[:, start_analysis:stop_analysis:step_analysis, :]
    a,b,c,=x.shape
    chain=chain.reshape(a*b, c)

    # flatchain = util_chain.chain_thin(chain.copy(), start_analysis, stop_analysis, step_analysis, None)
    ref_coords = prob.get_parameters(bestfit_pars)
    analysis.compute_models_new(flatchain.copy(), headers, model, prob, models_filename, ncpus=ncpus)

if do_visibilities:

    from fit_main import setup
    obs, model, prob, fitter, imagers = setup(headers)
    ref_coords = prob.get_parameters(bestfit_pars)
    analysis.compute_visibilities_new(models_filename, uvmaps_filename, obs, imagers, ref_coords, ncpus=ncpus, uvbinsize_factor=uvbinsize)


if do_bestfit:

    bestfit_model_filename = 'bestfit_model_{0}'.format(STARNAME)
    call("mkdir -p "+bestfit_dir, shell=True)

    print("Bestfit parameters:\n{0} \n".format(repr(bestfit_pars)))

    from fit_main import setup
    obs, model, prob, fitter, imagers, = setup(headers)
    ref_coords = prob.get_parameters(bestfit_pars)
    # to compute the bestfit model (and the residuals) the model must be shifted in the right disk center
    ref_coords['delta_delta'] = [0.]
    ref_coords['delta_alpha'] = [0.]

    analysis.compute_models_new(None, headers, model, prob, bestfit_model_filename, ncpus=1, single_model_pars=bestfit_pars)
    analysis.compute_visibilities_new(bestfit_model_filename, bestfit_uvmaps_filename, obs, imagers, ref_coords, ncpus=1, uvbinsize_factor=uvbinsize, working_dir=bestfit_dir)


if do_plot_surface_density:
    call("mkdir -p "+plots_dir, shell=True)

    loaded = util.load_models(models_filename)
    models = loaded['models']

    figure = util_plot.plot_physical_quantity(models, 'surf_dens', True, [-2., 2.], ylabel= r'$\Sigma(R)$ (g/cm$^{2}$)')
    figure.savefig(plots_dir + "_surface_density.pdf")


if do_plot_temp_mid:
    call("mkdir -p "+plots_dir, shell=True)

    loaded = util.load_models(models_filename)
    models = loaded['models']
    figure = util_plot.plot_physical_quantity(models, 'temp_mid', True, [0.93, 2], loglin=True, ylabel= r'$T_{\mathrm{mid}}(R)$ (K)',
                                              yticks=[1., 1.301, 1.477, 1.6989, 1.845, 2.], yticks_labels=[10, 20, 30, 50, 70, 100],
                                              nbins=(200, 300))

    figure.savefig(plots_dir + "_temp_mid.pdf")


if do_plot_uv:
    call("mkdir -p "+plots_dir, shell=True)

    loaded = util.load_models(uvmaps_filename)
    models = loaded['models']
    headers = loaded['headers']

    from fit_main import setup
    obs, model, prob, fitter, imagers, = setup(headers)
    print("Bestfit pars from medians: ", bestfit_pars)

    ref_coords = prob.get_parameters(bestfit_pars)
    figure = util_plot.plot_visibilities(models, headers, uvmaps_filename, obs, headers['analysis']['uvbinsize'], headers['analysis']['uvdist_bin'], plots_dir, ref_coords,
                      Jylims=Jylims, Jyunit=Jyunit, Jyticks=Jyticks, Jyticklabels=Jyticklabels, uvupperlim=uvupperlim,
                      export_ascii=False)


if do_plot_uv_bestfit:

    from fit_main import setup
    obs, model, prob, fitter, imagers, = setup(headers)

    ref_coords = prob.get_parameters(bestfit_pars)

    uvtables = [os.path.join(uvdata_dir, headers['obs']['data_filenames'][0]), os.path.join(bestfit_dir, headers['obs']['data_filenames'][0][:-4]+"_mod.txt")]
    unit_multiplier=1.
    wavelength = [headers['obs']['wle_mm'], headers['obs']['wle_mm']]
    ismodel = [False, True]
    colors = ['k', 'b']
    linestyle = ['.', '-']
    Jyticks_im = [-0.003, 0., 0.003]

    inc = np.array([ref_coords['inc'], ref_coords['inc']])
    PA = np.array([ref_coords['PA'], ref_coords['PA']])

    apply_shift = dict(delta_alpha=[ref_coords['delta_alpha'][0], ref_coords['delta_alpha'][0]],
                       delta_delta=[ref_coords['delta_delta'][0], ref_coords['delta_delta'][0]])

    plot_filename = os.path.join(plots_dir, "bestfit_uvplot_{0}.pdf".format(STARNAME))
    from pyvfit.utils.util_plot import plot_uvplot
    fig = plot_uvplot(uvtables, [uvbinsize[0] for i in xrange(2)], wavelength, inc, PA, np.array(Jylims[0])*unit_multiplier, np.array(Jyticks_im)*unit_multiplier, linestyle, colors,
                uvlim=[0., uvupperlim[0]], ismodel=ismodel,
                fontsize=30, fontsize_ticklabels=20, unit_multiplier=unit_multiplier,
                apply_shift=apply_shift)

    ax = fig.get_axes()
    ax[0].text(uvupperlim[0]*0.85, ax[0].get_ylim()[1]*unit_multiplier*0.85, SOURCE_LABEL, fontsize=30, fontweight='bold', ha='right')
    for tick in ax[1].xaxis.get_major_ticks():
                tick.label.set_fontsize(25)
    for tick in ax[1].yaxis.get_major_ticks():
                tick.label.set_fontsize(25)
    for tick in ax[0].yaxis.get_major_ticks():
                tick.label.set_fontsize(25)

    fig.tight_layout()
    fig.savefig(plot_filename)





if do_2d_contours:

    flatchain = util_chain.chain_thin(chain.copy(), start_mcmc, stop_mcmc, step_mcmc, logpars)
    from corner import hist2d
    import matplotlib as mpl

    import matplotlib.pyplot as plt
    # print(mpl.matplotlib_fname())  # matplotlibrc
    mpl.rcParams['axes.linewidth'] = 1.7
    mpl.rcParams['xtick.major.size'] = 14
    mpl.rcParams['xtick.minor.size'] = 7
    mpl.rcParams['xtick.major.width'] = 1.7
    mpl.rcParams['xtick.minor.width'] = 1.2
    mpl.rcParams['xtick.major.pad'] = 8
    mpl.rcParams['xtick.labelsize'] = 'large'
    mpl.rcParams['ytick.major.size'] = 14
    mpl.rcParams['ytick.minor.size'] = 7
    mpl.rcParams['ytick.major.width'] = 1.7
    mpl.rcParams['ytick.minor.width'] = 1.2
    mpl.rcParams['ytick.major.pad'] = 8
    mpl.rcParams['ytick.labelsize'] = 'large'
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.weight'] = 500
    mpl.rcParams['font.size'] = 16
    # if we want to plot gamma-Rc and p-Rout in two panels
    # from matplotlib import gridspec
    # fig = plt.figure(figsize=(10, 5))
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    # ax = plt.subplot(gs[0]), plt.subplot(gs[1])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    levels = 1.-np.exp(-0.5*np.arange(1., 3.1, 1)**2.)  # correspond to 1, 2 and 3 sigma for a 2D hystogram

    hist2d(flatchain[:, 0], flatchain[:, 2], ax=ax, bins=30, levels=levels, plot_datapoints=False,
            smooth=1.5, plot_contours=True, fill_contours=True, plot_density=False, color='#2D89E5')

    # add point for best-fit model
    ax.plot(bestfit_pars[0], bestfit_pars[2], marker='*',
     markerfacecolor='#FFFF66', markeredgecolor='#FFFF66', markersize=16)

    if headers['fit']['sigma_prescription'] == 'g':
        ax.set_xlabel(r'$\gamma$', fontsize=20)
        ax.set_ylabel(r'$R_c$', fontsize=20)
    elif headers['fit']['sigma_prescription'] == 'pl':
        ax.set_xlabel(r'$p$', fontsize=20)
        ax.set_ylabel(r'$R_{\mathrm{out}}$', fontsize=20)

    # set axes limits
    ax.set_ylim(0, 300.)
    yticks = np.arange(0, 301, 100, dtype='int')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    if headers['fit']['sigma_prescription'] == 'g':
        xticks = np.arange(-2, 2.1, 1, dtype='int')
        ax.set_xlim(-2., 2.)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
    elif headers['fit']['sigma_prescription'] == 'pl':
        xticks = np.arange(-2, 4.1, 1, dtype='int')
        ax.set_xlim(-2., 4.)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)

    ax.minorticks_on()
    sourcelabel_xcoord = (xticks[-1]-xticks[0])*0.15 + xticks[0]
    sourcelabel_ycoord = yticks[-1]*0.85
    ax.text(sourcelabel_xcoord, sourcelabel_ycoord, SOURCE_LABEL, fontsize=20, fontweight='bold')

    fig.tight_layout()

    fig.savefig(plots_dir+"dist_2D_"+STARNAME+".pdf")


















