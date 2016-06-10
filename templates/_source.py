from numpy import array, argsort, min

# starname, identifier of the fit
STARNAME = '<starname>'

# label used in plots
SOURCE_LABEL = '<source_label>'

# threshold for the acceptance ratio.
# If None, the worse 50% of walkers are automatically excluded
ACC_RATIO_THRS = None

# numer of steps to consider for the MCMC analysis (counting from the end)
NUM_STEPS = 3000

# number of models to compute/plot (among the models selected in NUM_STEPS)
NUM_MODELS = 3000

# parameters of the bestfit model
# If None, the parameters are determined automatically using get_bestfit_model()
BESTFIT_PARS = None

# Filename of the origina MS table (usually the one containing the observations)
ORIGINAL_MS_FILENAME = '<starname>.ms'


# *** PLOTS ***

# * MCMC plot *
# range of the MCMC triangle plot.
# If None, it is set automatically
mcmc_range = None

# * UV plots *
# For all the arrays the first dimension has length nwle
# x-axis limits, for each wle: [min(uv), max(uv)]
UVLIM = [[0., 800.]]

# size of the uv bins. Used for plotting, but also when computing visibilities
UVBINSIZE = [30.e3]

# y-axis limits, for each wle: [(Re_lower, Re_upper), (Im_lower, Im_upper)]
Jylims = [[(-0.01, 0.1), (-0.01, 0.01)]]

# y-axis ticks, for each wle: [[Re_ticks], [Im_ticks]]
Jyticks = [[[0., 0.02, 0.04, 0.06, 0.08], [-0.005, 0., 0.005]]]

# y-axis ticklabels, for each wle: [[Re_ticklabels], [Im_ticklabels]]
Jyticklabels = Jyticks

# y-axis units, for each wle: [[Re_unit], [Im_unit]]
Jyunit = [['(Jy)'], ['(Jy)']]


def get_bestfit_model(flatchain, blobs):
    # sort by increasing chisquare
    idx = argsort(blobs.flatten())
    flatchain = flatchain[idx, :]

    # uncomment to define criteria
    # criteria = (flatchain[:, 1] < 2.) # example: & (flatchain[:, 0] < 1.) & (...)

    try:
        # apply the criteria, if given
        flatchain = flatchain[criteria]
        print("Bestfit model chi-square: {0}".format(blobs.flatten()[criteria][0]))
    except NameError:
        # criteria is not defined, leave flatchain as it is.
        pass

    print("Minimum chi-square in the chain: {0}".format(min(blobs)))

    return flatchain[0]

# TODO: what to do with nwle = 1? For the moment it is not used.
