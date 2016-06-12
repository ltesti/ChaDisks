import sys
sys.path.append('/home/ltesti/git_repository/')
from ltutils.uvdata import UVDataMS
from astropy.io import ascii as aio

# if `True run in simulation
debug = True

# Datafiles and sources
datadir = '../data_fits/'

tab = aio.read('table_wizard.txt')

datafiles_orig = []
datafiles_now0 = []
datafiles = []
outfiles = []
phasec_string = []
for src in tab:
    datafiles_orig.append(datadir+src['Name']+'_selfcal.ms')
    datafiles_now0.append(datadir+src['Name']+'_selfcal_now0.ms')
    datafiles.append(datadir+src['Name']+'_selfcal_now0_phc.ms')
    outfiles.append('../uvdata/alma_b7_'+src['ID']+'_now0_phc.txt')
    phasec_string.append('J2000 '+src['ra_c']+' '+src['dec_c'])


for i in range(len(datafiles)):
    #
    # split out removing w=0 data
    print("--> Splitting ms table: {0} into {1} with keepflags=False".format(datafiles_orig[i],datafiles_now0[i]))
    if debug:
        print('    split(vis={0},outputvis={1},datacolumn='data',keepflags=False)'.format(datafiles_orig[i],datafiles_now0[i]))
    else:
        split(vis=datafiles_orig[i],outputvis=datafiles_now0[i],datacolumn='data',keepflags=False)
    print("--> Adjusting phase center of: {0} into {1}".format(datafiles_now0[i],datafiles[i]))
    if debug:
        print('    fixvis(vis={0},outputvis={1},phasecenter={2})'.format(datafiles_now0[i],datafiles[i],phasec_string[i]))
    else:
        fixvis(vis=datafiles_now0[i],outputvis=datafiles[i],phasecenter=phasec_string[i])
    print("--> Importing ms table: {0}".format(datafiles[i]))
    if debug:
        print('    myuvdata = UVDataMS("dummy", ({0}, tb))'.format(datafiles[i]))
        print('    rat_re,rat_im = myuvdata.get_weight()')
        print('    myuvdata.we = myuvdata.we*(rat_re+rat_im)/2.')
    else:
        myuvdata = UVDataMS("dummy", (datafiles[i], tb))
        rat_re,rat_im = myuvdata.get_weight()
        myuvdata.we = myuvdata.we*(rat_re+rat_im)/2.
    print("<-- Exporting uv table: {0}".format(outfiles[i]))
    if debug:
        print('     myuvdata.write_uv_to_ascii({0})'.format(outfiles[i]))
    else:
        myuvdata.write_uv_to_ascii(outfiles[i])   


