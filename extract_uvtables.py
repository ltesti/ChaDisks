import sys
sys.path.append('/home/ltesti/git_repository/')
from ltutils.uvdata import UVDataMS
from astropy.io import ascii as aio

# Datafiles and sources
datadir = '../data_fits/'

tab = aio.read('table_wizard.txt')

datafiles = []
outfiles =[]
for src in tab:
    datafiles.append(datadir+src['Name']+'_selfcal.ms')
    outfiles.append('../uvdata/alma_b7_'+src['ID']+'.txt')


for i in range(len(datafiles)):
	print("--> Importing ms table: {0}".format(datafiles[i]))
	myuvdata = UVDataMS("dummy", (datafiles[i], tb))
	rat_re,rat_im = myuvdata.get_weight()
	myuvdata.we = myuvdata.we*(rat_re+rat_im)/2.
	print("<-- Exporting uv table: {0}".format(outfiles[i]))
	myuvdata.write_uv_to_ascii(outfiles[i])   


