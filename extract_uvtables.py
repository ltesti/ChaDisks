import sys
sys.path.append('/home/mtazzari/repos/')
from ltutils.uvdata import UVDataMS

# Datafiles and sources
datafiles = ['G/science_calibrated.ms',
            'M/science_calibrated.ms']
sources = [['9', '14', '28'],
          ['10']]
newnames_root = ['G', 'M']

# defines a list with all new names
newnames = []
for dat in range(len(datafiles)):
       newnames.append([])
       for ns in range(len(sources[dat])):
           newnames[dat].append(newnames_root[dat]+sources[dat][ns])

for dat in range(len(datafiles)):
	for i, filename in enumerate(newnames[dat]):
		print("--> Importing ms table: {0}".format(filename+".ms"))
		myuvdata = UVDataMS("dummy", (filename+".ms", tb))
		rat_re,rat_im = myuvdata.get_weight()
		myuvdata.we = myuvdata.we*(rat_re+rat_im)/2.
		print("<-- Exporting uv table: {0}".format("alma_b7_"+filename+".txt"))
		myuvdata.write_uv_to_ascii("alma_b7_"+filename+".txt”)


