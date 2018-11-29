#######
''' maxij_overlap_analysis_28nov2018.py
by A Townsend
'''
#######

import pandas as pd
# import numpy as np
from maxijdefs import *
from maxij_nicerdefs import *
import glob
# from natsort import natsorted
# from itertools import groupby
# from datetime import datetime
# from scipy.signal import fftconvolve
# import matplotlib
# import matplotlib.pyplot as plt
# from pylab import subplot, subplots_adjust, legend

#######

time0,msg = timestart()
print msg

#######

# info for reading x-ray data (directory and filter for night)
# for example, for 3/28 the NICER date is 113
# filename_filter = 'ni1200120113_0mpu7_silver_GTI*.lc' # for 3/28
xpathnam = 'delivery/'
filename_filter = '*_0mpu7_silver_GTI*.lc'

#######

## read in nicer xray data for night and make pandas df
# get filenames for night's xray data
print "reading nicer x-ray data"
ftest = xpathnam + filename_filter
flistD = glob.glob(ftest)

# all files
flist = [fnam.split('/')[1] for fnam in flistD]
flist = natsorted(flist)
# print len(flist)
ndat_t = []
ndat_cts = []
ndat_cts1 = []
ndat_cts2 = []
ndat_cts34 = []
ndat_fnam = []

for fileX in flist:
    t,cts=rd_nicer_lc(xpathnam,[fileX])
    t1,cts1=rd_nicer_lc1(xpathnam,[fileX])
    t2,cts2=rd_nicer_lc2(xpathnam,[fileX])
    t34,cts34=rd_nicer_lc34(xpathnam,[fileX])
    filename_column = [fileX for x in np.arange(len(t))]

    ndat_t.append(t)
    ndat_cts.append(cts)
    ndat_cts1.append(cts1)
    ndat_cts2.append(cts2)
    ndat_cts34.append(cts34)
    ndat_fnam.append(filename_column)

ndict = {'filename': [item for sublist in ndat_fnam for item in sublist],
          'time': [item for sublist in ndat_t for item in sublist],
          'cts': [item for sublist in ndat_cts for item in sublist],
          'cts1': [item for sublist in ndat_cts1 for item in sublist],
          'cts2': [item for sublist in ndat_cts2 for item in sublist],
          'cts34': [item for sublist in ndat_cts34 for item in sublist]
          }
print "creating data frame for nicer x-ray data"
nicerdata = pd.DataFrame(ndict, index=ndict['time'])
# print nicerdata

# write pandas df to file in night folder
print "Saving database as xray_data_all_silverGTI.pkl"

nicerdata.to_pickle('./' + 'xray_data_all_silverGTI.pkl')

print timefinish(time0)