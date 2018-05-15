#######
''' maxij_imalign.py
by A Townsend
1. read in imshifts file created by maxij_imshifts
2. open each image
3. shift image
4. rebin image
5. save new image in aligned folder
'''
#######

from maxijdefs import *
import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.signal import fftconvolve
from astropy.io import fits
import time
from datetime import datetime
# from multiprocessing import Pool

#######

night = '2018-03-28'
## path to the directory with the MAXIJ native-resolution fits files:
pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'
scipathnam = pathnam + 'science/' #folder with science images
apathnam = pathnam + 'aligned/' #folder for aligned images
time0 = time.time()
print "start: " + str(datetime.now())

fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam
ajt_shifts = pd.read_csv(pathnam + 'imshifts_' + night + '.txt',
                         header=None,
                         names = ["fname", "y", "x", "ampcorr"],
                         index_col='fname',
                         delimiter=r"\s+")
## to compare to steve's shifts (tbd)
# sse_shifts = pd.read_csv('/home/amanda/Desktop/maxij_photanalysis_2018-05-05/shifts328.txt',
#                          header=None,
#                          names=["iter", "fname", "y", "x", "ampcorr"],
#                          index_col='fname',
#                          delimiter=r"\s+")
#
# plt.plot(ajt_shifts.loc[fnames]["y"],sse_shifts.loc[fnames]["y"],'.')
# plt.show()

for img in fnames:
    im0 = fits.getdata(scipathnam+fname)






# print ajt_shifts.loc["maxij_1s_0027_03-25-31_reduced.fits"]["x"]


time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " m)"