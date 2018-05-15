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
import numpy as np
# import matplotlib.pyplot as plt
# from scipy.signal import fftconvolve
from scipy.ndimage.interpolation import shift
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
imshifts = pd.read_csv(pathnam + 'imshifts_' + night + '.txt',
                         header=None,
                         names = ["fname", "y", "x", "ampcorr"],
                         index_col='fname',
                         delimiter=r"\s+")

# for img in fnames:

i = 0
im0 = fits.getdata(scipathnam+fnames[i])

im1 = shift(im0, (imshifts.loc[fnames[i]]["y"], imshifts.loc[fnames[i]]["x"]), mode='constant', cval=0.0)
im2 = rebin(im1, (im1.shape[0] / 4, im1.shape[1] / 4))
write_to_fits(apathnam + fnames[i].split('.')[0]+'_aligned.' + fnames[i].split('.')[1], im2)  # write to aligned subdir, add suffix "_aligned"

# next: turn above "for loop" into multiprocessing fns


# print ajt_shifts.loc["maxij_1s_0027_03-25-31_reduced.fits"]["x"] ## how to use indexes in pd df


time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " m)"