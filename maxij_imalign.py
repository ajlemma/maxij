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
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.signal import fftconvolve
from scipy.ndimage.interpolation import shift
from astropy.io import fits
import time
from datetime import datetime
from multiprocessing import Pool

#######

night = '2018-03-28'
## path to the directory with the MAXIJ native-resolution fits files:
pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'
scipathnam = pathnam + 'science/' #folder with science images
apathnam = pathnam + 'aligned/' #folder for aligned images
time0 = time.time()
print "start: " + str(datetime.now())
print "Processing images for " + night

print "getting list of filenames from " + night + "/science ..."
fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam

print "getting image shift data from " + night + '/imshifts_' + night + ".txt ..."
imshifts = pd.read_csv(pathnam + 'imshifts_' + night + '.txt',
                         header=None,
                         names = ["fname", "y", "x", "ampcorr"],
                         index_col='fname',
                         delimiter=r"\s+")


## first and last images to stack, and total no. of images
n1 = 0                          # first
n2 = len(fnames)                # last
n_tot = n2 - n1                 # total

## Parallelized version:
def align_rebin(filename):
    # i = 0
    # print filename
    im0 = fits.getdata(scipathnam+filename)
    im1 = shift(im0, (imshifts.loc[filename]["y"], imshifts.loc[filename]["x"]), mode='constant', cval=0.0)
    im2 = rebin(im1, (im1.shape[0] / 4, im1.shape[1] / 4))
    write_to_fits(apathnam + filename.split('.')[0]+'_aligned.' + filename.split('.')[1], im2)  # write to aligned subdir, add suffix "_aligned"
    return

p = Pool(7)
print "shifting images..."
# for img in fnames[n1:n2]:
#     align_rebin(img)
# p.imap_unordered(align_rebin,fnames[n1:n2],10)

p.map(align_rebin,fnames[n1:n2],10)

print "Aligned & rebinned images are saved at " + night + '/aligned with the suffix "_aligned"'


time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " m)"