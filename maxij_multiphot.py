#######
''' maxij_multiphot.py
by A Townsend
1.
'''
#######

from maxijdefs import *
import pandas as pd
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

print "getting list of reference star positions..."
refstars = pd.read_csv(pathnam + 'ref_stars_xy.txt', delimiter=r"\s+")

print refstars['X(FITS)'][0]

time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " m)"
