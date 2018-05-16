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
# scipathnam = pathnam + 'science/' #folder with science images
apathnam = pathnam + 'aligned/' #folder for aligned images

time0 = time.time()
print "start: " + str(datetime.now())
print "Processing images for " + night

print "getting list of filenames from " + night + "/science ..."
fnames = get_filelist_maxi(apathnam)  # list of sorted filenames for aligned images in pathnam

print "getting list of reference star positions..."
refstars = pd.read_csv(pathnam + 'ref_stars_xy.txt', delimiter=r"\s+")

print "getting os timestamps..."
ostime = pd.read_csv(pathnam + 'os_timestamps.txt',
                     skiprows=2,
                     header=None,
                     names=['filename1','filename2','epochtime(s)'],
                     delimiter=r"\s+")
# print len(ostime)
ostime = ostime[ostime['filename1'].str.contains('maxij_1s')]
# print ostime['filename1'][0],ostime['filename2'][0]
for n in ostime['filename1']:
    print n.split('_')[-1]


# ## first and last images to photometer, and total no. of images
# n1 = 0                          # first
# n2 = 1 #len(fnames)                # last
# n_tot = n2 - n1                 # total

# for filename in fnames[n1:n2]:
#     im0 = fits.getdata(apathnam + filename)
#     tsec = parse_time(filename)
#     print filename
#     print tsec
#     # print ostime.loc[filename]

time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " m)"
