#######
''' maxij_multiphot.py
by A Townsend
1. this code was never finished--now is maxij_initdb and maxij_dophot
'''
#######

from maxijdefs import *
import pandas as pd
from natsort import natsorted
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
ostime = pd.read_csv(pathnam + 'timestamps_' + night + '.txt',
                     skiprows=2,
                     header=None,
                     names=['filename1','filename2','epochtime(s)'],
                     delimiter=r"\s+")
# print len(ostime)
ostime = ostime[ostime['filename1'].str.contains('maxij_1s')]
fileID_os = [n.split('_')[-1] for n in ostime['filename1']]
ostime['fileID'] = fileID_os
ostime = ostime.set_index('fileID')

fileID_os = natsorted(fileID_os)

## shifted database to it's own code (maxi_database_init) so read that in here instead & then do photometry
# print "initializing database..."
# ## set up a pandas dataframe for collection of data to output
# fileID = [n.split('_')[2] for n in fnames]
# maxidat = {'fileID': fileID,
#            'filename': fnames,
#
#            'filetime': np.nan,
#            'filetime_s': np.nan,
#            'os_time': np.nan,
#
#            'shift_x': np.nan,
#            'shift_y': np.nan,
#            'shift_corr_amplitude': np.nan,
#
#            'gauss_params_0_amplitude': np.nan,
#            'gauss_params_1_x0': np.nan,
#            'gauss_params_2_y0': np.nan,
#            'gauss_params_3_sigma_x': np.nan,
#            'gauss_params_4_sigma_y': np.nan,
#            'gauss_params_5_theta': np.nan,
#            'gauss_params_6_offset': np.nan,
#            'gauss_sigma_avg': np.nan,   #avg of sigma_x and sigma_y
#
#            'phot_tyc': np.nan,
#            'phot_maxij': np.nan,
#            'phot_ref2': np.nan,
#            'phot_ref3': np.nan,
#            'phot_ref4': np.nan,
#            'phot_ref6': np.nan
#            }
# maxiframe = pd.DataFrame(maxidat,index=fileID)
#
#
# print maxiframe

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
