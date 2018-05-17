#######
''' maxij_database_init.py
by A Townsend
1.
'''
#######

import pandas as pd
from maxijdefs import *
import matplotlib.pyplot as plt
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


pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'
apathnam = pathnam + 'aligned/' #folder for aligned images

fnames = get_filelist_maxi(apathnam)  # list of sorted filenames for aligned images in pathnam
fnames = fnames
ostime = pd.read_csv(pathnam + 'os_timestamps.txt',
                     skiprows=2,
                     header=None,
                     names=['filename1','filename2','epochtime(s)'],
                     delimiter=r"\s+")

ostime = ostime[ostime['filename1'].str.contains('maxij_1s')]
fileID_os = [n.split('_')[-1] for n in ostime['filename1']]
ostime['fileID'] = fileID_os
ostime = ostime.set_index('fileID')

# fileID_os = natsorted(fileID_os)


## set up a pandas dataframe for collection of data to output
fileID = [n.split('_')[2] for n in fnames]
maxidat = {'fileID': fileID,
           'filename': fnames,

           'filetime': np.nan,
           'filetime_s': np.nan,
           'os_time': np.nan,

           # 'shift_x': np.nan,
           # 'shift_y': np.nan,
           # 'shift_corr_amplitude': np.nan,
           #
           # 'gauss_params_0_amplitude': np.nan,
           # 'gauss_params_1_x0': np.nan,
           # 'gauss_params_2_y0': np.nan,
           # 'gauss_params_3_sigma_x': np.nan,
           # 'gauss_params_4_sigma_y': np.nan,
           # 'gauss_params_5_theta': np.nan,
           # 'gauss_params_6_offset': np.nan,
           # 'gauss_sigma_avg': np.nan,  # avg of sigma_x and sigma_y
           #
           # 'phot_tyc': np.nan,
           # 'phot_maxij': np.nan,
           # 'phot_ref2': np.nan,
           # 'phot_ref3': np.nan,
           # 'phot_ref4': np.nan,
           # 'phot_ref6': np.nan
           }
maxiframe = pd.DataFrame(maxidat, index=fileID)

#parallelize this
def assign_times(id):
    maxiframe.loc[id,'filetime'] = parse_time(maxiframe.loc[id,'filename'])
    maxiframe.loc[id,'filetime_s'] = get_sec(maxiframe.loc[id,'filetime'])
    maxiframe.loc[id,'os_time'] = ostime.loc[id,'epochtime(s)']

p = Pool(6)
p.map(assign_times,fnames,25)
    # print ostime.loc[id,'epochtime(s)']
    # print maxiframe.loc[id,'os_time']
# print maxiframe

# plt.plot(maxiframe['fileID'],maxiframe['os_time'],'r.')
# plt.plot(maxiframe['fileID'],maxiframe['filetime_s'],'b.')
# plt.show()