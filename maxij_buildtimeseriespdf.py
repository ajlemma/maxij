#######
''' maxij_buildtimeseriespdf.py
by A Townsend
'''
#######

import pandas as pd
from maxijdefs import *
from maxij_nicerdefs import *

#######

time0,msg = timestart()
print msg

#######

## time calibration stuff:
# computer clock is ~51.5s fast (per lab tests + SSE calculation)
# (((POSSIBLY 8.5 s slow, but this seems to work?)))
clocklag = -51.5
# clocklag = 8.5
# clocklag -=60

# do not change this! it is absolutely calibrated to the 28th.
# unix epoch timestamp - nicer absolute time
# both at 2018-03-28 08:05:21.000 UTC
dt = 1522224321-133689923
# print dt

#######
# info to read rho optical data (directory and nights)

opathnam = './maxij_data_copies/'
onights = ['2018-03-28','2018-03-29','2018-04-06','2018-04-12',
           '2018-04-13','2018-04-17','2018-04-18',
           '2018-04-25','2018-04-26','2018-04-30',
           '2018-05-01','2018-05-02','2018-05-03','2018-05-04',
           '2018-05-05','2018-05-07','2018-05-09']

#######
print
#read in x-ray data as pandas dataframe:
print "reading in all x-ray data"
nicerdata = pd.read_pickle('./' + 'xray_data_all_silverGTI.pkl')
# print data
print "there are " + str(len(nicerdata)) + " total x-ray datapoints"

print "converting nicer time values to unix epoch time"
nicerdata['unix_time'] = nicerdata['time'] + dt

## pull out just datasets we are going to use
xt = nicerdata['unix_time']
xlc = nicerdata['cts']

#######
print
print "setting up a pandas dataframe for the timeseries data"
cols = ['nights' , 'unix_epoch_time' , 'xray_maxij_ts' , 'rho_maxij_ts' , 'rho_tycho_ts' ]
tseries_df = pd.DataFrame(columns = cols, index = onights)
tseries_df['nights'] = onights
# print tseries_df
#######
print
print "reading in optical data..."
print

for night in onights:
    #read in optical data as pandas dataframe:
    print "reading in optical data for " + night
    rhodata = pd.read_pickle(opathnam+'data_'+night+'.pkl')
    # print data
    print "there are " + str(len(rhodata)) + " datapoints on " + night

    print
    print "calibrating time values with offset to unix epoch time"
    rhodata['fixed_time'] = rhodata['os_time'] + clocklag

    print
    print "cleaning up data"

    ## clean up data
    shiftgood = rhodata['shift_corr_amplitude'] < 10.
    gaussgood = rhodata['gauss_sigma_avg'] < 3.
    gaussnonzero = rhodata['gauss_sigma_avg'] > 0.
    gaussgood2 = np.abs(rhodata['gauss_params_3_sigma_x'] - rhodata['gauss_params_4_sigma_y']) < 1.
    rhoclean = rhodata[shiftgood & gaussgood & gaussgood2 & gaussnonzero]


    print
    print "pulling out time and photometry datasets and creating 1-s timeseries"

    ## pull out just datasets we are going to use

    ot = rhoclean['fixed_time']
    olc = rhoclean['phot_maxij']
    tlc = rhoclean['phot_tyc']

    #######

    # ## define timeframe based on optical data
    # ## make 1-s timeseries' for whole night (with zeroes where no data)
    optical_time_seconds = np.arange(np.floor(np.min(rhoclean['fixed_time'])),np.ceil(np.max(rhoclean['fixed_time'])))
    tseries_df.at[night,'unix_epoch_time'] = optical_time_seconds

    xts = make_1s_ts(optical_time_seconds, np.array(xt),np.array(xlc))
    tseries_df.at[night,'xray_maxij_ts'] = xts

    ots = make_1s_ts(optical_time_seconds, np.array(ot),np.array(olc))
    tseries_df.at[night, 'rho_maxij_ts'] = ots

    tts = make_1s_ts(optical_time_seconds, np.array(ot),np.array(tlc))
    tseries_df.at[night, 'rho_tycho_ts'] = tts

# print tseries_df['unix_epoch_time']

#######
print
# write pandas df to file in night folder
print "Saving database to ./timeseries_data_all.pkl"
tseries_df.to_pickle('./timeseries_data_all.pkl')

#######
print
print timefinish(time0)