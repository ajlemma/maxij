#######
''' maxij_initdb.py
by A Townsend
- check files
- initialize database as pandas dataframe
- parse file IDs and times and adds them to the datbase
- save dataframe to pickle file to be read into later code
'''
#######

import pandas as pd
from maxijdefs import *
from multiprocessing import Pool
from functools import partial
import os
import shutil

#######

def initdb(night, path='/media/amanda/demeter/maxi_j1820_070/'):
    loglist = [] #initialize log

    # housekeeping stuff
    msg = "Running 'maxij_initdb.py' code on data in '" + night + "' folder..."
    loglist = addlog(msg,loglist)

    time0,msg = timestart()
    loglist = addlog(msg, loglist)

    pathnam = path + night + '/'
    scipathnam = pathnam + 'science/'  # folder for science images

    # check to make sure all necessary files are available before starting
    msg = "Checking files..."
    loglist = addlog(msg,loglist)

    if not os.path.isdir(path+night):
        msg = "No folder exists for this night. Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isfile(pathnam + 'timestamps_' + night + '.txt'):
        msg = "Missing os timestamps file. Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isdir(path+night+"/science"):
        msg = "No 'science' folder found, please rename folder containing science images. \n Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isdir(path+night+"/aligned"):
        msg = "No 'aligned' folder found, creating 'aligned' folder..."
        loglist = addlog(msg, loglist)
        os.mkdir(path+night+"/aligned")

    if os.path.isdir(path+night+"/science/aligned"):
        msg = "Reference image 'aligned' folder still exists in 'science' directory."
        loglist = addlog(msg, loglist)
        msg = "Deleting '/science/aligned' folder..."
        loglist = addlog(msg, loglist)
        shutil.rmtree(path+night+"/science/aligned")

    if not os.path.isfile(path+night+"/ref_stack.fits"):
        msg = "Warning, no reference image found. \n" \
              "Please save a stacked reference image as 'ref_stack.fits' in the night folder. \n " \
              "Continuing database initialization... \n"
        loglist = addlog(msg, loglist)

    if not os.path.isfile(path+night+"/ref_stars_xy.txt"):
        msg = "Warning, no reference star position file found. \n" \
              "Please save an AIJ measurement file as 'ref_stars_xy.txt' in the night folder. \n" \
              "Continuing database initialization... \n"
        loglist = addlog(msg, loglist)


    msg = "Initializing database for " + night + "..."
    loglist = addlog(msg, loglist)

    # read in list of sorted filenames in the science folder
    msg = "Getting list of filenames from " + night + "/science ..."
    loglist = addlog(msg, loglist)
    fnames = get_filelist_maxi(scipathnam)
    fnames = fnames

    fileID = [n.split('_')[2] for n in fnames]

    # set up a dictionary of column headers for the pandas dataframe
    maxidat = {'fileID': fileID,
               'filename': fnames,

               # time values
               'filetime': np.nan,
               'filetime_s': np.nan,
               'os_time': np.nan,

               # image shifts & correlation amplitude => flag bad images
               'shift_x': np.nan,
               'shift_y': np.nan,
               'shift_corr_amplitude': np.nan,

               # guassian parameters (~fwhm) => aperture size, flagging bad data
               'gauss_params_0_amplitude': np.nan,
               'gauss_params_1_x0': np.nan,
               'gauss_params_2_y0': np.nan,
               'gauss_params_3_sigma_x': np.nan,
               'gauss_params_4_sigma_y': np.nan,
               'gauss_params_5_theta': np.nan,
               'gauss_params_6_offset': np.nan,

               'gauss_sigma_avg': np.nan,  # avg of sigma_x and sigma_y
               'gauss_offset_dx': np.nan,  # actual delta x offset of tyc star from initial meas.
               'gauss_offset_dy': np.nan,  # actual delta y offset of tyc star from initial meas.

               # maxij + refstar photometry
               'phot_tyc': np.nan,
               'phot_maxij': np.nan,
               'phot_ref2': np.nan,
               'phot_ref3': np.nan,
               'phot_ref4': np.nan,
               'phot_ref5': np.nan,
               'phot_ref6': np.nan,

               # sky photometry
               'sky_tyc': np.nan,
               'sky_maxij': np.nan,
               'sky_ref2': np.nan,
               'sky_ref3': np.nan,
               'sky_ref4': np.nan,
               'sky_ref5': np.nan,
               'sky_ref6': np.nan

               # add flags here ?
               }
    maxiframe = pd.DataFrame(maxidat, index=fileID) #create the df

    # read in os times data
    msg = "Getting os timestamp data..."
    loglist = addlog(msg, loglist)
    ostime = pd.read_csv(pathnam + 'timestamps_' + night + '.txt',
                         skiprows=2,
                         header=None,
                         names=['filename1', 'filename2', 'epochtime(s)'],
                         delimiter=r"\s+")

    ostime = ostime[ostime['filename1'].str.contains('maxij_1s')]
    fileID_os = [n.split('_')[-1] for n in ostime['filename1']]
    ostime['fileID'] = fileID_os
    ostime = ostime.set_index('fileID')

    # add times to the pandas df
    msg = "Adding times to database..."
    loglist = addlog(msg, loglist)
    p = Pool(6)
    for results in p.imap_unordered(partial(assign_times, maxiframe=maxiframe, ostime=ostime), fileID,25):
        # print results
        maxiframe.loc[results[0], 'filetime'] = results[1]
        maxiframe.loc[results[0], 'filetime_s'] = results[2]
        maxiframe.loc[results[0], 'os_time'] = results[3]

    # write pandas df to file in night folder
    msg = "Saving database to " + pathnam + 'data_'+night+'.pkl'
    loglist = addlog(msg, loglist)
    maxiframe.to_pickle(pathnam + 'data_'+night+'.pkl')

    #get stop time and save screen output to log
    msg = timefinish(time0)
    loglist = addlog(msg, loglist)

    # write/append screen output to logfile
    msg = "Writing screen output to logfile..."
    loglist = addlog(msg, loglist)
    writelog(loglist,night,pathnam)

    return maxiframe


def assign_times(id, maxiframe, ostime):
    ftime = parse_time(maxiframe.loc[id, 'filename'])
    ftime_s = get_sec(ftime)
    ost = ostime.loc[id, 'epochtime(s)']
    return id, ftime, ftime_s, ost


if __name__ == "__main__":
    initdb('test', path='./')
    # print initdb('2018-03-28')
    # initdb('2018-04-13')
    # initdb('brokentest')