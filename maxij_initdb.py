#######
''' maxij_initdb.py
by A Townsend
- checks to make sure all files are available in night folder
- initializes database as pandas dataframe
- parses file IDs and times and adds them to the datbase
- saves to pickle file to be read into later code
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
    time0 = timestart()
    pathnam = path + night + '/'
    spathnam = pathnam + 'science/'  # folder for science images


    print "Checking files..."

    if not os.path.isdir(path+night):
        print "No folder exists for this night. Exiting..."
        return None

    if not os.path.isfile(pathnam + 'timestamps_' + night + '.txt'):
        print "Missing os timestamps file. Exiting..."
        return None

    if not os.path.isdir(path+night+"/science"):
        print "No 'science' folder found, please rename folder containing science images."
        print "Exiting..."
        return None

    if not os.path.isdir(path+night+"/aligned"):
        print "No 'aligned' folder found, creating 'aligned' folder..."
        os.mkdir(path+night+"/aligned")

    if os.path.isdir(path+night+"/science/aligned"):
        print "Reference image 'aligned' folder still exists in 'science' directory."
        print "Deleting '/science/aligned' folder..."
        shutil.rmtree(path+night+"/science/aligned")

    if not os.path.isfile(path+night+"/ref_stack.fits"):
        print "Warning, no reference image found."
        print "Please save a stacked reference image as 'ref_stack.fits' in the night folder."
        print "Continuing database initialization..."
        print

    if not os.path.isfile(path+night+"/ref_stars_xy.txt"):
        print "Warning, no reference star position file found."
        print "Please save an AIJ measurement file as 'ref_stars_xy.txt' in the night folder."
        print "Continuing database initialization..."
        print


    print "Initializing database for " + night + "..."

    print "Getting list of filenames from " + night + "/science ..."
    fnames = get_filelist_maxi(spathnam)  # list of sorted filenames for aligned images in pathnam
    fnames = fnames

    fileID = [n.split('_')[2] for n in fnames]
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

               # add flags here
               }
    maxiframe = pd.DataFrame(maxidat, index=fileID)

    print "Getting os timestamp data..."
    ostime = pd.read_csv(pathnam + 'timestamps_' + night + '.txt',
                         skiprows=2,
                         header=None,
                         names=['filename1', 'filename2', 'epochtime(s)'],
                         delimiter=r"\s+")

    ostime = ostime[ostime['filename1'].str.contains('maxij_1s')]
    fileID_os = [n.split('_')[-1] for n in ostime['filename1']]
    ostime['fileID'] = fileID_os
    ostime = ostime.set_index('fileID')

    print "Adding times to database..."
    p = Pool(6)
    for results in p.imap_unordered(partial(assign_times, maxiframe=maxiframe, ostime=ostime), fileID,25):
        # print results
        maxiframe.loc[results[0], 'filetime'] = results[1]
        maxiframe.loc[results[0], 'filetime_s'] = results[2]
        maxiframe.loc[results[0], 'os_time'] = results[3]

    print "Saving database to " + pathnam + 'data_'+night+'.pkl'
    maxiframe.to_pickle(pathnam + 'data_'+night+'.pkl')

    timefinish(time0)
    return maxiframe


def assign_times(id, maxiframe, ostime):
    ftime = parse_time(maxiframe.loc[id, 'filename'])
    ftime_s = get_sec(ftime)
    ost = ostime.loc[id, 'epochtime(s)']
    return id, ftime, ftime_s, ost


if __name__ == "__main__":
    # initdb('test', path='./')
    # print initdb('2018-03-28')
    initdb('2018-04-13')
    # initdb('brokentest')