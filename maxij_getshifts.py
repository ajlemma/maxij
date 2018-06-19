#######
''' maxij_getshifts.py
Original Code by S Eikenberry, modified by A Townsend
- check files
- read in pandas dataframe
- read in reference image
- cross correlate ref img with itself to get a reference baseline
- cross correlate ref img with each individual science image for the night
---> uses multiple processors
- add image shifts to dataframe & resave it
'''
#######

from maxijdefs import *
from scipy.signal import fftconvolve
from astropy.io import fits
from multiprocessing import Pool
from functools import partial
import pandas as pd
import os
import shutil
from datetime import datetime

#######

def getshifts(night, n1 = 0, n2 = 'max',path='/media/amanda/demeter/maxi_j1820_070/'):
    # n2 can also be an integer

    loglist = []  # initialize log

    # housekeeping stuff
    msg = "Running 'maxij_getshifts.py' code on data in '" + night + "' folder..."
    loglist = addlog(msg, loglist)

    time0, msg = timestart()
    loglist = addlog(msg, loglist)

    pathnam = path+night+'/'
    scipathnam = pathnam + 'science/'  # folder with science images

    # check to make sure all necessary files are available before starting
    msg = "Checking files..."
    loglist = addlog(msg,loglist)

    if not os.path.isdir(path+night):
        msg = "No folder exists for this night. Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isfile(pathnam+"data_"+night+".pkl"):
        msg = "No database called '" + "data_"+night+".pkl" + "' exists. \n " \
            "Please run code to initialize the database. \n " \
            "Exiting..."
        loglist = addlog(msg, loglist)
        return None

    fdate = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    msg = "Making backup of database as 'data_" + night + "_bck_"+fdate+".pkl"
    loglist = addlog(msg, loglist)
    shutil.copyfile(pathnam+"data_"+night+".pkl",pathnam+"data_"+night+"_bck_"+fdate+".pkl")

    if not os.path.isfile(path+night+"/ref_stack.fits"):
        msg = "No reference image found. \n" \
              "Please save a stacked reference image as 'ref_stack.fits' in the night folder. \n" \
              "Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isfile(path+night+"/ref_stars_xy.txt"):
        msg = "Warning, no reference star position file found. \n" \
              "Please save an AIJ measurement file as 'ref_stars_xy.txt' in the night folder. \n" \
              "Continuing shift calculations..."
        loglist = addlog(msg, loglist)


    msg = "Loading database for " + night + "..."
    loglist = addlog(msg, loglist)
    dataf = pd.read_pickle(pathnam + 'data_'+night+'.pkl')
    fileID = dataf['fileID']

    fnames = dataf['filename']
    msg = "Total science frames:  " + str(len(fnames))
    loglist = addlog(msg, loglist)

    ## get rid of low-level wave information for reference image & cross-correlate it with itself:
    msg = "Getting correlation reference..."
    loglist = addlog(msg, loglist)
    imref = fits.getdata(pathnam+'ref_stack.fits') #read in stacked image
    imref[imref < 1e3] = 0
    corr_ref = fftconvolve(imref, imref[::-1, ::-1])  # cross-correlate them

    ## define last image
    # to stack
    if n2 == 'max':
        n2 = len(fnames)
    else:
        n2 = int(n2)


    p = Pool(6)
    msg = "calculating individual image shifts..."
    loglist = addlog(msg, loglist)
    for results in p.imap_unordered(partial(get_shifts, imref = imref,corr_ref = corr_ref,scipathnam = scipathnam,dataf=dataf),
                                  fileID[n1:n2], 25):
        fid = results[0]
        dataf.loc[fid, 'shift_y'] = results[1]
        dataf.loc[fid, 'shift_x'] = results[2]
        dataf.loc[fid, 'shift_corr_amplitude'] = results[3]

    msg = "Image shifts are saved in data_"+night+".pkl"
    loglist = addlog(msg, loglist)


    #get stop time and save screen output to log
    msg = timefinish(time0)
    loglist = addlog(msg, loglist)

    # write/append screen output to logfile
    msg = "Writing screen output to logfile..."
    loglist = addlog(msg, loglist)
    writelog(loglist,night,pathnam)


def get_shifts(fID,imref,corr_ref,scipathnam,dataf):
    # finds the shift of science_image with respect to reference_image
    # uses pre-calculated correlation_reference baseline (so you don't have to recalculate it every time)
    science_image = dataf.loc[fID,'filename']
    reference_image = imref
    correlation_reference = corr_ref
    fnam = scipathnam + science_image  # make the full filename
    im0 = fits.getdata(fnam)

    im0 = im0 - np.median(im0)
    im0[im0 < 1e3] = 0

    correlation1 = fftconvolve(im0, reference_image[::-1, ::-1])  # cross-correlate them
    correlation2 = correlation1 - correlation_reference

    xy = zip(*np.where(correlation2 == np.max(correlation2)))  # get index of max)
    xc = xy[0][1]
    yc = xy[0][0]

    # print xc, yc
    return (fID, (correlation2.shape[0] / 2) - int(yc), (correlation2.shape[1] / 2) - int(xc), np.max(correlation2) / 1e10)  # required pixel shifts


if __name__ == "__main__":
    getshifts('test', n2 = 10, path = './')
