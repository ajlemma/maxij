#######
''' maxij_savealigned.py
by A Townsend
- check files
- read in pandas dataframe with shift data & filenames
- open each image
- shift image
- rebin image
- save new image in aligned folder
'''
#######

from maxijdefs import *
import pandas as pd
from scipy.ndimage.interpolation import shift
from astropy.io import fits
from multiprocessing import Pool
from functools import partial
import os
import gc

#######

def savealigned(night, n1 = 0, n2 = 'max', path='/media/amanda/demeter/maxi_j1820_070/'):

    loglist = []  # initialize log

    # housekeeping stuff
    msg = "Running 'maxij_savealigned.py' code on data in '" + night + "' folder..."
    loglist = addlog(msg, loglist)

    time0, msg = timestart()
    loglist = addlog(msg, loglist)

    pathnam = path+night+'/'
    scipathnam = pathnam + 'science/'  # folder with science images
    apathnam = pathnam + 'aligned/'  # folder for aligned images


    # check to make sure all necessary files are available before starting
    msg = "Checking files..."
    loglist = addlog(msg,loglist)

    if not os.path.isdir(path+night):
        msg = "No folder exists for this night. Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isdir(path+night+"/aligned"):
        msg = "No 'aligned' folder found, creating 'aligned' folder..."
        loglist = addlog(msg, loglist)
        os.mkdir(path+night+"/aligned")

    if not os.path.isfile(pathnam+"data_"+night+".pkl"):
        msg = "No database called '" + "data_"+night+".pkl" + "' exists. \n " \
            "Please run initdb() to initialize the database. \n " \
            "Exiting..."
        loglist = addlog(msg, loglist)
        return None

    # load pandas df
    msg = "Loading database for " + night + "..."
    loglist = addlog(msg, loglist)
    dataf = pd.read_pickle(pathnam+"data_"+night+".pkl")

    fileID = dataf['fileID']
    fnames = dataf['filename']
    msg = "Total science frames:  " + str(len(fnames))
    loglist = addlog(msg, loglist)

    ## define last image
    if n2 == 'max':
        n2 = len(fnames)
    else:
        n2 = int(n2)

    # shift & rebin images and save them
    p = Pool(6) #uses 6 cores
    msg = "Shifting images..."
    loglist = addlog(msg, loglist)
    p.map(partial(align_rebin, scipathnam = scipathnam,apathnam = apathnam,dataf = dataf), fileID[n1:n2],35)

    msg = "Aligned & rebinned images are saved at " + night + '/aligned with the suffix "_aligned"'
    loglist = addlog(msg, loglist)


    #get stop time and save screen output to log
    msg = timefinish(time0)
    loglist = addlog(msg, loglist)

    # write/append screen output to logfile
    msg = "Writing screen output to logfile..."
    loglist = addlog(msg, loglist)
    writelog(loglist,night,pathnam)

    gc.collect()
    return None


def align_rebin(fID,scipathnam,apathnam,dataf):
    # print filename
    filename = dataf.loc[fID,'filename']
    im0 = fits.getdata(scipathnam+filename)
    im1 = shift(im0, (dataf.loc[fID,"shift_y"], dataf.loc[fID,"shift_x"]), mode='constant', cval=0.0)
    im2 = rebin(im1, (im1.shape[0] / 4, im1.shape[1] / 4))
    write_to_fits(apathnam + filename.split('.')[0]+'_aligned.' + filename.split('.')[1], im2)  # write to aligned subdir, add suffix "_aligned"
    im0 = None
    im1 = None
    im2 = None
    gc.collect()
    return


if __name__ == "__main__":
    savealigned('test', n2 = 10, path = './')
