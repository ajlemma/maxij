#######
''' maxij_dophot.py
by A Townsend, based on/using original code by SSE
- check files
- read in pandas dataframe with filenames
- read in list of reference star positions
- first star in list is the bright tycho reference star; (2nd is maxij1820+070)
- adjust aperture positions to bright star position (minor adjustment)
- & choose aperture size based on bright star gaussian fit
- photometer each reference star (aperture photometry with a sky annulus)
- add photometry & all parameters to the df and save it
- append/save screen output to log file

takes about 2.5 min for 10,000 images
'''
#######

import pandas as pd
import numpy as np
from maxijdefs import *
from multiprocessing import Pool
from functools import partial
from astropy.io import fits
import os
from datetime import datetime
import shutil
import gc
#######

## note: if dophot throws lots of errors and fails, you probably need a better reference image!!!

def dophot(night, path='/media/amanda/demeter/maxi_j1820_070/'):

    loglist = []  # initialize log

    # housekeeping stuff
    msg = "Running 'maxij_dophot.py' code on data in '" + night + "' folder..."
    loglist = addlog(msg, loglist)

    time0, msg = timestart()
    loglist = addlog(msg, loglist)

    pathnam = path + night + '/'
    apathnam = pathnam + 'aligned/'  # folder for aligned images

    # check to make sure all necessary files are available before starting
    msg = "Checking files..."
    loglist = addlog(msg,loglist)

    if not os.path.isdir(path+night):
        msg = "No folder exists for this night. Exiting..."
        loglist = addlog(msg, loglist)
        return None

    if not os.path.isdir(path+night+"/aligned"):
        msg = "No 'aligned' folder found, please re-run savealigned() function for this night."
        loglist = addlog(msg, loglist)

    if not os.path.isfile(path + night + "/ref_stars_xy.txt"):
        if not os.path.isfile(pathnam + "MeasurementsXY.txt"):
            msg = "No reference star position file found. \n" \
                  "Please save an AIJ measurement file as 'ref_stars_xy.txt' in the night folder. \n" \
                  "Exiting... \n"
            loglist = addlog(msg, loglist)
            return None
        else:
            msg = "Converting 'MeasurementsXY.txt' file to a 'ref_stars_xy.txt' file (rebinning)."
            loglist = addlog(msg, loglist)
            refsLG = pd.read_csv(pathnam + 'MeasurementsXY.txt',
                                 delimiter=r"\s+")
            refsLG2 = refsLG / 4
            msg = "Saving 'ref_stars_xy.txt' file to night folder."
            loglist = addlog(msg, loglist)
            refsLG2.to_csv(pathnam + 'ref_stars_xy.txt', index=None, sep=' ')

    if not os.path.isfile(pathnam+"data_"+night+".pkl"):
        msg = "No database called '" + "data_"+night+".pkl" + "' exists. \n " \
            "Please run initdb() to initialize the database, and run getshifts() and savealigned() before " \
            "running this code. \n " \
            "Exiting..."
        loglist = addlog(msg, loglist)
        return None

    fdate = datetime.now().strftime("%Y%m%d%H%M%S")
    msg = "Making backup of database as 'data_" + night + "_bck" + fdate + ".pkl"
    loglist = addlog(msg, loglist)
    shutil.copyfile(pathnam+"data_"+night+".pkl",pathnam+"data_"+night+"_bck"+fdate+".pkl")


    # load pandas df
    msg = "Loading database for " + night + "..."
    loglist = addlog(msg, loglist)
    dataf = pd.read_pickle(pathnam+"data_"+night+".pkl")

    fileID = dataf['fileID']
    fnames = dataf['filename']

    msg = "Total science frames:  " + str(len(fnames))
    loglist = addlog(msg, loglist)


    # read in AIJ ref star positions file. must have X(FITS) and Y(FITS) positions
    # (other columns should not cause errors though)
    msg = "Loading list of reference star positions..."
    loglist = addlog(msg, loglist)

    refs = pd.read_csv(pathnam + 'ref_stars_xy.txt',
                delimiter=r"\s+")

    # do aperture photometry on all reference stars,
    # using multiprocessing
    msg = "Doing photometry..."
    loglist = addlog(msg, loglist)
    p = Pool(6) # uses 6 cpu cores
    for results in p.imap_unordered(partial(multiphot,
                                            dataf = dataf,
                                            pathnam = pathnam,
                                            apathnam = apathnam,
                                            refs = refs),
                                    fileID, 25):
        # print "maxi j1820 flux at " + dataf.loc[results[0]]['filetime'] + ': '+ str(results[4][1])
        #assign results to things in dataf
        # results = id, pt, sg, [dx, dy], sky, flux
        fid = results[0]

        # guassian parameters
        dataf.loc[fid,'gauss_params_0_amplitude'] = results[1][0]
        dataf.loc[fid,'gauss_params_1_x0'] = results[1][1]
        dataf.loc[fid,'gauss_params_2_y0'] = results[1][2]
        dataf.loc[fid,'gauss_params_3_sigma_x'] = results[1][3]
        dataf.loc[fid,'gauss_params_4_sigma_y'] = results[1][4]
        dataf.loc[fid,'gauss_params_5_theta'] = results[1][5]
        dataf.loc[fid,'gauss_params_6_offset'] = results[1][6]
        dataf.loc[fid,'gauss_sigma_avg'] = results[2]

        dataf.loc[fid,'gauss_offset_dx'] = results[3][0]
        dataf.loc[fid,'gauss_offset_dy'] = results[3][1]

        # maxij + refstar photometry
        dataf.loc[fid,'phot_tyc'] = results[5][0]
        dataf.loc[fid,'phot_maxij'] = results[5][1]
        dataf.loc[fid,'phot_ref2'] = results[5][2]
        dataf.loc[fid,'phot_ref3'] = results[5][3]
        dataf.loc[fid,'phot_ref4'] = results[5][4]
        dataf.loc[fid,'phot_ref5'] = results[5][5]
        dataf.loc[fid,'phot_ref6'] = results[5][6]

        # sky photometry
        dataf.loc[fid,'sky_tyc'] = results[4][0]
        dataf.loc[fid,'sky_maxij'] = results[4][1]
        dataf.loc[fid,'sky_ref2'] = results[4][2]
        dataf.loc[fid,'sky_ref3'] = results[4][3]
        dataf.loc[fid,'sky_ref4'] = results[4][4]
        dataf.loc[fid,'sky_ref5'] = results[4][5]
        dataf.loc[fid,'sky_ref6'] = results[4][6]

    # pickle the pandas dataframe (should have ALL data now!)
    msg = "Saving photometry database to " + night + '/data_'+night+'.pkl'
    loglist = addlog(msg, loglist)
    dataf.to_pickle(pathnam + 'data_'+night+'.pkl')

    #get stop time and save screen output to log
    msg = timefinish(time0)
    loglist = addlog(msg, loglist)
    writelog(loglist, night, pathnam)

    gc.collect()
    return dataf

def multiphot(id, dataf, pathnam, apathnam, refs):
    filename = dataf.loc[id]['filename']
    afilename = filename.split('.')[0]+'_aligned.' + filename.split('.')[1]
    im0 = fits.getdata(apathnam + afilename)

    tx0 = np.rint(refs.loc[0]['X(FITS)'])
    ty0 = np.rint(refs.loc[0]['Y(FITS)'])

    win = 30.
    subim = get_subimage(im0,tx0,ty0,win) #tycho postage stamp subimage
    # plt.imshow(subim)
    # plt.show()

    # calculate gaussian parameters
    z, pt = gauss2dfit(subim, 2, 28, (2e4, 15.0, 14.0, 3., 3., 0., 0.))  # 2d Gaussian fit
    sg = (pt[3] + pt[4]) / 2.0  # get the sigma - average X, Y sigmas

    dx = pt[1]-win/2.
    dy = pt[2]-win/2.

    sky = np.zeros(len(refs.index))
    flux = np.zeros(len(refs.index))

    for i in refs.index:
        #for ref star positions rx, ry (with tycho star gaussian shifts):
        rx = np.rint(refs.loc[i]['X(FITS)']) + dx
        ry = np.rint(refs.loc[i]['Y(FITS)']) + dy

        sky[i], flux[i], skyplot = measure_star(im0, rx, ry, sg)


    return id, pt, sg, [dx, dy], sky, flux



if __name__ == "__main__":
    try:
        print dophot('test', path='./')
    except:
        print('An error occured.')
    # dophot('2018-03-28')