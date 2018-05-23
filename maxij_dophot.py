#######
''' maxij_dophot.py
by A Townsend
1.
'''
#######

import pandas as pd
import numpy as np
from maxijdefs import *
from multiprocessing import Pool
from functools import partial
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
#######

def dophot(night, path='/media/amanda/demeter/maxi_j1820_070/'):
    time0 = timestart()
    print "Analyzing images for " + night + "..."

    pathnam = path + night + '/'
    apathnam = pathnam + 'aligned/'  # folder for aligned images

    print "Loading database for " + night + "..."
    dataf = pd.read_pickle(pathnam + 'data_'+night+'.pkl')
    fileID = dataf['fileID']

    print "Loading list of reference star positions..."
    refs = pd.read_csv(pathnam + 'ref_stars_xy.txt',
                delimiter=r"\s+")
    # print refs

    print "Doing photometry..."
    p = Pool(6)
    for results in p.imap_unordered(partial(multiphot, dataf = dataf, apathnam = apathnam, refs = refs), fileID, 25):
        print results[1]
        #assign results to things in dataf

    #dataf.to_pickle(pathnam + 'photdata_'+night+'.pkl')

    timefinish(time0)
    return dataf

def multiphot(id, dataf, apathnam, refs):
    im0 = fits.getdata(apathnam + dataf.loc[id]['filename'])

    tx = np.rint(refs.loc[0]['X(FITS)'])
    ty = np.rint(refs.loc[0]['Y(FITS)'])

    win = 30.
    subim = get_subimage(im0,tx,ty,win) #tycho postage stamp subimage
    # plt.imshow(subim)
    # plt.show()

    # calculate gaussian parameters
    z, pt = gauss2dfit(subim, 2, 28, (2e4, 15.0, 14.0, 3., 3., 0., 0.))  # 2d Gaussian fit
    sg = (pt[3] + pt[4]) / 2.0  # get the sigma - average X, Y sigmas

    # print pt[1]-win/2., pt[2]-win/2.

    ax = plt.gca()
    ax.cla()  # clear things for fresh plot
    plt.imshow(im0,cmap='gray')
    ax.add_artist(plt.Circle((pt[1]-win/2.+tx, pt[2]-win/2.+ty), 3.*sg, color='y', fill=False))
    ax.set_xlim((tx-win, tx+win))
    ax.set_ylim((ty-win, ty+win))
    plt.show()



    return pt, sg

if __name__ == "__main__":
    dophot('test', path='./')
