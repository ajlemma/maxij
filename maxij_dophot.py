#######
''' maxij_dophot.py
by A Townsend
1.
'''
#######

import pandas as pd
# import numpy as np
from maxijdefs import *
from multiprocessing import Pool
from functools import partial
from astropy.io import fits
from matplotlib import pyplot as plt

plt.register_cmap(name='SkyCMAP', data=cdict4)
#######



def dophot(night, path='/media/amanda/demeter/maxi_j1820_070/',savecheck = 100):
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
    for results in p.imap_unordered(partial(multiphot,
                                            dataf = dataf,
                                            pathnam = pathnam,
                                            apathnam = apathnam,
                                            refs = refs,
                                            savecheck = savecheck),
                                    fileID, 25):
        print "maxi j1820 flux at " + dataf.loc[results[0]]['filetime'] + ': '+ str(results[4][1])
        #assign results to things in dataf

    #dataf.to_pickle(pathnam + 'photdata_'+night+'.pkl')

    timefinish(time0)
    return dataf

def multiphot(id, dataf, pathnam, apathnam, refs, savecheck):
    im0 = fits.getdata(apathnam + dataf.loc[id]['filename'])

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

    ## plot tycho star with an aperture drawn on for each image
    fig = plt.figure(facecolor = 'w')
    ax = plt.gca()
    ax.cla()  # clear things for fresh plot
    plt.imshow(np.arcsinh(im0+10000),cmap='gray')
    # ax.add_artist(plt.Circle((tx0 + dx, ty0 + dy), 3.*sg, color='y', fill=False))


    sky = np.zeros(len(refs.index))
    flux = np.zeros(len(refs.index))

    for i in refs.index:
        #for ref star positions rx, ry (with tycho star gaussian shifts):
        rx = np.rint(refs.loc[i]['X(FITS)']) + dx
        ry = np.rint(refs.loc[i]['Y(FITS)']) + dy

        sky[i], flux[i], skyplot = measure_star(im0, rx, ry, sg)

        #plot sky apertures:
        plt.imshow(skyplot, cmap='SkyCMAP')

        # to plot apertures:
        if i == 1:
            # maxij in blue
            ax.add_artist(plt.Circle((rx, ry), 3. * sg, color='b', fill=False))

            plt.text(rx, ry+10*sg, 'maxi j1820+070', ha="center", family='sans-serif', size=14, color = 'blue')
            plt.text(10, 235, 'BH f=' + str(np.rint(flux[i])) + ', sky=' + str(np.rint(sky[i])), family='sans-serif', size=14, color='blue')
        else:
            # the reference stars in red
            # name = 'refstar'+str(int(i)+1)
            ax.add_artist(plt.Circle((rx, ry), 3. * sg, color='r', fill=False))
            plt.text(rx, ry-5*sg, 'r' + str(i+1), ha="center", family='sans-serif', size=14, color = 'red')


    ## save checkimages with filenames:
    plt.text(10, 220, dataf.loc[id]['filename'], family='sans-serif', size=14, color='yellow')
    if int(id)%savecheck == 0:
        # plt.show()
        fig.savefig(pathnam + 'checkapertures_' + id + '.png')

    return id, pt, sg, [dx, dy], sky, flux



if __name__ == "__main__":
    dophot('test', path='./',savecheck = 5)
