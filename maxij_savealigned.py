#######
''' maxij_savealigned.py
by A Townsend
1. read in imshifts file created by maxij_getshifts
2. open each image
3. shift image
4. rebin image
5. save new image in aligned folder
'''
#######

from maxijdefs import *
import pandas as pd
from scipy.ndimage.interpolation import shift
from astropy.io import fits
from multiprocessing import Pool
from functools import partial
#######

def savealigned(night, n1 = 0, n2 = 'max', path='/media/amanda/demeter/maxi_j1820_070/'):

    time0 = timestart()
    print "Processing images for " + night

    pathnam = path+night+'/'
    scipathnam = pathnam + 'science/'  # folder with science images
    apathnam = pathnam + 'aligned/'  # folder for aligned images

    print "getting list of filenames from " + night + "/science ..."
    fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam
    print "total science frames:  " + str(len(fnames))

    print "getting image shift data from " + night + '/imshifts_' + night + ".txt ..."
    imshifts = pd.read_csv(pathnam + 'imshifts_' + night + '.txt',
                             header=None,
                             names = ["fname", "y", "x", "ampcorr"],
                             index_col='fname',
                             delimiter=r"\s+")


    ## define last image to stack
    if n2 == 'max':
        n2 = len(fnames)
    else:
        n2 = int(n2)

    p = Pool(6)
    print "shifting images..."
    p.map(partial(align_rebin, scipathnam = scipathnam,apathnam = apathnam,imshifts = imshifts), fnames[n1:n2],35)

    print "Aligned & rebinned images are saved at " + night + '/aligned with the suffix "_aligned"'

    timefinish(time0)


def align_rebin(filename,scipathnam,apathnam,imshifts):
    # print filename
    im0 = fits.getdata(scipathnam+filename)
    im1 = shift(im0, (imshifts.loc[filename]["y"], imshifts.loc[filename]["x"]), mode='constant', cval=0.0)
    im2 = rebin(im1, (im1.shape[0] / 4, im1.shape[1] / 4))
    write_to_fits(apathnam + filename.split('.')[0]+'_aligned.' + filename.split('.')[1], im2)  # write to aligned subdir, add suffix "_aligned"
    return


if __name__ == "__main__":
    savealigned('test', n2 = 10, path = './')
