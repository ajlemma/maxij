#######
''' maxij_getshifts.py
Original Code by S Eikenberry, modified by A Townsend
1. read in reference image
2. cross correlate it with itself to get a reference baseline
3. cross correlate the zref with each individual science image for the night
using parallel processing
4. save image shifts as imshifts_date.txt file
'''
#######

from maxijdefs import *
from scipy.signal import fftconvolve
from astropy.io import fits
from multiprocessing import Pool
from functools import partial
#######

def getshifts(night, n1 = 0, n2 = 'max',path='/media/amanda/demeter/maxi_j1820_070/'):
    # n2 can also be an integer

    time0 = timestart()
    print "Processing images for " + night

    pathnam = path+night+'/'
    scipathnam = pathnam + 'science/'  # folder with science images

    print "getting list of filenames from " + night + "/science ..."
    fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam
    print "total science frames:  " + str(len(fnames))

    ## get rid of low-level wave information for reference image & cross-correlate it with itself:
    print "getting correlation reference..."
    imref = fits.getdata(pathnam+'ref_stack.fits') #read in stacked image
    imref[imref < 1e3] = 0
    corr_ref = fftconvolve(imref, imref[::-1, ::-1])  # cross-correlate them

    ## define last image to stack
    if n2 == 'max':
        n2 = len(fnames)
    else:
        n2 = int(n2)


    p = Pool(6)
    print "calculating individual image shifts..."
    with open(pathnam + 'imshifts_' + night + '.txt', 'w') as f:
        for a in p.imap_unordered(partial(get_shifts, imref = imref,corr_ref = corr_ref,scipathnam = scipathnam),
                                  fnames[n1:n2], 25):
            output = (a[0] + ' %f ' % a[1] + ' %f' % a[2] + ' %f' % a[3] + '\n')
            # print output
            f.write(output)
    print "image shifts are saved at " + night + 'imshifts_' + night + '.txt'

    timefinish(time0)


def get_shifts(science_image,imref,corr_ref,scipathnam):
    # finds the shift of science_image with respect to reference_image
    # uses pre-calculated correlation_reference baseline (so you don't have to recalculate it every time)
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
    return (science_image, (correlation2.shape[0] / 2) - int(yc), (correlation2.shape[1] / 2) - int(xc), np.max(correlation2) / 1e10)  # required pixel shifts


if __name__ == "__main__":
    getshifts('test', n2 = 10, path = './')
