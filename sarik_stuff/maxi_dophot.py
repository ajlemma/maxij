#! /usr/bin/env python

# The code below was written by Dr. Steve Eikenberry with some modifications by Sarik Jeram.

import multiprocessing as mp
import numpy as np
from astropy.io import fits
from glob import glob


def write_to_fits(f_name, data):
    fits.PrimaryHDU(data).writeto(f_name, overwrite=True)


def fix_328(im0):
    # replaces the 2 bad pixels in 28-March-2018 data on MAXI J1820
    im0[437, 916] = 0.0
    im0[893, 353] = 0.0

    return im0


def get_filelist_maxi(pathname):
    # returns the list of maxij1820 fits files in the directory pathname
    return glob('%s/maxij_1s*.fits' % pathname)


def stack_maxi(pathname, n1, n2):
    # given a pathname, stacks the images from n1 to n2 in the index

    maxi_files = get_filelist_maxi(pathname)

    n = n2 - n1  # number to stack
    nct = 0
    ict = 0

    for i in range(n1, n2):
        fnam = pathname + maxi_files[i]  # make the full filename
        # print fnam
        fID = fits.open(fnam)
        im0 = fID[0].data
        # print im0.shape
        im1 = fix_328(im0)
        fID.close()
        if i == n1:  # first one
            imf = np.zeros(im1.shape)
        imf = imf + im1

    return imf / n


def stack_maxi_rebinned(pathname, n1, n2):
    # given a pathname, stacks the images from n1 to n2 in the index

    maxi_files = get_filelist_maxi(pathname)

    n = n2 - n1  # number to stack

    for i in range(n1, n2):
        fnam = pathname + maxi_files[i]  # make the full filename
        # print fnam
        im0 = fits.getdata(fnam)
        # print im0.shape
        # im1=fix_328(im0)
        if i == n1:  # first one
            imf = np.zeros(im0.shape)
        imf += im0

    return imf / n


def find_shift_maxi(im0, imref):  # OBSOLETE
    # finds the shift of im0 with respect to imref
    from scipy.signal import fftconvolve

    dum0 = fftconvolve(imref, imref[::-1, ::-1])
    dum1 = fftconvolve(im0, imref[::-1, ::-1])  # cross-correlate them
    dum1 = dum1 - dum0
    xy = zip(*np.where(dum1 == np.max(dum1)))  # get index of max)
    xc = xy[0][1]
    yc = xy[0][0]

    return (dum1.shape[0] / 2) - yc, (dum1.shape[1] / 2) - xc, np.max(dum1) / 1e10  # required pixel shifts


def find_shift_maxi_wcorr(im0, imref, corr_ref):
    # finds the shift of im0 with respect to imref
    # variant uses pre-calculated correlation baseline (so you don't have to recalculate it every time)
    from scipy.signal import fftconvolve

    ima = np.zeros(im0.shape)
    imb = im0 - np.median(im0)
    ima[imb > 1e3] = imb[imb > 1e3]

    imrefa = np.zeros(imref.shape)
    imrefa[imref > 1e3] = imref[imref > 1e3]
    dum1 = fftconvolve(ima, imrefa[::-1, ::-1])  # cross-correlate them
    dum1 = dum1 - corr_ref
    # plt.imshow(dum1)
    # plt.colorbar()
    # plt.show()
    xy = zip(*np.where(dum1 == np.max(dum1)))  # get index of max)
    xc = xy[0][1]
    yc = xy[0][0]

    return (dum1.shape[0] / 2) - yc, (dum1.shape[1] / 2) - xc, np.max(dum1) / 1e10  # required pixel shifts


def do_shift_maxi(im0, imref):  # OBSOLETE - needs updating
    # finds the shift of im0 with respect to imref
    # then SHIFTS im0 to match imref
    from scipy.ndimage.interpolation import shift

    return shift(im0, find_shift_maxi(im0, imref)[0:2], mode='constant', cval=0.0)


def shift_stack_maxi(pathname, n1, n2):  # OBSOLETE - needs updating
    # given a pathname, stacks the images from n1 to n2 in the index
    # variant to SHIFT the images while stacking on second pass

    imref = stack_maxi(pathname, n1, n2)  # first pass
    maxi_files = get_filelist_maxi(pathname)

    n = n2 - n1  # number to stack
    nct = 0
    ict = 0

    for i in range(n1, n2):
        fnam = pathname + maxi_files[i]  # make the full filename
        # print fnam
        fID = fits.open(fnam)
        im0 = fID[0].data
        im1 = fix_328(im0)
        fID.close()
        if i == n1:  # first one
            imf = np.zeros(im1.shape)
        im2 = do_shift_maxi(im1, imref)
        imf = imf + im2

    return imf / n


def write_shifts_output_maxi(pathname, n1, n2, imref, onam):
    # given a pathname, indices of files, and reference image
    # prints name, shifts, correlation strength
    from scipy.signal import fftconvolve

    maxi_files = get_filelist_maxi(pathname)
    fout = open(onam, "w")

    n = n2 - n1  # number to stack
    nct = 0
    ict = 0

    imrefa = np.zeros(imref.shape)
    imrefa[imref > 1e3] = imref[imref > 1e3]
    dum1 = fftconvolve(imrefa, imrefa[::-1, ::-1])  # cross-correlate them
    for i in range(n1, n2):
        fnam = pathname + maxi_files[i]  # make the full filename
        # print fnam
        im0 = fits.getdata(fnam)
        im1 = fix_328(im0)
        a = find_shift_maxi_wcorr(im1, imref, dum1)
        print a
        fout.write('%d ' % i + maxi_files[i] + ' %f ' % a[0] + ' %f' % a[1] + ' %f' % a[2] + '\n')

    fout.close()
    return n


def rebin(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def rd_shifts_maxi(fnam):
    # reads the shifts file

    f = open(fnam, "r")
    filenams = list()
    xlist = list()
    ylist = list()
    clist = list()

    while True:
        line = f.readline()
        if not line: break
        jnk = line.split()

        if len(jnk) == 5:
            filenams.append(jnk[1])
            xlist.append(float(jnk[3]))
            ylist.append(float(jnk[2]))
            clist.append(float(jnk[4]))

    f.close()

    return filenams, np.array(ylist), np.array(xlist), np.array(clist)


def get_grid(imshape):
    xx = np.arange(imshape[1])
    yy = np.arange(imshape[0])
    # XX, YY = np.meshgrid(xx, yy)
    return np.meshgrid(xx, yy)


def mkrad(imshape, xc, yc):
    x, y = get_grid(imshape)
    return ((x - xc) ** 2 + (y - yc) ** 2) ** 0.5


def centroid(subim):
    x, y = get_grid(subim.shape)
    xc = np.sum(x * subim) / np.sum(subim)
    yc = np.sum(y * subim) / np.sum(subim)

    return xc, yc


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    # g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
    # + c*((y-yo)**2)))

    g = amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                              + c * ((y - yo) ** 2)))
    return g.ravel()


def gauss2dfit(data, lim1, lim2, p_tmp):
    # data = 2d square array
    # lim2, lim2 = subarray bounds for fitting (also square)
    # guess = initial guess (amplitude, xo, yo, sigma_x, sigma_y, theta, offset)

    import scipy.optimize as opt

    x = np.linspace(0, data.shape[0] - 1, data.shape[0])
    y = np.linspace(0, data.shape[0] - 1, data.shape[0])
    x, y = np.meshgrid(x, y)
    print x.shape
    subx = x[lim1:lim2, lim1:lim2]
    suby = y[lim1:lim2, lim1:lim2]
    subdata = data[lim1:lim2, lim1:lim2]
    print p_tmp
    pt, pc = opt.curve_fit(twoD_Gaussian, (subx, suby), subdata.ravel(), p0=p_tmp)

    z = twoD_Gaussian((x, y), *pt)
    data_fitted = z.reshape(data.shape)

    return data_fitted, pt


def apply_shifts_maxi(pathname, fnam):
    # SHIFTS im0 to match imref
    from scipy.ndimage.interpolation import shift

    flist, y, x, cc = rd_shifts_maxi(pathname + fnam)

    for i in range(len(flist)):
        print i
        fID = fits.open(pathname + flist[i])
        jnk = fID[0].data
        im0 = fix_328(jnk)
        im1 = shift(im0, (y[i], x[i]), mode='constant', cval=0.0)
        im2 = rebin(im1, (im1.shape[0] / 4, im1.shape[1] / 4))
        write_to_fits(pathname + 'Aligned/A' + flist[i], im2)  # write to Aligned subdir, add prefix "A"

    return i


def lc_detrend(photin):
    photout = np.zeros(len(photin))
    for i in range(len(photin)):
        i1 = i - 500
        i2 = i + 500
        if i1 < 0:
            i1 = 0
        if i2 > len(photin):
            i2 = len(photin)
        photout[i] = photin[i] - np.median(photin[i1:i2])

    return photout


def maxi_phot(flist, onam):
    import matplotlib.pyplot as plt
    pathname = 'D:/NICER/MAXI_J1820/2018-03-28/reduction/Aligned/'
    # does photometry on Tycho and MAXIJ only
    # uses Gaussian fit to Tycho to get sigma and offset from nominal aligned position for each image
    # uses a pre-calculated offset to MAXI J

    # writes output to file onam

    gparams = np.zeros([len(flist), 7])  # Gaussian params array
    phot = np.zeros(len(flist))  # Tycho phot array
    phot2 = np.zeros(len(flist))  # MAXI phot array
    fout = open(onam, "w")

    for i in range(len(flist)):  # loop over files
        jnk = fits.open(pathname + 'A' + flist[i])  # read in the file
        im0 = jnk[0].data
        jnk.close()

        print i
        subim = im0[25:55, 95:125]  # select sub-image around nominal pre-calc'd Tycho postn
        z, pt = gauss2dfit(subim, 2, 28, (2e4, 15.0, 14.0, 3., 3., 0., 0.))  # 2d Gaussian fit
        gparams[i, :] = pt[:]  # save the Gaussian parameters
        sg = (pt[3] + pt[4]) / 2.0  # get the sigma - average X, Y sigmas
        r = mkrad(im0.shape, (pt[1] + 95.), (pt[2] + 25.))  # make an array of distances from Tycho in this image
        r2 = r - 10 * sg  # subtract 10 from this
        sky = np.median(
            im0[np.abs(r2) < (2 * sg)])  # get median sky value in annulus centerd 10 pixels from Tycho +-2-sigma
        phot[i] = np.sum(im0[r < (3. * sg)]) - sky * len(
            im0[r < (3. * sg)].ravel())  # photometer 3-sigma aperture, subtracting sky

        subim = im0[140:170, 92:122]  # MAXI J1820
        rm = mkrad(im0.shape, (pt[1] + 95.0 - 3.4), (pt[2] + 25.0 + 116.9))  # MAXI offsets
        r2m = r - 10 * sg
        sky = np.median(im0[np.abs(r2m) < (2 * sg)])
        phot2[i] = np.sum(im0[rm < (3. * sg)]) - sky * len(im0[rm < (3. * sg)].ravel())
        if i == 63:
            rr = np.zeros(r.shape)
            rr[r < (3. * sg)] = 1.0
            imr = im0 * rr
            plt.imshow(imr)
            plt.colorbar()
            plt.show()
            rr = np.zeros(r.shape)
            rr[rm < (3. * sg)] = 1.0
            imr = im0 * rr
            plt.imshow(imr)
            plt.colorbar()
            plt.show()

        fout.write('%d ' % i + flist[i] + ' %f ' % float(phot[i]) + ' %f ' % float(phot2[i]) + ' %f ' % gparams[
            i, 0] + ' %f ' % gparams[i, 1] + ' %f ' % gparams[i, 2] + ' %f ' % gparams[i, 3] + ' %f ' % gparams[
                       i, 4] + ' %f ' % gparams[i, 5] + ' %f \n' % gparams[i, 6])
    fout.close()
    return gparams, phot, phot2
