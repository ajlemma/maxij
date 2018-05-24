# maxij image analysis module

from astropy.io import fits
from glob import glob
from os.path import basename
from natsort import natsorted
import numpy as np
import time
from datetime import datetime
import scipy.optimize as opt

def write_to_fits(f_name, data):
    # writes file as fits
    fits.PrimaryHDU(data).writeto(f_name, overwrite=True)

def get_filelist_maxi(pathname):
    # returns sorted list of filnames in pathname
    f = [basename(x) for x in glob(pathname + '*')]
    return natsorted(f)

def fix_bad_pixels(im0, badpix_list):
    # replaces the bad pixels for all pixels in list
    for pix in badpix_list:
        im0[pix] = 0.0
    return im0

def stack_maxi(pathname, n1, n2):
    # given a pathnam, stacks the images from n1 to n2 in the index
    maxi_files = get_filelist_maxi(pathname)
    n = n2 - n1  # number to stack

    for i in range(n1, n2):
        fnam = pathname + maxi_files[i]  # make the full filename
        im0 = fits.getdata(fnam)

        if i == n1:  # first one
            imf = np.zeros(im0.shape)
        imf = imf + im0

    return imf / n

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def get_sec(time_str):
    h, m, s = time_str.split('-')
    return int(h) * 3600 + int(m) * 60 + int(s)

def parse_time(fname):
    ss = fname.find('-')
    return fname[ss - 2:ss + 6]

def timestart():
    print "start: " + str(datetime.now())
    return time.time()

def timefinish(t_initial):
    t_now = time.time()
    print "finish: " + str(datetime.now())
    elapsed_time = t_now - t_initial
    print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time / 60.) + " min)"


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

    x = np.linspace(0, data.shape[0] - 1, data.shape[0])
    y = np.linspace(0, data.shape[0] - 1, data.shape[0])
    x, y = np.meshgrid(x, y)
    # print x.shape
    subx = x[lim1:lim2, lim1:lim2]
    suby = y[lim1:lim2, lim1:lim2]
    subdata = data[lim1:lim2, lim1:lim2]
    # print p_tmp
    pt, pc = opt.curve_fit(twoD_Gaussian, (subx, suby), subdata.ravel(), p0=p_tmp)

    z = twoD_Gaussian((x, y), *pt)
    data_fitted = z.reshape(data.shape)

    return data_fitted, pt

def get_subimage(im0, x,y, winsize):

    d = winsize/2. # half of subimage size

    xlow = int(x-d)
    xhigh = int(x+d)
    ylow = int(y-d)
    yhigh = int(y+d)

    return im0[ylow:yhigh, xlow:xhigh]  # select sub-image around nominal pre-calc'd Tycho pos'n

def get_grid(imshape):

    xx = np.arange(imshape[1])
    yy = np.arange(imshape[0])
    #XX, YY = np.meshgrid(xx, yy)
    return np.meshgrid(xx, yy)


def mkrad(imshape,xc,yc):

    x,y=get_grid(imshape)
    return ((x-xc)**2. + (y-yc)**2.)**0.5


def measure_star(im0,x,y,sg):

    # array of distances from star
    r = mkrad(im0.shape, x,y)

    # get median sky value in annulus centered 10 pixels from Tycho +-2-sigma
    r2 = r - 10 * sg
    sky = np.median(im0[np.abs(r2) < (2 * sg)])

    #sky annulus for plotting:
    skyplot = r
    skyplot[np.abs(r2) > (2 * sg)] = 0

    #photometer 3-sigma aperture & subtract sky
    flux = np.sum(im0[r < (3. * sg)]) - sky * len(im0[r < (3. * sg)].ravel())


    return sky, flux, skyplot


# sky annulus transparency colormap
cdict = {'red':   ((0.0,  1.0, 1.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  0.0, 0.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  0.8, 0.8)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  0.4, 0.4))}

cdict4 = cdict.copy()
a = .45
cdict4['alpha'] = ((0.0, 0.0, 0.0),
                   (0.5, a, a),
                   (1.0, a, a))
