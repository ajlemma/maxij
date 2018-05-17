# maxij image analysis module

from astropy.io import fits
from glob import glob
from os.path import basename
from natsort import natsorted
import numpy as np


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