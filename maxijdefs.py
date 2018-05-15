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

def stack_maxi(pathname, n1, n2, badpix_list):
    # given a pathnam, stacks the images from n1 to n2 in the index
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
        im1 = fix_bad_pixels(im0, badpix_list)
        fID.close()
        if i == n1:  # first one
            imf = np.zeros(im1.shape)
        imf = imf + im1

    return imf / n

# def rd_shifts_maxi(fnam):
#     # reads the shifts file
#
#     f = open(fnam, "r")
#     filenams = list()
#     xlist = list()
#     ylist = list()
#     clist = list()
#
#     while True:
#         line = f.readline()
#
#         if not line: break
#         jnk = line.split()
#
#         if len(jnk) == 4:
#             filenams.append(jnk[0])
#             xlist.append(float(jnk[2]))
#             ylist.append(float(jnk[1]))
#             clist.append(float(jnk[3]))
#
#         if len(jnk) == 5:
#             filenams.append(jnk[1])
#             xlist.append(float(jnk[3]))
#             ylist.append(float(jnk[2]))
#             clist.append(float(jnk[4]))
#
#     f.close()
#
#     return filenams, np.array(ylist), np.array(xlist), np.array(clist)
