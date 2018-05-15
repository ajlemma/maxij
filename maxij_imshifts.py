#######
''' maxij_imshifts.py
Original Code by S Eikenberry, modified by A Townsend
1. make & save a 1-minute zref image (check this!!)
2. cross correlate it with itself to get a reference baseline
3. cross correlate the zref with each individual science image for the night
using parallel processing
4. save image shifts as imshifts_date.txt file
'''
#######

from maxijdefs import *
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
import time
from datetime import datetime
from multiprocessing import Pool

#######
night = '2018-05-05'
## path to the directory with the MAXIJ native-resolution fits files:
pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'
scipathnam = pathnam + 'science/' #folder with science images

time0 = time.time()
print "start: " + str(datetime.now())

# ## list of bad pixels; this code is obsolete
# ## however, as these pixels are saturated in EVERY image
# ## it seems prudent to leave
# fix_list = [[437, 916], [893, 353]]
#


fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam
# print len(fnames)

##create & plot (to check) a 60s reference image
## should check this by-eye prior to beginning shifts & save it to
zref = stack_maxi(scipathnam,120,180) #, fix_list)
write_to_fits(pathnam+'ref_stack.fits',zref)

## to show/save/check reference image:
# f1 = plt.figure()
# plt.imshow(zref,vmin=1e3,vmax=3e3)
# plt.colorbar()
# plt.show()
# f1.savefig('zref'+night+'.png')


## get rid of low-level wave information for reference image & cross-correlat it with itself:
imref = zref
imref[imref < 1e3] = 0
corr_ref = fftconvolve(imref, imref[::-1, ::-1])  # cross-correlate them

## first and last images to stack, and total no. of images
n1 = 0                          # first
n2 = len(fnames)                # last
n_tot = n2 - n1                 # total

def get_shifts(science_image):
    # finds the shift of science_image with respect to reference_image
    # uses pre-calculated correlation_reference baseline (so you don't have to recalculate it every time)
    reference_image = imref
    correlation_reference = corr_ref
    fnam = scipathnam + science_image  # make the full filename
    im0 = fits.getdata(fnam)
    # im0 = fix_bad_pixels(im0, fix_list)

    im0 = im0 - np.median(im0)
    im0[im0 < 1e3] = 0

    correlation1 = fftconvolve(im0, reference_image[::-1, ::-1])  # cross-correlate them
    correlation2 = correlation1 - correlation_reference

    xy = zip(*np.where(correlation2 == np.max(correlation2)))  # get index of max)
    xc = xy[0][1]
    yc = xy[0][0]

    # print xc, yc
    return (science_image, (correlation2.shape[0] / 2) - int(yc), (correlation2.shape[1] / 2) - int(xc), np.max(correlation2) / 1e10)  # required pixel shifts

## Non-parallelized shifts version:
## for i in range(n_tot):
##     a = get_shifts(scipathnam,fnames[i],fix_list,imref,corr_ref)
##     print ('%d ' % i + ' ' + fnames[i] + ' %f ' % a[0] + ' %f' % a[1] + ' %f' % a[2] + '\n')

## Parallelized version--still takes a bit over an hour on 6 processors,
## So consider carefully before uncommenting:
# p = Pool(6)
# with open(pathnam+'imshifts_'+night+'.txt','w') as f:
#     for a in p.imap_unordered(get_shifts,fnames[n1:n2],25):
#         output = (a[0] + ' %f ' % a[1] + ' %f' % a[2] + ' %f' % a[3] + '\n')
#         # print output
#         f.write(output)

time1 = time.time()
print "finish: " + str(datetime.now())
elapsed_time = time1 - time0
print "time elapsed: " + str(elapsed_time) + " s (" + str(elapsed_time/60.) + " min)"
