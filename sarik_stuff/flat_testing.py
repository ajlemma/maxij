#! /usr/bin/env python

import numpy as np
from astropy.io import fits
from glob import glob
from PIL import Image

master_dark = fits.getdata('masters/dark_1s_master.fits')
master_flat = fits.getdata('masters/flat_1s_master.fits')

files = sorted(glob('flat_1s_11*.tiff'))
sets_of = 11

for i in range(len(files)/sets_of):
    subset = files[(i*sets_of):(i*sets_of)+sets_of]
    print('Subset:', subset)
    final_array = []
    for j in range(len(subset)):
        final_array.append(np.array(Image.open(subset[j]), dtype=np.uint16) - master_dark)
    final_array = np.array(final_array, dtype=np.float32)
    subset_combined = np.median(final_array, axis=0)
    subset_normalized = subset_combined / np.median(subset_combined)
    pHDU = fits.PrimaryHDU(subset_normalized)
    pHDU.writeto("flat_%dcombined_set%d.fits" % (sets_of, i+1), overwrite=True)
    print('Max: %.2f' % np.max(master_flat - subset_normalized))
    print('Min: %.2f' % np.min(master_flat - subset_normalized))
    print('Std: %.2f' % np.std(master_flat - subset_normalized))
