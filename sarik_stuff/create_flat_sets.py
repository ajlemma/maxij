#! /usr/bin/env python

import numpy as np
import sys
from astropy.io import fits
from glob import glob
from PIL import Image


# Create the master dark array
dark_image_names = glob('dark_1s*.tiff')
dark_image_array = np.zeros([960, 1280, len(dark_image_names)], dtype=np.uint16)
for i in range(len(dark_image_names)):
    dark_image_array[:, :, i] = np.array(Image.open(dark_image_names[i]), dtype=np.uint16)
master_dark = np.array(np.median(dark_image_array, axis=-1), dtype=np.float32)

# Create flats by combining every 'sets_of' flats
files = sorted(glob('flat_1s_*_*.tiff'))
if len(sys.argv) == 2:
    sets_of = int(sys.argv[-1])
else:
    sets_of = 5
median_of_flats_sets = []
for i in range(len(files)/sets_of):
    subset = files[(i*sets_of):(i*sets_of)+sets_of]
    print('Subset:', subset)
    final_array = []
    for j in range(len(subset)):
        final_array.append(np.array(Image.open(subset[j]), dtype=np.uint16) - master_dark)
    final_array = np.array(final_array, dtype=np.float32)
    median_combined_subset = np.median(final_array, axis=0)
    median_combined_normalized = median_combined_subset / np.median(median_combined_subset)
    median_of_flats_sets.append(median_combined_normalized)
    pHDU = fits.PrimaryHDU(median_combined_normalized)
    pHDU.writeto("flat_%dcombined_set%d.fits" % (sets_of, i+1), overwrite=True)

flat_median_of_sets = np.median(np.array(median_of_flats_sets), axis=0)
pHDU = fits.PrimaryHDU(flat_median_of_sets)
pHDU.writeto('flat_1s_median_of_sets.fits', overwrite=True)

print('done!')
