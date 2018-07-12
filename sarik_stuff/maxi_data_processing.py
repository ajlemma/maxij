#! /usr/bin/env python2
# Written by Alma and Sarik

import multiprocessing as mp
import numpy as np
import os
import sys
from astropy.io import fits
from glob import glob
from PIL import Image


if len(sys.argv) == 1:
    # Path to all the image files (i.e. /home/sarik/maxi/2018-05-07_observation/data); ensure the path is valid, that
    # the timestamps file exists, and that the data_location path doesn't end in a slash
    data_location = raw_input('Enter the full path to the data folder '
                              '(ex. /home/sarik/maxi/2018-05-07_observation/data): ')
    if not os.path.exists(data_location):
        raise Exception('The data folder you entered does not exist. Double-check the path!')
    if not os.path.exists('%s/os_timestamps.txt' % data_location):
        raise Exception('The timestamps file "os_timestamps.txt" was not found in the data directory!')
    if data_location[-1] == '/':  # make sure file data folder path doesn't in a slash
        data_location = data_location[:-1]

    # Path to root directory where final image folder is going to be created
    final_location = raw_input('Enter the full path to the final folder '
                               '(ex. /home/sarik/bast/MAXI_J1820_070_REDUCED): ')
    if not os.path.exists(final_location):
        raise Exception('The destination folder you entered does not exist. Double-check the path!')

elif len(sys.argv) == 3:
    data_location, final_location, exp_time = sys.argv[1], sys.argv[2], '1s'

elif len(sys.argv) == 4:
    data_location, final_location, exp_time = sys.argv[1], sys.argv[2], sys.argv[3]

else:
    raise Exception('Incorrect number of command-line arguments specified. Please double check your command!')

# Output folder name; create the necessary folders
folder_date = data_location.split('/data')[0].split('/')[-1].replace('_observation', '')
if not os.path.exists('%s/%s' % (final_location, folder_date)):
    os.mkdir('%s/%s' % (final_location, folder_date))
    masters_folder = '%s/%s/masters' % (final_location, folder_date)
    os.mkdir(masters_folder)
    science_folder = '%s/%s/science' % (final_location, folder_date)
    os.mkdir(science_folder)
    os.system('/bin/cp %s/os_timestamps.txt %s/%s/timestamps_%s.txt' % (data_location, final_location, folder_date,
                                                                        folder_date))
else:
    raise Exception('A folder already exists for this date in the destination! Back up, delete, or rename that folder '
                    'before you continue!')

# Remove the space from the filenames in that folder and replace it with an underscore; >/dev/null to suppress output
print('Replacing the space in the filenames with an underscore...')
os.system('/bin/bash -c \'for i in %s/*tiff; do mv "$i" "${i//\ /_}"; done\' >/dev/null 2>&1' % data_location)

# Separates images by file name/type
dark_image_names = glob('%s/dark_%s*.tiff' % (data_location, exp_time))
flat_image_names = glob('%s/flat_1s*.tiff' % data_location)
maxij_image_names = glob('%s/maxij_%s*.tiff' % (data_location, exp_time))

# Create log file to flag bad pixels
logfile_name = '%s/%s/badpixels.txt' % (final_location, folder_date)
logfile = open(logfile_name, 'a')
logfile.write('# Image Name                  X   Y  Flag\n')

# Create master dark for 1s images
print("Creating master dark...")
dark_image_array = np.zeros([960, 1280, len(dark_image_names)], dtype=np.uint16)
for idx in range(len(dark_image_names)):
    dark_image_array[:, :, idx] = np.array(Image.open(dark_image_names[idx]), dtype=np.uint16)
dark_combined = np.array(np.median(dark_image_array, axis=-1), dtype=np.float32)
dark_final_name = 'dark_1s_master.fits'
try:
    y_indcs, x_indcs = np.where(dark_combined == 65535.0)
    for idx in range(len(y_indcs)):
        logfile.write('%-30s  %4d %4d SAT.\n' % (dark_final_name, y_indcs[idx], x_indcs[idx]))
except ValueError:
    pass
pHDU = fits.PrimaryHDU(dark_combined)
pHDU.writeto('%s/%s' % (masters_folder, dark_final_name), overwrite=True)
print("done!")

# Create master flat by median combining the median combined sets of 'sets_of' flat images
print("Creating master flat...")
median_of_flats_sets = []
sets_of = 5
for idx in range(len(flat_image_names) / sets_of):
    subset = flat_image_names[(idx * sets_of):(idx * sets_of) + sets_of]
    final_array = []
    for j in range(len(subset)):
        final_array.append(np.array(Image.open(subset[j]), dtype=np.uint16) - dark_combined)
    final_array = np.array(final_array, dtype=np.float32)
    median_combined_subset = np.median(final_array, axis=0)
    median_combined_normalized = median_combined_subset / np.median(median_combined_subset)
    median_of_flats_sets.append(median_combined_normalized)
flat_combined = np.median(np.array(median_of_flats_sets), axis=0)
flat_final_name = 'flat_1s_master.fits'
try:
    y_indcs, x_indcs = np.where(flat_combined == 0.0)
    for idx in range(len(y_indcs)):
        logfile.write('%-30s  %4d %4d ZERO\n' % (flat_final_name, y_indcs[idx], x_indcs[idx]))
        flat_combined[y_indcs[idx], x_indcs[idx]] = 1.0  # change values of 0 in the flat image to 1
except ValueError:
    pass
pHDU = fits.PrimaryHDU(flat_combined)
pHDU.writeto('%s/%s' % (masters_folder, flat_final_name), overwrite=True)
print("done!")
logfile.close()

# Total number of science images to create
image_total = len(maxij_image_names)


# Subtract master dark and divide by normalized master flat
def create_science(i):
    maxij_image_array = np.array(Image.open(maxij_image_names[i]), dtype=np.float32)  # load image data
    try:
        y_indcs, x_indcs = np.where(maxij_image_array == 65535.0)
        for idx in range(len(y_indcs)):
            os.system("/usr/bin/printf '%-30s  %4d %4d SAT.\n' >> %s" % (maxij_image_names[i].split('/')[-1], y_indcs[idx],
                                                                      x_indcs[idx], logfile_name))
    except ValueError:
        pass
    maxij_image_array -= dark_combined  # subtract master dark
    maxij_image_array /= flat_combined  # divide master flat
    pHDU = fits.PrimaryHDU(maxij_image_array)
    final_image_path = '%s/%s' % (science_folder, maxij_image_names[i].split('/')[-1].replace('.tiff', '_reduced.fits'))
    pHDU.writeto(final_image_path, overwrite=True)
    print('\033[FProgress: %2d%% (%05d/%05d)' % (int(100 * (i+1.) / image_total), (i+1), image_total))
    return


# Parallelize the create_science task for speed increase!
print("Creating final science images...\n")
pool = mp.Pool()
pool.map(create_science, range(image_total), chunksize=1)
logfile.close()

print("done! Created %d final images for %s!" % (len(maxij_image_names), folder_date))
