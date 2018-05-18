from maxij_makeref import *
from maxij_getshifts import *

night = '2018-05-05'
# night_path = '/media/amanda/demeter/maxi_j1820_070/'

#demeter is the default but should be changed for other computers, this is the folder with the night folder in it.
#night folder should contain flat-fielded & dark subtracted images in a folder called "science"
# as well as masters in a folder called "masters" and an os timestamps file

# make a refence image for correlation & finding star positions
# reference gets saved as ref_stack.fits in night directory
makeref(night,1000, 1015)

# calculates shifts via cross-correlation for each image wrt to the reference image created above
getshifts(night)
