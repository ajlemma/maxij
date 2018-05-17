#######
''' maxij_makeref.py
create a reference image for the night
'''
#######

from maxijdefs import *
import matplotlib.pyplot as plt

#######

def makeref(night,im0=60,im1=120,path='/media/amanda/demeter/maxi_j1820_070/'):
    time0 = timestart()
    print "Processing images for " + night

    pathnam = path+night+'/'
    scipathnam = pathnam + 'science/'  # folder with science images

    print "getting list of filenames from " + night + "science ..."
    fnames = get_filelist_maxi(scipathnam)  # list of sorted filenames for sci images in pathnam
    print "total science frames:  " + str(len(fnames))

    ##create & plot (to check) a 60s reference image
    ## should check this by-eye prior to beginning shifts & save it to
    print "stacking reference image..."
    zref = stack_maxi(scipathnam,im0,im1) #, fix_list)
    write_to_fits(pathnam+'ref_stack.fits',zref)
    print "reference image is saved at " + night + "ref_stack.fits"

    # to show/save/check reference image:
    f1 = plt.figure()
    plt.imshow(zref,vmin=1e3,vmax=3e3)
    plt.colorbar()
    plt.show()
    f1.savefig('zref'+night+'.png')

    timefinish(time0)

if __name__ == "__main__":
    makeref('2018-03-28')
