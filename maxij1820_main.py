#####
from maxij_getshifts import *
from maxij_savealigned import *
from maxij_initdb import *
from maxij_dophot import *
import time
import gc
#####


def runall(night):
    # note: need ref image for getshifts (AIJ)
    # and get x,y measurements before running dophot (AIJ)
    # (don't forget to save a check/reference png)
    print
    print
    print "Running code for " + night
    print
    try:
        initdb(night)  # initialize the df
        gc.collect()
        print
        getshifts(night)  # get shifts & save them to the df
        gc.collect()
        time.sleep(60)
        print
        savealigned(night)  # align & rebin individual ims, save to aligned folder
        gc.collect()
        time.sleep(60)
        print
        dophot(night) # do photometry & save values to the df
        gc.collect()
        time.sleep(60)
        print
    except:
        print
        print night + " failed, check log for details"
        gc.collect()



# #example
# night = 'test'
# runall(night)


night = '2018-03-27'

# night = '2018-03-28' # this night is done!
# runall(night)

night = '2018-03-29' # 1s ims for this night are done but it ought to be redone for 500ms ims
# runall(night)
# time.sleep(60)

night = '2018-03-30' # 1s ims are done for this night but it ought to be redone for 2s ims
# runall(night)
# time.sleep(60)

# night = '2018-04-06' # this night is done!
# runall(night)
# time.sleep(60)
# dophot(night)

# night = '2018-04-12' # this night is done!
# runall(night)
# time.sleep(60)

# night = '2018-04-13' # this night is done!
# runall(night)
# time.sleep(60)
# savealigned(night)
# dophot(night)

night = '2018-04-14'

night = '2018-04-17' # running getshifts, savealigned, dophot
runall(night)
time.sleep(60)

# night = '2018-04-18' # this night is done!
# runall(night)
# time.sleep(60)
# initdb(night)
# getshifts(night, n1 = 0, n2 = 100)
# dophot(night)

night = '2018-04-19'

night = '2018-04-20'

night = '2018-04-25' # rerunning to see if shifts are better (else always have the backup)
runall(night)
# time.sleep(60)

# night = '2018-04-26' # this night is done!
# runall(night)
# time.sleep(60)
# initdb(night)
# getshifts(night, n1 = 0, n2 = 10)
# dophot(night)

# night = '2018-04-30' # this night is done!
# runall(night)
# time.sleep(60)
#
# night = '2018-05-01' # this night is done!
# runall(night)
# time.sleep(60)
# savealigned(night)
# dophot(night)

# night = '2018-05-02' # this night is done!
# runall(night)
# time.sleep(60)
# getshifts(night)
# savealigned(night)
# try:
# dophot(night)
# except:
#     print night + " dophot failed"
#     print
# time.sleep(60)

# night = '2018-05-03' # this night is done!
# runall(night)
# time.sleep(60)
# getshifts(night, n1 = 100, n2 = 110)
# dophot(night)

# night = '2018-05-04' # this night is done!
# runall(night)
# time.sleep(60)
# getshifts(night, n1 = 0, n2 = 10)
# dophot(night)

# night = '2018-05-05' # this night is done!
# runall(night)
# time.sleep(60)
# initdb(night)
# getshifts(night)

# night = '2018-05-07' # this night is done!
# runall(night)
# time.sleep(60)
# dophot(night)


# night = '2018-05-09' # this night is done!
# runall(night)
# time.sleep(60)
# dophot(night)

night = '2018-05-10'   # redo ff/ds, has missing images???
# runall(night)
# time.sleep(60)
# dophot(night)

night = '2018-05-11'  # i can't get these to align at ALL wtf
# runall(night)
# time.sleep(60)
# getshifts(night, n1 = 2000, n2 = 2010)
# getshifts(night)
# savealigned(night)
# dophot(night)



gc.collect()
####################################################################3

# try:
#     initdb(night)  # initialize the df
#     print
#     getshifts(night)  # get shifts & save them to the df
#     # time.sleep(60)
#     print
#     savealigned(night)  # align & rebin individual ims, save to aligned folder
#     # time.sleep(60)
#     print
#     dophot(night) # do photometry & save values to the df
#     # time.sleep(60)
#     print
# except:
#     print
#     print night + " failed, check log for details"

# getshifts('2018-04-26', path='/media/amanda/demeter/maxi_j1820_070/')
# savealigned('2018-04-26', path='/media/amanda/demeter/maxi_j1820_070/')
# dophot('2018-04-26', path='/media/amanda/demeter/maxi_j1820_070/')


# # attemping to use 3/28 references for these nights -- this worked! -- it did not, in fact, work.
# night = '2018-04-13'
# runall(night) # done
#
# night = '2018-04-25'
# runall(night) # done

# night = '2018-03-28'
# ## makeref(night) #done
# ## getshifts(night) #done
# ## savealigned(night) #done
# # initdb(night)
# # dophot(night)
#
#
# night = '2018-03-29' # 5953 images
# # ## makeref(night,330,350) #redone with aij & more images
# # getshifts(night) # done; ~ 16 min
# # savealigned(night) # done; ~ 11.5 min
# # initdb(night) # done; ~ 11 seconds
# # dophot(night) # done; ~ 1 min
#
#
# night = '2018-03-30' #9668 images -- deleted a few b/c super bad weather meant bad shifts.
# # ## makeref(night,6010,6030) #redone with aij & more images
# # getshifts(night) # done; ~ 27 min !!
# ## but the aligned images are still saved so it should be okay
# # savealigned(night) # done; ~ 23 min
# ## can make a ref_stars_xy.txt file after you've got the first few aligned images saved
# ## which you need for dophot (done)
# # initdb(night) # done; ~ 25 s
# # dophot(night) # done; ~ 2 min -- need to not kill process for not calculating gaussian
# # also consider switching phot to grab fileIDs from images in the folder instead of fIDs in the df?
#
#
# night = '2018-05-05'
# # ## makeref(night,1000, 1015) #done
# # ## getshifts(night) #done
# # ## savealigned(night) #done
# # initdb(night)
# # dophot(night)
#
# night = '2018-04-06' #12449 images
# # makeref(night) # done w/ aij
# # getshifts(night) # done; ~36 min
# savealigned(night) #next
# # initdb(night)
# # dophot(night)