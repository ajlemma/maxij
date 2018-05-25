#####
from maxij_makeref import *
from maxij_getshifts import *
from maxij_savealigned import *
from maxij_initdb import *
from maxij_dophot import *
#####


night = '2018-03-28'
## makeref(night) #done
## getshifts(night) #done
## savealigned(night) #done
# initdb(night)
# dophot(night)


night = '2018-03-29' # 5953 images
# ## makeref(night,330,350) #redone with aij & more images
# getshifts(night) # done; ~ 16 min
# savealigned(night) # done; ~ 11.5 min
# initdb(night) # done; ~ 11 seconds
# dophot(night) # done; ~ 1 min


night = '2018-03-30' #9668 images -- deleted a few b/c super bad weather meant bad shifts.
# ## makeref(night,6010,6030) #redone with aij & more images
# getshifts(night) # done; ~ 27 min !!
## but the aligned images are still saved so it should be okay
# savealigned(night) # done; ~ 23 min
## can make a ref_stars_xy.txt file after you've got the first few aligned images saved
## which you need for dophot (done)
# initdb(night) # done; ~ 25 s
# dophot(night) # done; ~ 2 min -- need to not kill process for not calculating gaussian
# also consider switching phot to grab fileIDs from images in the folder instead of fIDs in the df?


night = '2018-05-05'
# ## makeref(night,1000, 1015) #done
# ## getshifts(night) #done
# ## savealigned(night) #done
# initdb(night)
# dophot(night)

night = '2018-04-06' #12449 images
# makeref(night) # done w/ aij
# getshifts(night) # done; ~36 min
savealigned(night) #next
# initdb(night)
# dophot(night)