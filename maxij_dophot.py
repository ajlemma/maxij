#######
''' maxij_dophot.py
by A Townsend
1.
'''
#######

import pandas as pd
from maxijdefs import *
from multiprocessing import Pool
from functools import partial


#######

def dophot(night, path='/media/amanda/demeter/maxi_j1820_070/'):
    time0 = timestart()
    print "Initializing database for " + night + "..."

    pathnam = path + night + '/'
    apathnam = pathnam + 'aligned/'  # folder for aligned images

    dataf = pd.read_pickle(pathnam + 'data_'+night+'.pkl')





    timefinish(time0)



if __name__ == "__main__":
    dophot('test', path='./')
