import pandas as pd
from matplotlib import pyplot as plt
# import numpy as np


# path='/media/amanda/demeter/maxi_j1820_070/'
# night = '2018-05-05'

path = '/home/amanda/PycharmProjects/maxij/'
night = "test"

#read in data as pandas dataframe:
data = pd.read_pickle(path+night+'/photdata_'+night+'.pkl')


shiftdata = pd.read_csv(path+night+'/imshifts_test.txt',
                        header = None,
                        delimiter = r"\s+",
                        names=['filename', 'xshift', 'yshift','shift_amplitude'],)

# use maxij_initdb code to put shifts values into big db and then resave them

print shiftdata