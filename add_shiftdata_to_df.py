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
fileID = [n.split('_')[2] for n in fnames]
# use maxij_initdb code to put shifts values into big db and then resave them

fileID = data['fileID']

print fileID

for results in p.imap_unordered(partial(assignshifts,
                                        shiftdata=shiftdata),
                                fileID, 25):

    fid = results[0]
    data.loc[fid, 'shift_x'] = results[1]
    data.loc[fid, 'shift_y'] = results[2]
    data.loc[fid, 'shift_corr_amplitude'] = results[3]


def assignshifts(fileID,shiftdata):
    shift_x = shiftdata.iloc[fileID,'xshift']
    shift_y = shiftdata.loc[fileID,'yshift']

    return [fileID,shift_x,shift_y,shift_corr_amplitude]