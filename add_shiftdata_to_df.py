import pandas as pd
from matplotlib import pyplot as plt
from multiprocessing import Pool
from functools import partial

# import numpy as np



path='/media/amanda/demeter/maxi_j1820_070/'
night = '2018-03-30'

# path = '/home/amanda/PycharmProjects/maxij/'
# night = "test"

pathnam = path + night

def assignshifts(fname,shiftdata):
    fileID = fname.split('_')[2]
    shift_x = shiftdata.loc[fname,'xshift']
    shift_y = shiftdata.loc[fname,'yshift']
    shift_corr_amplitude = shiftdata.loc[fname,'shift_amplitude']

    return [fileID,shift_x,shift_y,shift_corr_amplitude]


#read in data as pandas dataframe:
data = pd.read_pickle(path+night+'/photdata_'+night+'.pkl')


shiftdata = pd.read_csv(path+night+'/imshifts_' + night + '.txt',
                        header = None,
                        delimiter = r"\s+",
                        names=['filename', 'xshift', 'yshift','shift_amplitude'],)

# use maxij_initdb code to put shifts values into big db and then resave them

filenames = shiftdata['filename']
shiftdata = shiftdata.set_index('filename')


p = Pool(6)
for results in p.imap_unordered(partial(assignshifts,
                                        shiftdata=shiftdata),
                                filenames, 25):
    # print results
    fid = results[0]
    data.loc[fid, 'shift_x'] = results[1]
    data.loc[fid, 'shift_y'] = results[2]
    data.loc[fid, 'shift_corr_amplitude'] = results[3]

# print data
data.to_pickle(pathnam + '/photdata_' + night + '_shifts.pkl')