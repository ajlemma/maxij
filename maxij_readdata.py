import pandas as pd
from matplotlib import pyplot as plt
# import numpy as np


path='/media/amanda/demeter/maxi_j1820_070/'
# path = './'
night = '2018-05-05'

#read in data as pandas dataframe:
data = pd.read_pickle(path+night+'/photdata_'+night+'.pkl')
print data

#access data with column headers, here's a full list of the ones available:
print list(data.columns.values)
## (shift data isn't actually saved in here yet so that's all going to be NaNs if you try to access it)

# for example, here I'm assigning the flux of each of these things from the dataframe to it's own list:
tyc_phot = data['phot_tyc']
maxij_phot = data['phot_maxij']
ref2_phot = data[ 'phot_ref2']
ref3_phot = data['phot_ref3']
ref4_phot = data[ 'phot_ref4']
ref5_phot = data['phot_ref5']
ref6_phot = data['phot_ref6']

time = data['os_time']
# for i in range(len(tyc_phot)):
#     print('tyc flux = %5.2f' % tyc_phot[i])
#     print('MAXIJ flux = %5.2f' % maxij_phot[i])

# plot lightcurves:
plt.plot(time,maxij_phot,'.-')
# plt.plot(time,tyc_phot/5,'g-.')
# plt.plot(time,ref3_phot+ref4_phot+ref6_phot+ref5_phot,'r-.')

plt.plot(time,tyc_phot/9,'g.-')
plt.plot(time,ref3_phot+ref4_phot,'r.-')

plt.show()