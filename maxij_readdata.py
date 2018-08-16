import pandas as pd
from matplotlib import pyplot as plt
# import numpy as np
import gc


path='/media/amanda/demeter/maxi_j1820_070/'
# path = './'
night = '2018-05-08'

#read in data as pandas dataframe:
# data = pd.read_pickle(path+night+'/photdata_'+night+'_shifts.pkl')
data = pd.read_pickle(path+night+'/data_'+night+'.pkl')
# print data.loc['16000']
print len(data)
#access data with column headers, here's a full list of the ones available:
# print list(data.columns.values)
## (shift data isn't actually saved in here yet so that's all going to be NaNs if you try to access it)

# for example, here I'm assigning the flux of each of these things from the dataframe to it's own list:
tyc_phot = data['phot_tyc']
maxij_phot = data['phot_maxij']
ref2_phot = data[ 'phot_ref2']
ref3_phot = data['phot_ref3']
ref4_phot = data[ 'phot_ref4']
ref5_phot = data['phot_ref5']
ref6_phot = data['phot_ref6']
#
time = data['os_time']
# time = data['fileID']
# # for i in range(len(tyc_phot)):
# #     print('tyc flux = %5.2f' % tyc_phot[i])
# #     print('MAXIJ flux = %5.2f' % maxij_phot[i])
#
# # plot lightcurves:
plt.plot(time,tyc_phot/5,'r.-')
# plt.plot(time,ref3_phot+ref4_phot+ref6_phot+ref5_phot,'r-.')
#
plt.plot(time,tyc_phot/3-50000,'m.-')
#
# # plt.plot(time,(ref2_phot),'r.-') # star at the bottom
# plt.plot(time,(ref3_phot-10000),'y.-')
# plt.plot(time,(ref4_phot-20000),'k.-')
# # plt.plot(time,(ref4_phot+ref3_phot)-10000,'.-m')
# plt.plot(time,(ref5_phot-30000),'c.-')
# # plt.plot(time,(ref6_phot-40000),'m.-') # star at the top
#
# plt.plot(time,(ref2_phot),'r.-') # star at the bottom
plt.plot(time,(ref3_phot),'y.-')
plt.plot(time,(ref4_phot),'k.-')
plt.plot(time,(ref5_phot),'c.-')
# plt.plot(time,(ref6_phot),'m.-') # star at the top

plt.plot(time,maxij_phot,'.-b')

# #
# # plt.plot(time,(maxij_phot/(ref3_phot+ref4_phot+ref5_phot)),'r.-')
plt.title(night)
# plt.show()
#
# plt.plot(time,maxij_phot/(ref4_phot+ref3_phot),'.')
plt.show()

gc.collect()