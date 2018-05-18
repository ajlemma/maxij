#######

import pandas as pd
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
#######
night = '2018-05-16'
pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'

led_df=pd.read_pickle(pathnam+'timedata.pkl')

# print led_df.iloc[0]['os_time']-led_df.iloc[0]['filetime_s']
# plt.plot(led_df['os_time']-1526443200.,led_df['ntp1min'],'.-')
# plt.plot(led_df['filetime_s'],led_df['ntp1min'],'r.-')
# plt.show()

bright = led_df[led_df['ntp1min']>1000000]
# plt.plot(bright['os_time']-1526443200.-bright['filetime_s'],'.-')
# plt.show()

print led_df.iloc[0]

# tosave = led_df.filter(['filetime', 'filetime_s', 'os_time', 'led10s', 'ntp1min'], axis=1)
# print tosave[0:10]
# tosave.to_csv(pathnam+'ntp_data.txt',sep=' ')

acorr = fftconvolve(led_df['ntp1min'],led_df['ntp1min'][::-1])

plt.plot(acorr)
plt.show()