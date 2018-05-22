#######

from maxijdefs import *
import pandas as pd
import numpy as np

#######
night = '2018-05-16'
pathnam = '/media/amanda/demeter/maxi_j1820_070/' + night + '/'

time0 = timestart()

print "getting os timestamps..."
ostime = pd.read_csv(pathnam + 'timestamps_' + night + '.txt',
                     skiprows=2,
                     header=None,
                     names=['filename1','filename2','epochtime(s)'],
                     delimiter=r"\s+")
# print len(ostime)
ostime = ostime[ostime['filename1'].str.contains('LED_1s')]
fileID_os = [n.split('_')[-1] for n in ostime['filename1']]
ostime['fileID'] = fileID_os
ostime = ostime.set_index('fileID')

fileID_os = natsorted(fileID_os)
print len(fileID_os)



print "getting LED photometry measurements..."
led_df = pd.read_csv(pathnam + 'Measurements.txt',
                      skiprows=1,
                      header=None,
                      names=['filename1','filename2','junk1','junk2','junk3','led10s','ntp1min'],
                     delimiter=r"\s+")
# ledphot['filename1'] = ledphot['filename1']+'_'
# ledphot['filename'] = ledphot[['filename1', 'filename2']].apply(lambda x: ''.join(x), axis=1)
del led_df['junk1']
del led_df['junk2']
del led_df['junk3']
fileID = [n.split('_')[-1] for n in led_df['filename1']]
led_df['fileID'] = fileID
led_df = led_df.set_index('fileID')
fileID = natsorted(fileID)


led_df['filetime'] = np.nan
led_df['filetime_s'] = np.nan
led_df['os_time'] = np.nan

def assign_times(id):
    led_df.loc[id, 'filetime'] = parse_time(led_df.loc[id,'filename2'])
    led_df.loc[id,'filetime_s'] = get_sec(led_df.loc[id,'filetime'])
    led_df.loc[id,'os_time'] = ostime.loc[id,'epochtime(s)']

for id in fileID:
    assign_times(id)

# print led_df[0:10]

led_df.to_pickle(pathnam+'timedata.pkl')


timefinish(time0)