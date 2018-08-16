#######
''' maxij_rhonicerxcorr.py
by A Townsend
'''
#######

import numpy as np

import pandas as pd
import glob
import time
from natsort import natsorted

from maxij_nicerdefs import *

from scipy.signal import fftconvolve


import matplotlib
import matplotlib.pyplot as plt
from pylab import subplot, subplots_adjust

#######


## user inputs:

# info for reading x-ray data (directory and filter for night)
# for example, for 3/28 the NICER date is 113
xpathnam = 'delivery/'
filename_filter = 'ni1200120113_0mpu7_silver_GTI*.lc' # for 3/28

# infor to read rho optical data (directory and night)
opathnam = './maxij_data_copies/'
night = '2018-03-28'

# computer clock is ~51.5s fast (per lab tests + SSE calculation)
# (((POSSIBLY 8.5 s slow, but this seems to work?)))
# clocklag = -51.5
clocklag = 8.5
# clocklag -=60*3
##################
# get filenames for night's xray data
ftest = xpathnam + filename_filter
flistD = glob.glob(ftest)

# all files
flist = [fnam.split('/')[1] for fnam in flistD]
flist = natsorted(flist)
print len(flist)
ndat_t = []
ndat_cts = []
ndat_cts1 = []
ndat_cts2 = []
ndat_cts34 = []
ndat_fnam = []

for fileX in flist:
    t,cts=rd_nicer_lc(xpathnam,[fileX])
    t1,cts1=rd_nicer_lc1(xpathnam,[fileX])
    t2,cts2=rd_nicer_lc2(xpathnam,[fileX])
    t34,cts34=rd_nicer_lc34(xpathnam,[fileX])
    filename_column = [fileX for x in np.arange(len(t))]

    ndat_t.append(t)
    ndat_cts.append(cts)
    ndat_cts1.append(cts1)
    ndat_cts2.append(cts2)
    ndat_cts34.append(cts34)
    ndat_fnam.append(filename_column)

ndict = {'filename': [item for sublist in ndat_fnam for item in sublist],
          'time': [item for sublist in ndat_t for item in sublist],
          'cts': [item for sublist in ndat_cts for item in sublist],
          'cts1': [item for sublist in ndat_cts1 for item in sublist],
          'cts2': [item for sublist in ndat_cts2 for item in sublist],
          'cts34': [item for sublist in ndat_cts34 for item in sublist]
          }

nicerdata = pd.DataFrame(ndict, index=ndict['time'])
# print nicerdata

#read in optical data as pandas dataframe:
rhodata = pd.read_pickle(opathnam+'data_'+night+'.pkl')
# print data
print len(rhodata)
rhodata['fixed_time'] = rhodata['os_time'] + clocklag


# do not change this! it is absolutely calibrated to the 28th.
# unix epoch timestamp - nicer absolute time
# both at 2018-03-28 08:05:21.000 UTC
dt = 1522224321-133689923
print dt
nicerdata['unix_time'] = nicerdata['time'] + dt

t_startrho = np.min(rhodata['fixed_time'])
t_stoprho = np.max(rhodata['fixed_time'])
flist_overlap = []
# for each GTI:
i = 1
for fnam in flist:
    t = nicerdata.loc[nicerdata['filename']==fnam]['time']
    t_start = np.min(t) #t[0]
    t_start_unix = np.int(t_start + dt) #round down to the second

    t_stop = np.max(t)
    t_stop_unix = np.int(t_stop + dt)

    print "Time at beginning of Nicer GTI " + str(i)
    print str(t_start_unix) + " <-- unix epoch timestamp"
    print time.strftime('%m/%d/%Y %H:%M:%S',  time.gmtime(t_start_unix)) + " UTC"


    print "Time at end of GTI: " + time.strftime('%m/%d/%Y %H:%M:%S', time.gmtime(t_stop_unix)) + " UTC"

    print "Length of GTI: " + str(t_stop_unix - t_start_unix) + " seconds"

    if (t_start_unix >= t_startrho) and (t_stop_unix <= t_stoprho):
        print "GTI Overlaps with RHO data"
        flist_overlap.append(fnam)


    print
    i += 1

mar28a_nfiles = flist[0:3]
# mar28a_nfiles = flist[3:6] #mar28b


mar28a_nicerdata = nicerdata.loc[nicerdata['filename'].isin(mar28a_nfiles)]
mar28a_rhodata = rhodata.loc[(rhodata['fixed_time'] > np.floor(np.min(mar28a_nicerdata['unix_time']))) & (rhodata['fixed_time'] < np.ceil(np.max(mar28a_nicerdata['unix_time'])))]
print np.max(mar28a_nicerdata['unix_time']) - np.min(mar28a_nicerdata['unix_time'])
print np.max(mar28a_rhodata['fixed_time']) - np.min(mar28a_rhodata['fixed_time'])

xt = mar28a_nicerdata['unix_time']
xlc = mar28a_nicerdata['cts']
ot = mar28a_rhodata['fixed_time']
olc = mar28a_rhodata['phot_maxij']
tlc = mar28a_rhodata['phot_tyc']





# ## full lightcurve plot for data part
# fig1,ax=plt.subplots(figsize=[15,10],facecolor='w')
# ax.plot(xt,xlc-np.median(xlc),'.-b')
# ax.plot(ot,olc-np.median(olc),'.-r')
# ax.plot(ot,tlc-np.median(tlc),'.-k',alpha=.3)
# # ax.plot(xt,xlc,'.-b')
# # ax.plot(ot,olc,'.-r')
# # ax.plot(ot,tlc,'.-k',alpha=.3)
# # ax.set_ylim(0,40000)
# ax.legend(['MaxiJ@Nicer','MaxiJ@RHO',"Tycho@RHO"])
# ax.set_xlabel('Time (unix epoch, seconds)')
# ax.set_ylabel('Flux, arbitrary')
# # plt.show()

t = np.arange(np.floor(np.min(mar28a_nicerdata['unix_time'])),np.ceil(np.max(mar28a_nicerdata['unix_time'])))
tyc_min = 300000.
xts = make_1s_ts(t, np.array(xt),np.array(xlc)-np.median(xlc))
ots = make_1s_ts(t, np.array(ot[tlc>tyc_min]),np.array(olc[tlc>tyc_min])-np.median(olc[tlc>tyc_min]))
tts = make_1s_ts(t, np.array(ot[tlc>tyc_min]),np.array(tlc[tlc>tyc_min])-np.median(tlc[tlc>tyc_min]))
print len(t), len(ots)


ntot = 128
xts_c = chunks(xts,ntot)
ots_c = chunks(ots,ntot)
tts_c = chunks(tts,ntot)
tc = np.arange(ntot)
# print xts_c
print "number of 128s chunks: " + str(len(xts_c))
ctot = len(xts_c)

# # plot of time series over chunk of data
# fig2,ax1=plt.subplots(figsize=[15,15],facecolor='w')
# subplots_adjust(hspace=0.000)
# for c,v in enumerate(xrange(ctot)):
#     v = v+1
#     ax1 = subplot(ctot,1,v)
#     ax1.plot(tc, xts_c[c],'.-b')
#     ax1.plot(tc, ots_c[c],'.-r')
#     ax1.plot(tc, tts_c[c],'.-k',alpha=.3)
# ax1.legend(['MaxiJ@Nicer','MaxiJ@RHO',"Tycho@RHO"])

# plt.show()
## to get power spectra:
# tf=np.arange(len(t))
# fx,px,pxtot = pds_chunk128(tf,xts,0,np.max(tf))
# fo,po,potot = pds_chunk128(tf,ots,0,np.max(tf))

# to get cross corrl'ns:
auto_x = []
auto_o = []
auto_t = []

crosscorr_xo = []
crosscorr_xt = []
crosscorr_ot = []

for c in xrange(ctot):
    auto_x.append(fftconvolve(xts_c[c],xts_c[c][::-1],mode='same')) #x-ray autocorrelation
    auto_o.append(fftconvolve(ots_c[c],ots_c[c][::-1],mode='same')) #optical autocorrelation
    auto_t.append(fftconvolve(tts_c[c],tts_c[c][::-1],mode='same')) #tycho reference star autocorrelation

    crosscorr_xo.append(fftconvolve(xts_c[c],ots_c[c][::-1],mode='same')) # cross correlation
    crosscorr_xt.append(fftconvolve(xts_c[c],tts_c[c][::-1],mode='same')) # cross corr'ln w/ refstar
    crosscorr_ot.append(fftconvolve(ots_c[c],tts_c[c][::-1],mode='same'))

# #plot auto correlations
# fig3,ax=plt.subplots(figsize=[9,6])
# ax.plot(auto_x,'b')
# ax.plot(auto_o,'r')
# ax.plot(auto_t,'k')
# # ax.plot([63.,63.],[0,3.e29],'y',alpha=.7)
# # ax.set_ylim(0,9e24)
# # ax.set_xlim(40,80)
# # print np.argmax(auto_x)
print np.argmax(auto_o)
# plt.show()

# fig4,ax=plt.subplots(figsize=[15,10],facecolor='w')
# subplots_adjust(hspace=0.000)
# for c,v in enumerate(xrange(ctot)):
#     v = v+1
#     ax = subplot(ctot,1,v)
#     ax.plot(crosscorr_xo[c],'g')
#     ax.plot(crosscorr_xt[c],'c')
#     ax.plot(crosscorr_ot[c], 'm',alpha=.3)
# # ax.plot(crosscorr_xo - crosscorr_xt,'k')
# # ax.plot(auto_o/4,'r')
# # ax.plot(auto_x,'b')
# # ax.plot(auto_t/3,'k')
#     ax.plot([np.argmax(auto_x),np.argmax(auto_x)],[-4.e9,2.e10],'y',alpha=.7)
#     ax.set_ylim(-2.5e9,2.5e9)
#     ax.set_xlim(0,128)
# ax.set_title('Optical+X-ray Correlation (one chunk of 128-second overlapping data) 51.5s lead')
# ax.legend(['Cross-Correlation (MaxiJ Nicer/MaxiJ RHO)',
#            'Cross-Correlation (MaxiJ Nicer/Tycho)',
#            'Cross-Correlation (MaxiJ RHO/Tycho)',
#            # 'RHO Optical Auto-correlation',
#            # 'Nicer X-ray Auto-correlation',
#            # 'RHO Tycho Reference Auto-correlation',
#            'Autocorrelation Peak'])
print np.argmax(crosscorr_xo)




# fig5,ax=plt.subplots(figsize=[15,10],facecolor='w')
# ax.plot(np.mean(crosscorr_xo,axis=0), 'g')
# ax.plot(np.mean(crosscorr_xt,axis=0), 'c')
# ax.plot(np.mean(crosscorr_ot,axis=0), 'm', alpha=.3)
# ax.plot([np.argmax(auto_x), np.argmax(auto_x)], [-4.e9, 2.e10], 'y', alpha=.7)
# ax.set_ylim(np.min(np.mean(crosscorr_xo,axis=0)), np.max(np.mean(crosscorr_xo,axis=0)))
# ax.set_xlim(0,128)

# plt.show()
# raw_input()

ccxoa = crosscorr_xo
ccxta = crosscorr_xt
ccota = crosscorr_ot

# mar28a_nfiles = flist[0:3]
mar28a_nfiles = flist[3:6] #mar28b


mar28a_nicerdata = nicerdata.loc[nicerdata['filename'].isin(mar28a_nfiles)]
mar28a_rhodata = rhodata.loc[(rhodata['fixed_time'] > np.floor(np.min(mar28a_nicerdata['unix_time']))) & (rhodata['fixed_time'] < np.ceil(np.max(mar28a_nicerdata['unix_time'])))]
print np.max(mar28a_nicerdata['unix_time']) - np.min(mar28a_nicerdata['unix_time'])
print np.max(mar28a_rhodata['fixed_time']) - np.min(mar28a_rhodata['fixed_time'])

xt = mar28a_nicerdata['unix_time']
xlc = mar28a_nicerdata['cts']
ot = mar28a_rhodata['fixed_time']
olc = mar28a_rhodata['phot_maxij']
tlc = mar28a_rhodata['phot_tyc']





# ## full lightcurve plot for data part
# fig1,ax=plt.subplots(figsize=[15,10],facecolor='w')
# ax.plot(xt,xlc-np.median(xlc),'.-b')
# ax.plot(ot,olc-np.median(olc),'.-r')
# ax.plot(ot,tlc-np.median(tlc),'.-k',alpha=.3)
# # ax.plot(xt,xlc,'.-b')
# # ax.plot(ot,olc,'.-r')
# # ax.plot(ot,tlc,'.-k',alpha=.3)
# # ax.set_ylim(0,40000)
# ax.legend(['MaxiJ@Nicer','MaxiJ@RHO',"Tycho@RHO"])
# ax.set_xlabel('Time (unix epoch, seconds)')
# ax.set_ylabel('Flux, arbitrary')
# # plt.show()

t = np.arange(np.floor(np.min(mar28a_nicerdata['unix_time'])),np.ceil(np.max(mar28a_nicerdata['unix_time'])))
tyc_min = 300000.
xts = make_1s_ts(t, np.array(xt),np.array(xlc)-np.median(xlc))
ots = make_1s_ts(t, np.array(ot[tlc>tyc_min]),np.array(olc[tlc>tyc_min])-np.median(olc[tlc>tyc_min]))
tts = make_1s_ts(t, np.array(ot[tlc>tyc_min]),np.array(tlc[tlc>tyc_min])-np.median(tlc[tlc>tyc_min]))
print len(t), len(ots)


ntot = 128
xts_c = chunks(xts,ntot)
ots_c = chunks(ots,ntot)
tts_c = chunks(tts,ntot)
tc = np.arange(ntot)
# print xts_c
print "number of 128s chunks: " + str(len(xts_c))
ctot = len(xts_c)

# # plot of time series over chunk of data
# fig2,ax1=plt.subplots(figsize=[15,15],facecolor='w')
# subplots_adjust(hspace=0.000)
# for c,v in enumerate(xrange(ctot)):
#     v = v+1
#     ax1 = subplot(ctot,1,v)
#     ax1.plot(tc, xts_c[c],'.-b')
#     ax1.plot(tc, ots_c[c],'.-r')
#     ax1.plot(tc, tts_c[c],'.-k',alpha=.3)
# ax1.legend(['MaxiJ@Nicer','MaxiJ@RHO',"Tycho@RHO"])

# plt.show()
## to get power spectra:
# tf=np.arange(len(t))
# fx,px,pxtot = pds_chunk128(tf,xts,0,np.max(tf))
# fo,po,potot = pds_chunk128(tf,ots,0,np.max(tf))

# to get cross corrl'ns:
auto_x = []
auto_o = []
auto_t = []

crosscorr_xo = []
crosscorr_xt = []
crosscorr_ot = []

for c in xrange(ctot):
    auto_x.append(fftconvolve(xts_c[c],xts_c[c][::-1],mode='same')) #x-ray autocorrelation
    auto_o.append(fftconvolve(ots_c[c],ots_c[c][::-1],mode='same')) #optical autocorrelation
    auto_t.append(fftconvolve(tts_c[c],tts_c[c][::-1],mode='same')) #tycho reference star autocorrelation

    crosscorr_xo.append(fftconvolve(xts_c[c],ots_c[c][::-1],mode='same')) # cross correlation
    crosscorr_xt.append(fftconvolve(xts_c[c],tts_c[c][::-1],mode='same')) # cross corr'ln w/ refstar
    crosscorr_ot.append(fftconvolve(ots_c[c],tts_c[c][::-1],mode='same'))

# #plot auto correlations
# fig3,ax=plt.subplots(figsize=[9,6])
# ax.plot(auto_x,'b')
# ax.plot(auto_o,'r')
# ax.plot(auto_t,'k')
# # ax.plot([63.,63.],[0,3.e29],'y',alpha=.7)
# # ax.set_ylim(0,9e24)
# # ax.set_xlim(40,80)
# # print np.argmax(auto_x)
print np.argmax(auto_o)
# plt.show()

# fig4,ax=plt.subplots(figsize=[15,10],facecolor='w')
# subplots_adjust(hspace=0.000)
# for c,v in enumerate(xrange(ctot)):
#     v = v+1
#     ax = subplot(ctot,1,v)
#     ax.plot(crosscorr_xo[c],'g')
#     ax.plot(crosscorr_xt[c],'c')
#     ax.plot(crosscorr_ot[c], 'm',alpha=.3)
# # ax.plot(crosscorr_xo - crosscorr_xt,'k')
# # ax.plot(auto_o/4,'r')
# # ax.plot(auto_x,'b')
# # ax.plot(auto_t/3,'k')
#     ax.plot([np.argmax(auto_x),np.argmax(auto_x)],[-4.e9,2.e10],'y',alpha=.7)
#     ax.set_ylim(-2.5e9,2.5e9)
#     ax.set_xlim(0,128)
# ax.set_title('Optical+X-ray Correlation (one chunk of 128-second overlapping data) 51.5s lead')
# ax.legend(['Cross-Correlation (MaxiJ Nicer/MaxiJ RHO)',
#            'Cross-Correlation (MaxiJ Nicer/Tycho)',
#            'Cross-Correlation (MaxiJ RHO/Tycho)',
#            # 'RHO Optical Auto-correlation',
#            # 'Nicer X-ray Auto-correlation',
#            # 'RHO Tycho Reference Auto-correlation',
#            'Autocorrelation Peak'])
print np.argmax(crosscorr_xo)




# fig5,ax=plt.subplots(figsize=[15,10],facecolor='w')
# ax.plot(np.mean(crosscorr_xo,axis=0), 'g')
# ax.plot(np.mean(crosscorr_xt,axis=0), 'c')
# ax.plot(np.mean(crosscorr_ot,axis=0), 'm', alpha=.3)
# ax.plot([np.argmax(auto_x), np.argmax(auto_x)], [-4.e9, 2.e10], 'y', alpha=.7)
# ax.set_ylim(np.min(np.mean(crosscorr_xo,axis=0)), np.max(np.mean(crosscorr_xo,axis=0)))
# ax.set_xlim(0,128)

ccxob = crosscorr_xo
ccxtb = crosscorr_xt
ccotb = crosscorr_ot

ccxo = ccxoa + ccxob
ccxt = ccxta + ccxtb
ccot = ccota + ccotb

fig6,ax=plt.subplots(figsize=[15,10],facecolor='w')
ax.plot(np.mean(ccxo,axis=0), 'g')
ax.plot(np.mean(ccxt,axis=0), 'c')
ax.plot(np.mean(ccot,axis=0), 'm', alpha=.3)
ax.plot([np.argmax(auto_x), np.argmax(auto_x)], [-4.e9, 2.e10], 'y', alpha=.7)
ax.set_ylim(np.min(np.mean(crosscorr_xo,axis=0)), np.max(np.mean(crosscorr_xo,axis=0)))
ax.set_xlim(0,128)
# ax.legend(['Cross-Correlation (MaxiJ Nicer/MaxiJ RHO)',
#            'Cross-Correlation (MaxiJ Nicer/Tycho)',
#            'Cross-Correlation (MaxiJ RHO/Tycho)',
#            # 'RHO Optical Auto-correlation',
#            # 'Nicer X-ray Auto-correlation',
#            # 'RHO Tycho Reference Auto-correlation',
#            'Autocorrelation Peak'])



plt.show()
raw_input()

