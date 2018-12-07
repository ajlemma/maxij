#######
''' maxij_overlap_analysis_28nov2018.py
by A Townsend

note: make sure to edit the "deletethese" variable for each night (will have to run multiple times)
'''
#######

## variables that need to be edited for each night:
# deletethese = [] # start with this
deletethese = [6, 10, 14] # for 3/28

#######
import pandas as pd
# import numpy as np
from maxijdefs import *
from maxij_nicerdefs import *
# import glob
from itertools import groupby
# from datetime import datetime
from scipy.signal import fftconvolve
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import subplot, subplots_adjust, legend

#######

# #######



#######
time0,msg = timestart()
print msg

#######

#read in timeseries data as pandas dataframe:
print "reading in timeseries data as pandas df"
tsdata = pd.read_pickle('./timeseries_data_all.pkl')
# print data
print "there are " + str(len(tsdata)) + "nights of data:"
onights = tsdata['nights'].tolist()
print onights
night = onights[0]

#######
print
print "pulling out time and photometry datasets and creating 1-s timeseries"

## define timeframe based on optical data
## make 1-s timeseries' for whole night (with zeroes where no data)
times = np.asarray(tsdata.loc[night,'unix_epoch_time'])
xts = np.asarray(tsdata.loc[night, 'xray_maxij_ts'])
ots = np.asarray(tsdata.loc[night, 'rho_maxij_ts'])
tts = np.asarray(tsdata.loc[night, 'rho_tycho_ts'])

fig1,ax = plt.subplots(figsize=[8,12],facecolor='w')
fig1.suptitle('Raw Timeseries Data for ' + night)
ax2 = subplot(3,1,2)
ax2.plot(times[ots>0],ots[ots>0],'r')
ax2.set_title('RHO MAXI J1820+070 Data')
ax3 = subplot(3,1,3,sharex=ax2)
ax3.plot(times[tts>0],tts[tts>0],'k',alpha=.5)
ax3.set_title("RHO Reference Star Data")
ax1 = subplot(3,1,1,sharex=ax2)
ax1.plot(times[xts>0],xts[xts>0],'b')
ax1.set_title('X-ray MAXI J1820+070 Data')
# plt.show()
plt.savefig('overlap_analysis_plots/'+night+'_all_timeseries.pdf')

#######
print
print "splitting data into 64 & 128-s chunks to get power spectra..."
print '(ignore the  "Mean of empty slice" warning messages if they show up)'

## for all data
## split into 64-s chunks
ntot = 64
xts_c = chunks(xts, ntot)
ots_c = chunks(ots, ntot)
tts_c = chunks(tts, ntot)
time_c = chunks(times, ntot)
tc = np.arange(ntot)
ctot = len(ots_c)



datadict = {}
# overlap = []
for c in xrange(ctot):
    # cm = median-subtracted cleaned timeseries
    xts_cm = np.zeros(np.shape(xts_c[c]))
    ots_cm = np.zeros(np.shape(ots_c[c]))
    tts_cm = np.zeros(np.shape(tts_c[c]))

    # m is the median of the timeseries
    xts_m = np.median(xts_c[c][xts_c[c] > 0])
    ots_m = np.median(ots_c[c][ots_c[c] > 0])
    tts_m = np.median(tts_c[c][tts_c[c] > 0])

    # fixed = set zeroes equal to median
    xts_f = xts_cm + xts_m
    ots_f = ots_cm + ots_m
    tts_f = tts_cm + tts_m

    for i in xrange(ntot):
        if xts_c[c][i] > 0:
            xts_cm[i] = xts_c[c][i] - xts_m
            xts_f[i] = xts_c[c][i]
        if ots_c[c][i] > 0:
            ots_cm[i] = ots_c[c][i] - ots_m
            ots_f[i] = ots_c[c][i]
        if tts_c[c][i] > 0:
            tts_cm[i] = tts_c[c][i] - tts_m
            tts_f[i] = tts_c[c][i]

    chunkID = np.arange(ctot)
    chunkdict = {'unix_time': time_c[c],
                 'chunktime': tc,
                 'xray_maxij_clean_ts': xts_c[c],
                 'rho_maxij_clean_ts': ots_c[c],
                 'rho_tycho_clean_ts': tts_c[c],

                 'xray_maxij_fixed_ts': xts_f,
                 'rho_maxij_fixed_ts': ots_f,
                 'rho_tycho_fixed_ts': tts_f,

                 'xray_maxij_med_ts': xts_cm,
                 'rho_maxij_med_ts': ots_cm,
                 'rho_tycho_med_ts': tts_cm
                 }

    chunkframe = pd.DataFrame(chunkdict, index=tc)
    datadict[str(c)] = chunkframe


def powerspectra(pdkey):

    c128 = np.arange(0,ctot,2) # 0 to n by 2

    powerspec128 = []
    powerspec64 = []
    for c in c128:

        freq64, power = pds_single_chunk(datadict[str(c)][pdkey])
        powerspec64.append(power)
        freq64, power = pds_single_chunk(datadict[str(c + 1)][pdkey])
        powerspec64.append(power)

        fluxes_rho_maxij = np.ravel([list(datadict[str(c)][pdkey]),
                                     list(datadict[str(c + 1)][pdkey])])

        c128_times = np.ravel([list(datadict[str(c)]['unix_time']),
                               list(datadict[str(c + 1)]['unix_time'])])


    #
        freq128, power = pds_single_chunk(fluxes_rho_maxij)
        powerspec128.append(power)

    powerspec64 = np.asarray(powerspec64)
    powerspec128 = np.asarray(powerspec128)
    return freq64,powerspec64, freq128, powerspec128

print "calculating power spectrum for each chunk..."
freq64, rho_maxij_pspec64, freq128, rho_maxij_pspec128 = powerspectra('rho_maxij_fixed_ts')
freq64, rho_tycho_pspec64, freq128, rho_tycho_pspec128 = powerspectra('rho_tycho_fixed_ts')
freq64, xray_maxij_pspec64, freq128, xray_maxij_pspec128 = powerspectra('xray_maxij_fixed_ts')

print np.shape(rho_maxij_pspec128)

rho_maxij_ptot = np.sum(rho_maxij_pspec128, axis=0)
rho_maxij_ptotnorm = rho_maxij_ptot/(rho_maxij_ptot[0]**0.5)
rho_maxij_flog,rho_maxij_ptotlog = logbin(freq128,rho_maxij_ptotnorm,0.05)

rho_tycho_ptot = np.sum(rho_tycho_pspec128, axis=0)
rho_tycho_ptotnorm = rho_tycho_ptot/(rho_tycho_ptot[0]**0.5)
rho_tycho_flog,rho_tycho_ptotlog = logbin(freq128,rho_tycho_ptotnorm,0.05)

xray_maxij_ptot = np.sum(np.nan_to_num(xray_maxij_pspec128), axis=0)
xray_maxij_ptotnorm = xray_maxij_ptot/(xray_maxij_ptot[0]**0.5)
xray_maxij_flog,xray_maxij_ptotlog = logbin(freq128,xray_maxij_ptotnorm,0.05)


print "plotting sum of power spectra, 128s chunks..."
fig2,ax = plt.subplots(facecolor='w',figsize=[10, 12])
ax3 = subplot(2,1,1)
ax3.set_title('Sum over all 128s Power Spectra for ' + night)
# ax3.plot(10.**rho_maxij_flog[rho_maxij_ptotlog>0.], (10.**rho_maxij_flog[rho_maxij_ptotlog>0.])*rho_maxij_ptotlog[rho_maxij_ptotlog>0.],
#          '-r',drawstyle='steps-mid',lw=3)
# ax3.plot(10.**rho_tycho_flog[rho_tycho_ptotlog>0.], (10.**rho_tycho_flog[rho_tycho_ptotlog>0.])*rho_tycho_ptotlog[rho_tycho_ptotlog>0.],
#          '-k',drawstyle='steps-mid',lw=3,alpha=.5)
ax3.plot(10.**xray_maxij_flog[xray_maxij_ptotlog>0.], (10.**xray_maxij_flog[xray_maxij_ptotlog>0.])*xray_maxij_ptotlog[xray_maxij_ptotlog>0.],
         '-b',drawstyle='steps-mid',lw=3)
ax3.plot(10.**rho_maxij_flog[rho_maxij_ptotlog>0.],(10.**rho_maxij_flog[rho_maxij_ptotlog>0.])*rho_maxij_ptotlog[rho_maxij_ptotlog>0.] - (10.**rho_tycho_flog[rho_tycho_ptotlog>0.])*rho_tycho_ptotlog[rho_tycho_ptotlog>0.],
         '-r',drawstyle='steps-mid',lw=3)
ax3.set_xlabel('Frequency (Hz)')
ax3.set_ylabel('f x Power')
ax3.axvline(x=0.09, c='c', linestyle='--', lw=2)
ax3.axvline(x=0.045, c='m', linestyle='--', lw=2)
# ax3.axvline(x=0.0225,color='y',linestyle='--',lw =2)
ax3.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax3.set_xscale('log')
ax3.set_xlim(7e-3,0.5)
ax3.set_ylim(0,1.3e3)
plt.legend([
            # 'RHO Tycho',
            'NICER (X-ray) Maxi J1820+070',
            'RHO (Optical) Maxi J1820+070',
            # 'RHO MaxiJ-Tycho',
            '11s' ,
            '22s' ,
            # '44s',
            ],
          loc='upper center', bbox_to_anchor=(.5, -0.15),
          fancybox=True, shadow=True, ncol=2)
# plt.show()
plt.savefig('overlap_analysis_plots/'+night+'_power_spectra_sum.pdf')

print "creating 2D power spectrum"
logfreq128, logpsd_rho_maxij128 = p_2_logp(freq128, rho_maxij_pspec128,0.03,4)
rho_maxij_map128 = np.asarray((10.**logfreq128[:1])*(logpsd_rho_maxij128[:,1:]))
print len(logfreq128)
print 10**logfreq128

fig99,ax = plt.subplots(facecolor='w')
ax.plot(10**logfreq128, np.arange(len(logfreq128)),'.-')
# ax.plot(np.arange(len(logfreq128)),np.log10())
plt.savefig('overlap_analysis_plots/testfig.pdf')

print np.shape(rho_maxij_map128)
print "plotting median-subtracted & cleaned up timeseries"
time_clean = np.ravel([list(datadict[str(c)]['unix_time']) for c in np.arange(ctot)])
xray_maxij_clean = np.ravel([list(datadict[str(c)]['xray_maxij_med_ts']) for c in np.arange(ctot)])
rho_maxij_clean = np.ravel([list(datadict[str(c)]['rho_maxij_med_ts']) for c in np.arange(ctot)])
rho_tycho_clean = np.ravel([list(datadict[str(c)]['rho_tycho_med_ts']) for c in np.arange(ctot)])

fig3,ax = plt.subplots(facecolor='w')

ax1 = subplot(3,1,1)
ax1.plot(time_clean[(time_clean>0) & (xray_maxij_clean!=0)],xray_maxij_clean[(time_clean>0) & (xray_maxij_clean!=0)],'b-')

ax2 = subplot(3,1,2)
ax2.plot(time_clean[time_clean>0],rho_maxij_clean[time_clean>0],'r-')

ax3 = subplot(3,1,3)
ax3.imshow(rho_maxij_map128.T, vmax = 40, cmap='magma',interpolation='none', origin='lower')
# this one needs a gridspec to look nice though so do that later; add in the colorbar too
# plt.show()


fig4=plt.figure(figsize=(6,12), facecolor='w')
fig4.suptitle('RHO Optical MAXI J1820+070\n 2D Power Spectrum (128s) \n' + night, fontsize=16)

grd = gridspec.GridSpec(1,2, wspace = .25, hspace = 0.0, width_ratios=(10, 1))


ax2 = fig4.add_subplot(grd[0,0])
# ax2.set_axis_off()
im = ax2.imshow(rho_maxij_map128, vmax = 40, cmap='magma',interpolation='none', extent = [10**logfreq128[0],10**logfreq128[-1],1,0])
ax2.autoscale(False)
# ax2.set_xscale('log', basex=10)
# ax2.set_xticks([.01, .1])
# ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.get_xaxis().set_tick_params(which='both',direction='out', width=1)

# ax1 = ax2.twiny()
ax2.axvline(x=10**np.log10(0.09), c='c', linestyle='--', lw=2)
ax2.axvline(x=10**np.log10(0.045), c='r', linestyle='--', lw=2)
# ax1.autoscale(False)

cbar_ax = fig4.add_subplot(grd[0,1])
plt.colorbar(im,cbar_ax);
# plt.plot([4.1,4.1],[0,86],c='r',linestyle='--')
# ax.plot([10.25,10.25],[0,86],c='r',linestyle='--')
# ax.set_title('RHO MAXI J1820+070 - ' +night, y=1.015)
# plt.colorbar()


# plt.show()
plt.savefig('overlap_analysis_plots/'+night+'_2d_powerspectrum.pdf')

#######
print
print "finding x-ray/optical data overlap..."

## keep just overlapping data
ticks = np.arange(len(xts))
superchunk_lists= [list(zip(*group)) for k, group in groupby(zip(ticks,xts), lambda x: x[1] == 0) if not k]
n_superchunks = np.shape(superchunk_lists)[0]

t_init = [] #do i ever use t_init?
for x in superchunk_lists:
    t_init.append(tts[np.min(x[0])])
print "number of contiguous overlapping data chunks (any size): ", n_superchunks
print "total time overlapping data: ", np.sum([np.shape(x)[1] for x in superchunk_lists]), "s (", \
    np.sum([np.shape(x)[1] for x in superchunk_lists])/60., "min)"
## can look at times and combine some of them by hand if necessary????

#######
print
print "splitting data into 64-s chunks for just overlap data..."
## just dealing with overlap data
## split into 64s chunks
csize = 64
k = 0
datadict = {}
# fig5,ax = plt.subplots(facecolor='w')
# fig5.suptitle('X-ray chunk timeseries')
for i in np.arange(n_superchunks):
    print "contiguous x-ray data segment: ", i
    print "initial time: ", t_init[i]
    x = superchunk_lists[i]
    superchunk_size = np.shape(x)[1]
    tickdata = x[0]
    actualtimedata = [times[n] for n in tickdata]
    xdata = [xts[n] for n in tickdata]
    odata = [ots[n] for n in tickdata]
    tdata = [tts[n] for n in tickdata]
    chunkticks = np.arange(csize)
    #     plt.plot(actualtimedata,xdata,'.')

    xts_c = chunks(xdata, csize)
    ots_c = chunks(odata, csize)
    tts_c = chunks(tdata, csize)
    time_c = chunks(actualtimedata, csize)
    tc = np.arange(csize)
    ctot = len(ots_c)
    print "# of chunks in this segment: ", ctot
    print

    for c in xrange(ctot):
        # cm = median-subtracted cleaned timeseries
        xts_cm = np.zeros(np.shape(xts_c[c]))
        ots_cm = np.zeros(np.shape(ots_c[c]))
        tts_cm = np.zeros(np.shape(tts_c[c]))

        # m is the median of the timeseries
        xts_m = np.median(xts_c[c][xts_c[c] > 0])
        ots_m = np.median(ots_c[c][ots_c[c] > 0])
        tts_m = np.median(tts_c[c][tts_c[c] > 0])

        # fixed = set zeroes equal to median
        xts_f = xts_cm + xts_m
        ots_f = ots_cm + ots_m
        tts_f = tts_cm + tts_m

        for i in xrange(ntot):
            if xts_c[c][i] > 0:
                xts_cm[i] = xts_c[c][i] - xts_m
                xts_f[i] = xts_c[c][i]
            if ots_c[c][i] > 0:
                ots_cm[i] = ots_c[c][i] - ots_m
                ots_f[i] = ots_c[c][i]
            if tts_c[c][i] > 0:
                tts_cm[i] = tts_c[c][i] - tts_m
                tts_f[i] = tts_c[c][i]

        #         print xts_cm
        # ax.plot(tc, xts_cm)  # all 64-s chunks of data, overlapping
        chunkID = np.arange(ctot)
        chunkdict = {'unix_time': time_c[c],
                     'chunktime': tc,
                     'xray_maxij_clean_ts': xts_c[c],
                     'rho_maxij_clean_ts': ots_c[c],
                     'rho_tycho_clean_ts': tts_c[c],

                     'xray_maxij_fixed_ts': xts_f,
                     'rho_maxij_fixed_ts': ots_f,
                     'rho_tycho_fixed_ts': tts_f,

                     'xray_maxij_med_ts': xts_cm,
                     'rho_maxij_med_ts': ots_cm,
                     'rho_tycho_med_ts': tts_cm
                     }

        #     print c
        #     print ots_cm
        chunkframe = pd.DataFrame(chunkdict, index=tc)
        datadict[str(k)] = chunkframe
        #         print "k = ", k
        k += 1
print 'total # of overlapping chunks:', len(datadict)
ochunk_datadict = datadict
ochunk_keys = np.arange(k)

#

# # deletethese = [6, 10, 14] # for 3/28 THIS NEEDS TO BE EDITED FOR EACH --> moved to top
ochunk_kfix = np.setdiff1d(ochunk_keys,deletethese)
# print ochunk_kfix



# plot of time series over chunk of data
ctot = len(ochunk_kfix)

fig6,ax1=plt.subplots(figsize=[6,10],facecolor='w')
subplots_adjust(hspace=0.000)
for c,v in zip(ochunk_kfix,xrange(ctot)):
    ch = str(c)
    v = v+1
    ax1 = subplot(ctot,1,v)
    ax1.plot(tc, ochunk_datadict[ch]['xray_maxij_med_ts'],'.-b')
    ax1.plot(tc, ochunk_datadict[ch]['rho_maxij_med_ts'],'.-r')
    # ax1.plot(tc, ochunk_datadict[ch]['rho_tycho_med_ts'],'.-k',alpha=.3)
    ax1.plot([0,ntot],[0,0],'k')
    ax1.set_xlim(0,ntot)
ax1.legend(['MaxiJ@Nicer','MaxiJ@RHO',"Tycho@RHO"])
plt.savefig('overlap_analysis_plots/'+night+'_time_series_median_chunks.pdf')

#######
print
print "calculating ACFs and CCFs of overlapping data..."

## ACFs & CCFs:
acf_xray_maxij = {}
acf_rho_maxij = {}
acf_rho_tycho = {}

ccf_xray_rho_maxij = {}
ccf_xray_maxij_rho_tycho_control = {}
crosscorr_rho_maxij_tycho_weather = {}

for k in ochunk_kfix:
    ch = str(k)
    xts_cm = ochunk_datadict[ch]['xray_maxij_med_ts']
    ots_cm = ochunk_datadict[ch]['rho_maxij_med_ts']
    tts_cm = ochunk_datadict[ch]['rho_tycho_med_ts']

    acf_xray_maxij[ch] = fftconvolve(xts_cm, xts_cm[::-1], mode='same')  # x-ray autocorrelation
    acf_rho_maxij[ch] = fftconvolve(ots_cm, ots_cm[::-1], mode='same')  # optical autocorrelation
    acf_rho_tycho[ch] = fftconvolve(tts_cm, tts_cm[::-1], mode='same')  # tycho reference star autocorrelation

    ccf_xray_rho_maxij[ch] = fftconvolve(xts_cm, ots_cm[::-1], mode='same')  # cross correlation
    ccf_xray_maxij_rho_tycho_control[ch] = fftconvolve(xts_cm, tts_cm[::-1], mode='same')  # cross corr'ln w/ refstar
    crosscorr_rho_maxij_tycho_weather[ch] = fftconvolve(ots_cm, tts_cm[::-1], mode='same')


xray_maxij_ts = {}
rho_maxij_ts = {}
rho_tycho_ts = {}

overlapdict_fixed = {}
overlapdict_med = {}

for c in ochunk_kfix:
    ch = str(c)

    # overlapdict_fixed = dictionary for acf/ccf/fft with cleaned data (zeros set to median)
    xray_maxij_ts[ch] = ochunk_datadict[ch]['xray_maxij_fixed_ts']
    rho_maxij_ts[ch] = ochunk_datadict[ch]['rho_maxij_fixed_ts']
    rho_tycho_ts[ch] = ochunk_datadict[ch]['rho_tycho_fixed_ts']

    chunkdict = {'acf_xray_maxij': fftconvolve(xray_maxij_ts[ch], xray_maxij_ts[ch][::-1]),
                 'acf_rho_maxij': fftconvolve(rho_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 'acf_rho_tycho': fftconvolve(rho_tycho_ts[ch], rho_tycho_ts[ch][::-1]),

                 # ccf_xray_rho_maxij[ch] = fftconvolve(xts_cm, ots_cm[::-1], mode='same')  # cross correlation
                 'ccf_xray_rho_maxij': fftconvolve(xray_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 # cross-correlate = CCF of xray w/ optical (maxij)
                 'ccf_xray_maxij_rho_tycho_control': fftconvolve(xray_maxij_ts[ch], rho_tycho_ts[ch][::-1]),
                 # ccf of xray w/ tycho

                 'fft_xray_maxij': np.fft.fft(xray_maxij_ts[ch]),
                 'fft_rho_maxij': np.fft.fft(rho_maxij_ts[ch]),
                 'fft_rho_tycho': np.fft.fft(rho_tycho_ts[ch])
                 }
    overlapdict_fixed[ch] = chunkdict

    #overlapdict_med = dictionary for acf/ccf/fft with median-subtracted data
    xray_maxij_ts[ch] = ochunk_datadict[ch]['xray_maxij_med_ts']
    rho_maxij_ts[ch] = ochunk_datadict[ch]['rho_maxij_med_ts']
    rho_tycho_ts[ch] = ochunk_datadict[ch]['rho_tycho_med_ts']

    chunkdict = {'acf_xray_maxij': fftconvolve(xray_maxij_ts[ch], xray_maxij_ts[ch][::-1]),
                 'acf_rho_maxij': fftconvolve(rho_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 'acf_rho_tycho': fftconvolve(rho_tycho_ts[ch], rho_tycho_ts[ch][::-1]),

                 # ccf_xray_rho_maxij[ch] = fftconvolve(xts_cm, ots_cm[::-1], mode='same')  # cross correlation
                 'ccf_xray_rho_maxij': fftconvolve(xray_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 # cross-correlate = CCF of xray w/ optical (maxij)
                 'ccf_xray_maxij_rho_tycho_control': fftconvolve(xray_maxij_ts[ch], rho_tycho_ts[ch][::-1]),
                 # ccf of xray w/ tycho

                 'fft_xray_maxij': np.fft.fft(xray_maxij_ts[ch]),
                 'fft_rho_maxij': np.fft.fft(rho_maxij_ts[ch]),
                 'fft_rho_tycho': np.fft.fft(rho_tycho_ts[ch])
                 }
    #     chunkframe = pd.DataFrame(chunkdict)
    overlapdict_med[ch] = chunkdict



#######
print
print "plotting CCFs..."

fig7,ax=plt.subplots(figsize=[20,12],facecolor='w')
subplots_adjust(hspace=0.0000)
subplots_adjust(wspace=0.075000)
ctot = len(ochunk_kfix)
# fig7.tight_layout()
fig7.subplots_adjust(top=.95)
fig7.suptitle('CCFs of all overlapping chunks of 64-second data',fontsize=16)

for c,v in zip(ochunk_kfix,xrange(ctot)):
    v = v+1
    ax = subplot(np.ceil(ctot/2)+1,2,v)
    ch = str(c)
    ax.plot(overlapdict_med[ch]['ccf_xray_rho_maxij'], 'g', linewidth=2)
    ax.plot(overlapdict_med[ch]['ccf_xray_maxij_rho_tycho_control'], 'c', linewidth=2)
    ax.axvline(x=ntot/2.,color='y',alpha=.7)
    ax.axhline(y=0,color='k')
    titletime = ochunk_datadict[ch]['unix_time'][0]
    ax.set_title(datetime.utcfromtimestamp(titletime).strftime('%Y-%m-%d %H:%M:%S') + " UTC", x=.83, y=0)
    ax.set_xlim(0,ntot)
#     ax.set_xlim(20,40)

ax.legend(['MaxiJ@Nicer/MaxiJ@RHO',
           'MaxiJ@Nicer/Tycho@RHO'],
         loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=3)

plt.savefig('overlap_analysis_plots/'+night+'_overlap_CCFs.pdf')

print "plotting mean of all CCFs..."
# mean of all (good) CCFs

fig8,ax=plt.subplots(figsize=[10,8],facecolor='w')

ax.plot(xrange(-32,32),np.mean([ccf_xray_rho_maxij[str(c)] for c in ochunk_kfix],axis=0),'g', lw=2)
ax.plot(xrange(-32,32),np.mean([ccf_xray_maxij_rho_tycho_control[str(c)] for c in ochunk_kfix],axis=0),'k',alpha=.8,lw = 2)
ax.axvline(x=0,color='r',lw = 2)
ax.axhline(y=0,color='k')
ax.legend(['NICER BH : RHO BH','NICER BH : Reference Star (Noise)'],
         loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=1, fontsize=20)

ax.set_title("mean ccf of x-ray to optical in a bunch of one-minute data segments")
ax.set_xlim(-32,32)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('Seconds',fontsize=16)
ax.set_ylabel('Amplitude',fontsize=16)
plt.savefig('overlap_analysis_plots/'+night+'_mean_all_CCFs.pdf')

xray_maxij_ts = {}
rho_maxij_ts = {}
rho_tycho_ts = {}

datadict2 = {}

for c in ochunk_kfix:
    ch = str(c)
    xray_maxij_ts[ch] = ochunk_datadict[str(c)]['xray_maxij_fixed_ts']
    rho_maxij_ts[ch] = ochunk_datadict[str(c)]['rho_maxij_fixed_ts']
    rho_tycho_ts[ch] = ochunk_datadict[str(c)]['rho_tycho_fixed_ts']

    chunkdict = {'acf_xray_maxij': fftconvolve(xray_maxij_ts[ch], xray_maxij_ts[ch][::-1]),
                 'acf_rho_maxij': fftconvolve(rho_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 'acf_rho_tycho': fftconvolve(rho_tycho_ts[ch], rho_tycho_ts[ch][::-1]),

                 'ccf_xray_rho_maxij': fftconvolve(xray_maxij_ts[ch], rho_maxij_ts[ch][::-1]),
                 # cross-correlate = CCF of xray w/ optical (maxij)
                 'ccf_xray_maxij_rho_tycho_control': fftconvolve(xray_maxij_ts[ch], rho_tycho_ts[ch][::-1]),
                 # ccf of xray w/ tycho

                 'fft_xray_maxij': np.fft.fft(xray_maxij_ts[ch]),
                 'fft_rho_maxij': np.fft.fft(rho_maxij_ts[ch]),
                 'fft_rho_tycho': np.fft.fft(rho_tycho_ts[ch])
                 }
    #     chunkframe = pd.DataFrame(chunkdict)
    datadict2[ch] = chunkdict



ch = str(ochunk_kfix[5])
nel2=len(datadict2[ch]['fft_xray_maxij'])/2
f = np.arange(nel2)/64.

fig10,ax=plt.subplots(figsize=[6,15],facecolor='w')
subplots_adjust(hspace=0.0000)
ctot = len(ochunk_kfix)
# fig10.tight_layout()
# fig10.subplots_adjust(top=.5)
fig10.suptitle('Power Spectra - Maxi J1820+070 \n (64s each, ' + night + ')', fontsize=16)
for c,v in zip(ochunk_kfix,xrange(ctot)):
    v = v+1
    ax = subplot(ctot+1,1,v)
    ch = str(c)
    ax.plot(f,np.abs(((datadict2[ch]['fft_xray_maxij'][:nel2]))**2)/np.abs(((datadict2[ch]['fft_xray_maxij'][:nel2]))**2)[0],'b-', lw=2, drawstyle='steps-mid') #xray maxij
    ax.plot(f,10*np.abs(((datadict2[ch]['fft_rho_maxij'][:nel2]))**2)/np.abs(((datadict2[ch]['fft_rho_maxij'][:nel2]))**2)[0],'r-', lw=2, drawstyle='steps-mid') #optical maxij
    # ax.plot(f,10*np.abs(((datadict2[ch]['fft_rho_tycho'][:nel2]))**2)/np.abs(((datadict2[ch]['fft_rho_tycho'][:nel2]))**2)[0],'g.-', drawstyle='steps-mid') #tycho
    # ax.axvline(x=ntot/2.,color='y',alpha=.7)
    # ax.axhline(y=0,color='k')
    ax.axvline(x=0.09, c='k', linestyle='--', lw=2)
    ax.axvline(x=0.045, c='k', linestyle='--', lw=2)

    ax.get_yaxis().set_ticks([])
    if c != np.max(ochunk_kfix):
        ax.get_xaxis().set_ticks([])

    titletime = ochunk_datadict[ch]['unix_time'][0]
    ax.set_title(datetime.utcfromtimestamp(titletime).strftime('%Y-%m-%d %H:%M:%S') + " UTC", x=.6, y=.5)

# ax.plot(crosscorr_xo - crosscorr_xt,'k')
# ax.plot(auto_o/4,'r')
# ax.plot(auto_x,'b')
# ax.plot(auto_t/3,'k')
#     ax.plot([np.argmax(auto_x),np.argmax(auto_x)],[-4.e9,2.e10],'y',alpha=.7)
    ax.set_ylim(0,.01)
#     ax.set_xlim(0,ntot)
    ax.set_xlim(0,.5)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Normalized Power', y=6)
# ax.set_title('Optical+X-ray Correlation (one chunk of 128-second overlapping data) 51.5s lead')
ax.legend(['MaxiJ@Nicer',
           'MaxiJ@Nicer',
           'Previously Identified QPOs'
           # 'Tycho@RHO',
           # 'RHO Optical Auto-correlation',
           # 'Nicer X-ray Auto-correlation',
           # 'RHO Tycho Reference Auto-correlation',
#            'Autocorrelation Peak'
          ],
         loc='upper center', bbox_to_anchor=(.5, -0.6),
          fancybox=True, shadow=True, ncol=1)
plt.savefig('overlap_analysis_plots/' + night + '_fft_overlaps.pdf')
# print np.argmax(crosscorr_xo)


fig11, ax = plt.subplots(figsize=[6, 15], facecolor='w')
subplots_adjust(hspace=0.0000)
ctot = len(ochunk_kfix)
# fig11.tight_layout()
# fig11.subplots_adjust(top=.95)
fig11.suptitle('Crosspower Spectra - Maxi J1820+070 \n (64s each, ' + night + ')', fontsize=16)

for c, v in zip(ochunk_kfix, xrange(ctot)):
    v = v + 1
    ax = subplot(ctot+1, 1, v)
    ch = str(c)

    xps = np.abs((datadict2[ch]['fft_xray_maxij'][:nel2] * np.conj(datadict2[ch]['fft_rho_maxij'][:nel2]))) / \
          np.abs((datadict2[ch]['fft_xray_maxij'][:nel2] * np.conj(datadict2[ch]['fft_rho_maxij'][:nel2])))[0]
    tycho_xps = np.abs((datadict2[ch]['fft_xray_maxij'][:nel2] * np.conj(datadict2[ch]['fft_rho_tycho'][:nel2]))) / \
                np.abs((datadict2[ch]['fft_xray_maxij'][:nel2] * np.conj(datadict2[ch]['fft_rho_tycho'][:nel2])))[0]

    ax.plot(f[f > 0], xps[f > 0], 'g-', drawstyle='steps-mid', lw=3)
    ax.plot(f[f > 0], tycho_xps[f > 0], '-k', drawstyle='steps-mid', alpha=.5, lw=3)
    # ax.axvline(x=0.0225,color='y',linestyle='--')
    ax.axvline(x=0.09, c='k', linestyle='--', lw=2)
    ax.axvline(x=0.045, c='k', linestyle='--', lw=2)
    # ax.axvline(x=0.18,color='y',linestyle='--')

    ax.set_ylim(0, .002)
    ax.set_xlim(0, .5)
    ax.get_yaxis().set_ticks([])
    if c != np.max(ochunk_kfix):
        ax.get_xaxis().set_ticks([])

    titletime = ochunk_datadict[ch]['unix_time'][0]
    ax.set_title(datetime.utcfromtimestamp(titletime).strftime('%Y-%m-%d %H:%M:%S') + " UTC", x=.6, y=.5)


ax.legend(['Maxi J1820+070 X-ray:Optical', 'Maxi J1820+070 X-ray:Reference Star Optical', 'Previously Identified QPOs'],
          loc='upper center', bbox_to_anchor=(.5, -0.6),
          fancybox=True, shadow=True, ncol=1)


ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Normalized Power', y=6)
plt.savefig('overlap_analysis_plots/' + night + '_crosspower_spectra.pdf')

#######
print
print timefinish(time0)
# plt.show()