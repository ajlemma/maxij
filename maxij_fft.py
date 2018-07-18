#######
''' maxij_fft.py
by A Townsend
'''
#######

import pandas as pd
import numpy as np
from maxijdefs import *

# import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

#######

def pds_chunk128(t, flx, t1, t2):
    # takes time series - not lightcurve - time and flx
    # extracts time from t1 to t2
    # calculates power spectrum in 128s segments and returns

    t_span = t2 - t1
    n_chunk = np.round(t_span / 128.0).astype(int)

    print t_span, n_chunk

    dt2 = 128.0 / 2.0  # half-width of chunk

    for i in range(n_chunk):
        ty = t - ((i * dt2 * 2) + dt2 + t1)
        y = flx[np.abs(ty) <= dt2]
        if i == 0:
            flen = len(y) / 2
            print i, flen
            p = np.zeros([n_chunk, flen])
            f = np.arange(flen) / (dt2 * 2)
        af = np.fft.fft(y)
        p[i, :] = np.abs(af[:flen]) ** 2

    ptot = np.sum(p, axis=0)

    return f, p, ptot

def p_2_logp(f, p, df, n_avg):
    # takes output of pds_chunk128-style program
    # averages over n_avg pieces, and logbins the 2D array

    nt = p.shape[0]  # number of power spectra

    for i in range(nt):
        i1 = i - n_avg / 2
        i2 = i1 + n_avg / 2
        if i1 < 0:
            i1 = 0
        if i2 > (nt - 1):
            i2 = nt - 1
        pmean = np.mean(p[i1:i2 + 1, :], axis=0) / (np.mean(p[i1:i2 + 1, 0]) ** 0.5)
        dumf, dump = logbin(f, pmean, df)
        if i == 0:
            nf = len(dump[dump > 0])
            p2 = np.zeros([nt, nf])
            f2 = dumf[dump > 0]
        p2[i, :] = dump[dump > 0]

    return f2, p2

def logbin(f,p,df):
    #df = step size in log frequency
    fx=np.zeros(len(f))
    px=np.zeros(len(p))
    fx[:]=f[:]
    fx[0]=fx[1]
    p0=p[0]
    px[:]=p[:]
    px[0]=0.

    fl=np.log10(fx)
    nel=len(f)

    fmin=fl[1]
    fmax=np.max(fl)

    nf=np.round((fmax-fmin)/df).astype(int)
    f2=np.arange(nf)*df+fmin
    p2=np.zeros(nf)

    for i in range(nf):
        fa=10.**f2[i]
        fb=10.**(f2[i]+df)
        fm=(fa+fb)/2.0
        dff=fm-fa
        if len(f[np.abs(f-fm)<=dff])>0:
            p2[i]=np.mean(px[np.abs(f-fm)<=dff])

    return f2,p2

## set night
path = './maxij_data_copies/'
night = '2018-05-11'

## read in data as pandas dataframe:
data = pd.read_pickle(path+'data_'+night+'.pkl')
fID = data['fileID']
# print data

## make a new column for the data to poke at
data['maxij_flux'] = data['phot_maxij']
data['tyc_flux'] = data['phot_tyc']
data['ref_flux'] = data['phot_ref3']+data['phot_ref4']+data["phot_ref5"]

## check for NaNs and set them to 0 if they exist, ie:
i = data[np.isnan(data['maxij_flux'])]['fileID']
print 'maxij ', i
data.loc[i,'maxij_flux'] = 0.0

i = data[np.isnan(data['tyc_flux'])]['fileID']
print 'tyc ', i
data.loc[i,'tyc_flux'] = 0.0

i = data[np.isnan(data['ref_flux'])]['fileID']
print 'ref ', i
data.loc[i,'ref_flux'] = 0.0

## get FFTs for 128-s chunks of data
data_0 = 0
data_f =  len(data)
ts = np.arange(data_f-data_0)
tseries = fID[data_0:data_0+data_f]

mphot = data.loc[tseries,'maxij_flux']
f,pmax,pmaxtot = pds_chunk128(ts,mphot,data_0,data_f)

tphot = data.loc[tseries,'tyc_flux']
f,ptyc,ptyctot = pds_chunk128(ts,tphot,data_0,data_f)

rphot = data.loc[tseries,'ref_flux']
f,pref,preftot = pds_chunk128(ts,rphot,data_0,data_f)


f2,p2max=p_2_logp(f,pmax,0.03,4)
f2,p2tyc=p_2_logp(f,ptyc,0.03,4)
f2,p2ref=p_2_logp(f,pref,0.03,4)

# Set up figure and image grid
fig = plt.figure(figsize=(9,12), facecolor='w')
fig.suptitle('RHO 1Hz Photometry - %s' % night, fontsize=20)

grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,3),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="10%",
                 cbar_pad=0.15,
                 )

# Add data to image grid
ax1 = grid[0]
ax1.set_title('MAXIJ1820+070')
im1 = ax1.imshow((10.**f2[:1])*(p2max[:,1:]),vmax=40,cmap='viridis',interpolation='none')

ax2 = grid[1]
ax2.set_title('TYC 444-2244-1')
im2 = ax2.imshow((10.**f2[:1])*(p2tyc[:,1:]),vmax=40,cmap='viridis',interpolation='none')

ax3=grid[2]
ax3.set_title('Sum of 3 other Ref Stars')
im3 = ax3.imshow((10.**f2[:1])*(p2ref[:,1:]),vmax=40,cmap='viridis',interpolation='none')


# Colorbar
ax3.cax.colorbar(im1)
ax3.cax.toggle_label(True)

# plt.tight_layout()    # Works, but may still require rect paramater to keep colorbar labels visible
plt.show()
fig.savefig('2d_qpo_plots/%s_2d.pdf' % night)