import numpy as np

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


def rd_nicer_lc(pathnam, fnams):
    
    nfils=len(fnams)
    xlist=list()
    tlist=list()
    
    for i in range(nfils):
        fnam=pathnam+fnams[i]
    
        fID=open(fnam,'r')
    
    
        for j in range(16):
            line=fID.readline() # skip first 10 lines
        
        while True:
            line=fID.readline()
            if not line: break
            jnk=line.split()
        
            #print len(jnk)
            if len(jnk)==14:
                xlist.append(float(jnk[2]))
                tlist.append(float(jnk[1]))
            
        fID.close()
    
    return np.array(tlist).astype(float),np.array(xlist).astype(float)

def rd_nicer_lc34(pathnam, fnams):
    # bands 3 and 4
    
    nfils=len(fnams)
    xlist=list()
    tlist=list()
    
    for i in range(nfils):
        fnam=pathnam+fnams[i]
    
        fID=open(fnam,'r')
    
    
        for j in range(16):
            line=fID.readline() # skip first 10 lines
        
        while True:
            line=fID.readline()
            if not line: break
            jnk=line.split()
        
            #print len(jnk)
            if len(jnk)==14:
                xlist.append(float(jnk[6])+float(jnk[7]))
                tlist.append(float(jnk[1]))
            
        fID.close()
        
    return np.array(tlist).astype(float),np.array(xlist).astype(float)

def rd_nicer_lc0(pathnam, fnams):
    # band 0
    
    nfils=len(fnams)
    xlist=list()
    tlist=list()
    
    for i in range(nfils):
        fnam=pathnam+fnams[i]
    
        fID=open(fnam,'r')
    
    
        for j in range(16):
            line=fID.readline() # skip first 10 lines
        
        while True:
            line=fID.readline()
            if not line: break
            jnk=line.split()
        
            #print len(jnk)
            if len(jnk)==14:
                xlist.append(float(jnk[3]))
                tlist.append(float(jnk[1]))
            
        fID.close()
    return np.array(tlist).astype(float),np.array(xlist).astype(float)

def rd_nicer_lc1(pathnam, fnams):
    # band 1
    
    nfils=len(fnams)
    xlist=list()
    tlist=list()
    
    for i in range(nfils):
        fnam=pathnam+fnams[i]
    
        fID=open(fnam,'r')
    
    
        for j in range(16):
            line=fID.readline() # skip first 10 lines
        
        while True:
            line=fID.readline()
            if not line: break
            jnk=line.split()
        
            #print len(jnk)
            if len(jnk)==14:
                xlist.append(float(jnk[4]))
                tlist.append(float(jnk[1]))
            
        fID.close()

    return np.array(tlist).astype(float),np.array(xlist).astype(float)

def rd_nicer_lc2(pathnam, fnams):
    # band 2
    
    nfils=len(fnams)
    xlist=list()
    tlist=list()
    
    for i in range(nfils):
        fnam=pathnam+fnams[i]
    
        fID=open(fnam,'r')
    
    
        for j in range(16):
            line=fID.readline() # skip first 10 lines
        
        while True:
            line=fID.readline()
            if not line: break
            jnk=line.split()
        
            #print len(jnk)
            if len(jnk)==14:
                xlist.append(float(jnk[5]))
                tlist.append(float(jnk[1]))
            
        fID.close()
    
    return np.array(tlist).astype(float),np.array(xlist).astype(float)

def nicer_lc_2_1sts(t,cts):
    # makes 1s time series from NICER lightcurve
    
    t0=np.min(t)
    tf=np.max(t)
    t_span=np.round(tf-t0).astype(int)
    t2=t-t0
    
    yy=np.zeros(t_span)
    
    for i in range(t_span):
        if i%1000 == 0:
            print i, t_span
        j=[k for k in range(len(t)) if t2[k]>=i and t2[k]<i+1]
        if len(j)>0:
            # yy[i]=yy[i]+np.sum(cts[j])
            yy[i] = yy[i] + np.mean(cts[j])
    return yy

def nicer_lc_2_ts(t,cts):
    # makes 1/8s time series from NICER lightcurve
    
    t0=np.min(t)
    tf=np.max(t)
    t_span=8*np.round(tf-t0).astype(int)
    t2=t-t0
    
    yy=np.zeros(t_span+20)
    
    tindex=(8*t2).astype(int)
    yy[tindex]=cts[:]
    
    return yy



def photlc_2_ts(t5,maxiphot5):
    # takes t, phot and turns it into a time series
    
    dt=(np.max(t5)-np.min(t5)).astype(int)
    print dt
    f=np.zeros(dt+1)
    
    t0=(t5-np.min(t5)).astype(int)
    f[t0]=maxiphot5[:]
    
    return f

def pds_chunk128_wgaps(t,flx,t1,t2):
    # takes time series - not lightcurve - time and flx
    # extracts time from t1 to t2
    # calculates power spectrum in 128s segments and returns
    # skips gaps
    
    t_span=t2-t1
    n_chunk=np.round(t_span/128.0).astype(int)
    
    print t_span, n_chunk
    
    dt2=128.0/2.0 # half-width of chunk
    
    for i in range(n_chunk):
        ty=t-((i*dt2*2)+dt2+t1)
        y=flx[np.abs(ty)<=dt2]
        if np.sum(y) > 3e5:
            af=np.fft.fft(y)
            if i==0:
                flen=len(y)/2
                print i, flen
                p=np.zeros([n_chunk,flen])
                f=np.arange(flen)/(dt2*2)
            p[i,:]=np.abs(af[:flen])**2
        
    ptot=np.sum(p,axis=0)
    
    return f,p,ptot

def pds_chunk128(t,flx,t1,t2):
    # takes time series - not lightcurve - time and flx
    # extracts time from t1 to t2
    # calculates power spectrum in 128s segments and returns
    
    t_span=t2-t1
    n_chunk=np.round(t_span/128.0).astype(int)
    
    print t_span, n_chunk
    
    dt2=128.0/2.0 # half-width of chunk
    
    for i in range(n_chunk):
        ty=t-((i*dt2*2)+dt2+t1)
        y=flx[np.abs(ty)<=dt2]
        if i==0:
            flen=len(y)/2
            print i, flen
            p=np.zeros([n_chunk,flen])
            f=np.arange(flen)/(dt2*2)
        af=np.fft.fft(y)
        p[i,:]=np.abs(af[:flen])**2
        
    ptot=np.sum(p,axis=0)
    
    return f,p,ptot


def pds_single_chunk(flx):
    clength = len(flx)
    t = xrange(clength)
    freqlength = clength / 2
    freq = np.arange(freqlength) / (np.float(clength))
    af = np.fft.fft(flx)
    power = np.abs(af[:freqlength]) ** 2

    return freq, power

def make_1s_ts(tout,tin,cts):
    # makes 1s time series from any lightcurve
    # tout should be an array of integer seconds that the output ts will be matched to
    # tin and cts are the input times and fluxes in the lightcurve

    t0 = np.min(tout)
    tf = np.max(tout)
    t_span = np.round(tf - t0).astype(int)+1

    yy = np.zeros(t_span)

    for i in range(t_span):
        # print "second " + str(tout[i])
        j = (tin >= tout[i]) & (tin < tout[i]+1)

        if any(j):
            # print cts[j]
            yy[i] = yy[i] + np.median(cts[j])
    return yy


def chunks(l, n):
    outl = []
    j = np.ceil(len(l) / float(n))
    for i in range(int(j)):
        l0 = l[i * n:(i + 1) * n]
        l1 = np.pad(l0, [0, n - len(l0)], 'constant', constant_values=0)
        outl.append(l1)

    return outl


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

# f2,p2max=p_2_logp(f,pmax,0.03,8)