import blimpy as bl
from blimpy import Waterfall
from scipy import signal
import numpy as np
from scipy import interpolate
import pandas as pd
import peakutils
import BaselineRemoval
#This files require to add a get_power() function to blimpy/plotting/
def splbase(f,p,spl_order=None):
    if spl_order==None:
        spl_order=16
    knots = np.arange(f[0],f[-1],(f[-1]-f[0])//spl_order+1)
    spl = interpolate.splrep(f, p,t=knots[1:])
    chan_fit = interpolate.splev(f, spl)
    return chan_fit

def peakbase(f,p,type,deg=3):
    if type==1:
        rm=BaselineRemoval.BaselineRemoval(p).IModPoly()
    if type==0:
        rm=BaselineRemoval.BaselineRemoval(p).ZhangFit(itermax=100)
    rm=np.where(rm<0,0,rm)
    base=peakutils.baseline(rm,deg=deg)
    base=np.where(base<0,0,base)
    return rm, base

def base(w,deg=None,flip=False):
    """Acquire the possible RFI with a waterfall file
        w:waterfall file
        deg: the degree of polynoimal to do a baseline fit
        thres: threshold to be considered as RFI (power level)
        dist: to organized the intercept into the same range
        """
    if type(w) != Waterfall: #Check the file type
        print("wrong file type")
        return
    if w.container.f_start==856.0 or w.container.f_start==544.0: #remove the first index for a better fitting
        f=w.get_power()[0][1:]
        p=w.get_power()[1][1:]
    else:
        f=w.get_power()[0]
        p=w.get_power()[1]
    if flip==True:
        p=np.flip(p)
    #Different situation:
    snr=p.mean()/p.std()
    #Eliminating Edge of band
    if w.container.f_start==856.0:
        index=np.where(f>=910)
        f=f[index]
        p=p[index]
    if w.container.f_stop==1712.0:
        index=np.where(f<=1650)
        f=f[index]
        p=p[index]
    if w.container.f_start==544.0:
        index=np.where(f>=650)
        f=f[index]
        p=p[index]
    if w.container.f_stop==1088.0 or w.container.f_stop==1020.0:
        index=np.where(f<=980)
        f=f[index]
        p=p[index]
    #Choose fiiting based on SNR
    if snr>=15:
        base=splbase(f,p)
    elif snr>=6:
        r,base=peakbase(f,p,0)
        p=r
    elif snr>=3:
        base=splbase(f,p)
    elif snr>=1.6:
        r,base=peakbase(f,p,0,deg=6)
        p=r
    else:
        #base=splbase(f,p)
        r,base=peakbase(f,p,1,deg=10)
        p=r
    return f, p, base

def intersection(w,multi=None,flip=False):
    if type(w) != Waterfall: #Check the file type
        print("wrong file type")
        return
    f,p,b=base(w,flip=flip)
    #set threshold
    maxi=np.amax(p-b)
    if multi==None:
        multi=1.0
    if maxi>10000:
        multi=multi/10
    thres=multi*(p-b).std()
    rm=p-b
    index=np.where(rm>=thres)
    return f[index]

def rfi_range(w,multi=None):
    fi=intersection(w,multi=multi)
    test=[]
    flag=fi[-1]
    while True:
        start=fi[0]
        if len(fi)==1:
            end=fi[0]
            test.append([start,end])
            break
        for i in range(len(fi)-1):
            if abs(fi[i+1]-fi[i])>abs(w.header['foff']):
                end=fi[i]#+w.header['foff']
                fi=fi[i+1:]
                test.append([start,end])
                break
            else:
                end=fi[i+1]
        if end==flag:
            test.append([start,end])
    return test
