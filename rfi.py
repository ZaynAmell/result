import blimpy as bl
from blimpy import Waterfall
from scipy import signal
import numpy as np
from scipy import interpolate
import pandas as pd
import peakutils
import BaselineRemoval

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
        rm=BaselineRemoval.BaselineRemoval(p).ZhangFit()
    rm=np.where(rm<0,0,rm)
    base=peakutils.baseline(rm,deg=deg)
    base=np.where(base<0,0,base)
    return rm, base

def base(w,deg=None):
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
        r,base=peakbase(f,p)
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

def intersection(w,multi=None):
    if type(w) != Waterfall: #Check the file type
        print("wrong file type")
        return
    f,p,b=base(w)
    #set threshold
    maxi=np.amax(p-b)
    if multi==None:
        multi=1.0
    if maxi>10000:
        multi=multi/10
    thres=multi*(p-b).std()
    
    rm=p-b
    #Determin the range above threshold
#     intercept=[]
#     if rm[0]>=thres:
#         intercept.append(f[0])
#     for i in range(1,len(rm)-1):
#         if rm[i]==thres:
#             intercept.append(f[i])
#         elif (rm[i]-thres)*(rm[i+1]-thres)<=0:
#             x_i=(thres-rm[i])*(f[i+1]-f[i])/(rm[i+1]-rm[i])+f[i]
#             intercept.append(x_i)
    index=np.where(rm>=thres)
    return f[index]

def rfi_range(w,multi=None):
    intercept=intersection(w,multi=multi)
    #Organize the RFI frequency range
    i=0
    test=[]
    dist=w.header['foff']
    while i<len(intercept)-1:
        start=intercept[i]
        if(intercept[i+1]-start<=dist):
            end=intercept[i+1]
            length=len(intercept[i+1:])
            for j in range(length):
                if j==length-1:
                    end=intercept[i+1:][j]
                    i+=j+1
                    test.append([start,end])
                    break
                else:
                    if (abs(intercept[i+1:][j]-end)<=dist):
                        end=intercept[i+1:][j]
                    else:
                        i+=j+1
                        if(j==0):
                            test.append([start,start])
                        else:
                            test.append([start,end])
                        break
        else:
            test.append([start,start])
            i+=1
    return test