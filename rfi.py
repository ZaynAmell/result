import blimpy as bl
from scipy import signal
import numpy as np
from scipy import interpolate

def localmax(wf,pro=(20000,3.0e6),dist=30):
    """This funcation takes into a waterfall file
        and returns its local maxima.
        
        wf:the waterfall file
        pro: the difference between the peak to its base line
        dist:the distance between a freqency that are considered to be within a range
        """
    #acquire the power and frequency
    plot_f,plot_data=wf.get_power()
    peak=signal.find_peaks(plot_data[1:],prominence=pro) #find local maxima
    #calculate the frequency range
    i=0
    rfi=[]
    while i<len(peak[0])-1:
        start=plot_f[peak[0][i]]
        if(plot_f[peak[0][i+1]]-start<=dist):#if the frequency is alone
            end=plot_f[peak[0][i+1]]
            length=len(peak[0][i+1:])
            for j in range(length): #determine the range of RFI
                if j==length-1:
                    end=plot_f[peak[0][i+1:][j]]
                    i+=j+1
                    rfi.append([start,end])
                    break
                else:
                    if (abs(plot_f[peak[0][i+1:][j]]-end)<=dist):
                        end=plot_f[peak[0][i+1:][j]]
                    else:
                        i+=j+1
                        if(j==0):
                            rfi.append([start,start])
                        else:
                            rfi.append([start,end])
                        break
        else:
            rfi.append([start,start])
            i+=1
    return rfi

def baseline(wf,spl_order=None):
    channel_width=wf.container.n_channels_in_file
    if spl_order==None:
        spl_order=8
    x=wf.get_freqs()
    integrated_channel=wf.get_power()[1]
    knots = np.arange(0, channel_width, channel_width//spl_order+1)
    spl = interpolate.splrep(x, integrated_channel,t=wf.container.f_start+knots[1:]*wf.header['foff'])
    chan_fit = interpolate.splev(x, spl)
    return chan_fit
