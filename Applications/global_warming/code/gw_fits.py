import numpy as np
import matplotlib.pyplot as plt
import cPickle
import random
import os
import sys
import pandas as pd
import time
#from pydl.pydlutils.spheregroup import spherematch
import kali.carma 
import kali
import time
import carmcmc as cm
from astropy.io import ascii
from collections import defaultdict
#Initialization
base_dir='/home/jtrichm2/long_time_series/'
plot_dir=base_dir+'plots/'
data_dir=base_dir
outdir=base_dir

save_dir=base_dir
outfile_pi=base_dir

def plot_power_spectrum_drw(samples, time, p,percentile=68.0, nsamples=None, plot_log=True, color="b", alpha=0.5, sp=None,
                            doShow=True):
        """
        Plot the posterior median and the credibility interval corresponding to percentile of the CAR(1) PSD. This
        function returns a tuple containing the lower and upper PSD credibility intervals as a function of
        frequency, the median PSD as a function of frequency, and the frequencies.
        :rtype : A tuple of numpy arrays, (lower PSD, upper PSD, median PSD, frequencies). If no subplot axes object
            is supplied (i.e., if sp = None), then the subplot axes object used will also be returned as the last
            element of the tuple.
        :param percentile: The percentile of the PSD credibility interval to plot.
        :param nsamples: The number of MCMC samples to use to estimate the credibility interval. The default is all
                         of them. Use less samples for increased speed.
        :param plot_log: A boolean. If true, then a logarithmic plot is made.
        :param color: The color of the shaded credibility region.
        :param alpha: The transparency level.
        :param sp: A matplotlib subplot axes object to use.
        :param doShow: If true, call plt.show()
        """
        sigmas = np.array(samples['sigma'])
        l=len(sigmas)
        if p==1:
            log_omegas = np.log(samples['ar_coefs']).flatten()
        else:
            ar_coefs = np.array(samples['ar_coefs'])
            ma_coefs=np.ones((l,p))
        
        if nsamples is None:
            # Use all of the MCMC samples
            nsamples = l
        else:
            try:
                nsamples <= l
            except ValueError:
                "nsamples must be less than the total number of MCMC samples."
            
            nsamples0 = l
            index = np.arange(nsamples) * (nsamples0 / nsamples)
            sigmas = sigmas[index]
            if p==1:
                log_omegas = log_omegas[index]
            else:
                ar_coefs = ar_coefs[index,:]
                ma_coefs = ma_coefs[index,:]
        nfreq = 1000
        dt_min = time[1:] - time[0:time.size - 1]
        dt_min = dt_min.min()
        dt_max = time.max() - time.min()

        # Only plot frequencies corresponding to time scales a factor of 2 shorter and longer than the minimum and
        # maximum time scales probed by the time series.
        freq_max = 0.5 / dt_min
        freq_min = 1.0 / dt_max

        frequencies = np.linspace(np.log(freq_min), np.log(freq_max), num=nfreq)
        frequencies = np.exp(frequencies)
        psd_credint = np.empty((nfreq, 3))

        lower = (100.0 - percentile) / 2.0  # lower and upper intervals for credible region
        upper = 100.0 - lower

        numer = sigmas ** 2
        if p==1:
            omegasq = np.exp(log_omegas) ** 2
            
            for i in xrange(nfreq):
                denom = omegasq + (2. * np.pi * frequencies[i]) ** 2
                psd_samples = numer / denom
            # Now compute credibility interval for power spectrum
                psd_credint[i, 0] = np.percentile(psd_samples, lower, axis=0)
                psd_credint[i, 2] = np.percentile(psd_samples, upper, axis=0)
                psd_credint[i, 1] = np.median(psd_samples, axis=0)
        else:
      
            omega = 2.0 * np.pi * 1j * frequencies
            ar_poly = np.zeros((nfreq, nsamples), dtype=complex)
            ma_poly = np.zeros_like(ar_poly)
            for k in xrange(p):
            # Here we compute:
            #   alpha(omega) = ar_coefs[0] * omega^p + ar_coefs[1] * omega^(p-1) + ... + ar_coefs[p]
            # Note that ar_coefs[0] = 1.0.
                argrid, omgrid = np.meshgrid(ar_coefs[:, k], omega)
                ar_poly += argrid * (omgrid ** (p - k))
            ar_poly += ar_coefs[:, p-1]

            psd_samples = (np.squeeze(sigmas) ** 2 )/ np.abs(ar_poly)**2
            psd_credint[:, 0] = np.percentile(psd_samples, lower, axis=1)
            psd_credint[:, 2] = np.percentile(psd_samples, upper, axis=1)
            psd_credint[:, 1] = np.median(psd_samples, axis=1)
        # Plot the power spectra
        if sp == None:
            fig = plt.figure()
            sp = fig.add_subplot(111)

        if plot_log:
            # plot the posterior median first
            sp.loglog(frequencies, psd_credint[:, 1], color=color)
        else:
            sp.plot(frequencies, psd_credint[:, 1], color=color)

        sp.fill_between(frequencies, psd_credint[:, 2], psd_credint[:, 0], facecolor=color, alpha=alpha)
        sp.set_xlim(frequencies.min(), frequencies.max())
        sp.set_xlabel('Frequency')
        sp.set_ylabel('Power Spectrum')

        if doShow:
            plt.close()

        if sp == None:
            return (psd_credint[:, 0], psd_credint[:, 2], psd_credint[:, 1], frequencies, fig)
        else:
            return (psd_credint[:, 0], psd_credint[:, 2], psd_credint[:, 1], frequencies)



def plot_power_spectrum1(samples, time, p, percentile=68.0, nsamples=None, plot_log=True, color="b", alpha=0.5, sp=None,doShow=True):
    sigmas = np.array(samples['sigma'])
    l=len(sigmas)
    ar_coefs = np.array(samples['ar_coefs'])#np.transpose(samples['ar_coefs'])
    
    ma_coefs = np.array(samples['ma_coefs'])#np.transpose(samples['ma_coefs'])
   
    
    if nsamples is None:
            # Use all of the MCMC samples
        nsamples = l#.shape[0]
    else:
        try:
            nsamples <= l #.shape[0]
        except ValueError:
                "nsamples must be less than the total number of MCMC samples."

        nsamples0 = l #.shape[0]
        index = np.arange(nsamples) * (nsamples0 / nsamples)
        sigmas = sigmas[index]
        ar_coefs = ar_coefs[index,:]
        ma_coefs = ma_coefs[index,:]

    nfreq = 1000
    dt_min = time[1:] - time[0:time.size - 1]
    dt_min = dt_min.min()
    dt_max = time.max() - time.min()

        # Only plot frequencies corresponding to time scales a factor of 2 shorter and longer than the minimum and
        # maximum time scales probed by the time series.
    freq_max = 0.5 / dt_min
    freq_min = 1.0 / dt_max

    frequencies = np.linspace(np.log(freq_min), np.log(freq_max), num=nfreq)
    frequencies = np.exp(frequencies)
    psd_credint = np.empty((nfreq, 3))

    lower = (100.0 - percentile) / 2.0  # lower and upper intervals for credible region
    upper = 100.0 - lower

        # Compute the PSDs from the MCMC samples
    omega = 2.0 * np.pi * 1j * frequencies
    ar_poly = np.zeros((nfreq, nsamples), dtype=complex)
    ma_poly = np.zeros_like(ar_poly)
    for k in xrange(p):
            # Here we compute:
            #   alpha(omega) = ar_coefs[0] * omega^p + ar_coefs[1] * omega^(p-1) + ... + ar_coefs[p]
            # Note that ar_coefs[0] = 1.0.
        argrid, omgrid = np.meshgrid(ar_coefs[:, k], omega)
        ar_poly += argrid * (omgrid ** (p - k))
    ar_poly += ar_coefs[:, p-1]
    for k in xrange(ma_coefs.shape[1]):
            # Here we compute:
            #   delta(omega) = ma_coefs[0] + ma_coefs[1] * omega + ... + ma_coefs[q] * omega^q
        magrid, omgrid = np.meshgrid(ma_coefs[:, k], omega)
        ma_poly += magrid * (omgrid ** k)

    psd_samples = np.squeeze(sigmas) ** 2 * np.abs(ma_poly) ** 2 / np.abs(ar_poly)**2

        # Now compute credibility interval for power spectrum
    psd_credint[:, 0] = np.percentile(psd_samples, lower, axis=1)
    psd_credint[:, 2] = np.percentile(psd_samples, upper, axis=1)
    psd_credint[:, 1] = np.median(psd_samples, axis=1)

        # Plot the power spectra
    if sp == None:
        fig = plt.figure()
        sp = fig.add_subplot(111)

    if plot_log:
            # plot the posterior median first
        sp.loglog(frequencies, psd_credint[:, 1], color=color)
    else:
        sp.plot(frequencies, psd_credint[:, 1], color=color)

    sp.fill_between(frequencies, psd_credint[:, 2], psd_credint[:, 0], facecolor=color, alpha=alpha)
    sp.set_xlim(frequencies.min(), frequencies.max())
    sp.set_xlabel('Frequency')
    sp.set_ylabel('Power Spectrum')

    if doShow:
        plt.savefig(plot_dir+'gw_test.png')
        plt.show()

    if sp == None:
        return (psd_credint[:, 0], psd_credint[:, 2], psd_credint[:, 1], frequencies, fig)
    else:
        return (psd_credint[:, 0], psd_credint[:, 2], psd_credint[:, 1], frequencies)


def fit(t,y,yerr):
    #Read in light curve of s82_id
    p=6
    q=5
    lc=kali.lc.externalLC('gw','r',tIn=t,yIn=y,yerrIn=yerr,maskIn=np.ones(len(yerr)))
    
    newTask = kali.carma.CARMATask(p, q, nwalkers=50, nsteps=500)
    #newTask.plotpsd(doShow = False)
    #plt.savefig(plot_dir+'gw_kali_bicc_current_psd.png')
    #plt.close()
    #newTask.plotsf(LC = lc, doShow = False)
    #plt.savefig(plot_dir+'gw_kali_bicc_current_sf.png')
    #plt.close()
    
    newTask.fit(lc)
    chain=np.zeros((p+q+1,50*newTask.nsteps))
    for i in xrange(p+q+1):
        arr=[]
        for j in xrange(50):
            for k in xrange(newTask.nsteps):
                arr.append(newTask.Chain[i,j,k])
        chain[i]=arr
	d = defaultdict(list)
    for i in xrange(3):              
        if i==0:
            arr=[] 
            for j in xrange(len(chain[0])):
                arr.append(chain[:p,j])  
            d['ar_coefs']=arr          
        elif i==1:
            arr=[] 
            for j in xrange(len(chain[0])):
                arr.append(chain[p:-1,j])
            d['ma_coefs']=arr
        else:
            arr=[] 
            for j in xrange(len(chain[0])):
                arr.append(chain[-1,j]) 
            d['sigma']=arr
    if q!=0:
        psd_low, psd_hi, psd_mid, frequencies = plot_power_spectrum1(d,t,p,percentile=95.0, nsamples=2000)
    else:
        psd_low, psd_hi, psd_mid, frequencies = plot_power_spectrum_drw(d,t,p,percentile=95.0, nsamples=2000)
    fmt1 = {'frequencies': '%.8f', 'psd_low': '%.8f', 'psd_hi': '%.8f', 'psd_mid': '%.8f'} 
    result={'frequencies':frequencies,'psd_low': psd_low, 'psd_hi':psd_hi, 'psd_mid':psd_mid }
    ascii.write(result, plot_dir+ 'gw_test', formats=fmt1, names=['frequencies','psd_low','psd_hi','psd_mid'])
    #plot_power_spectrum1(d,t,p,percentile=95.0,nsamples=2000,doShow=True)
    return d
    
   # return 0
df=pd.read_csv(base_dir+'data/CO2_ts.csv',names=['t','y'])
df=df.groupby('t', as_index=False).mean()
#df.drop_duplicates(inplace=True)
t=np.array(df['t']).astype(float)
y=np.array(df['y'])
#idx=np.where(t>100)
#y=y[idx]
#t=t[idx]



for i in xrange(1,len(y)-1):
    if t[i+1]-t[i]==0:
		print t[i]


yerr=np.zeros(len(y))
for i in xrange(len(y)):
	yerr[i]+=0.1
d=fit(t,y,yerr)



   
