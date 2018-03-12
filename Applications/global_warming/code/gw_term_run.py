import numpy as np
import matplotlib.pyplot as plt
import cPickle
import random
import os
import sys
import pandas as pd
import time
import kali.carma 
import kali
import time
import kali.util.mcmcviz as mcmcviz
import kali.util.triangle as triangle

#Initialization
base_dir='/home/jtrichm2/long_time_series/'
plot_dir=base_dir
data_dir=base_dir
outdir=base_dir

save_dir=base_dir
outfile_pi=base_dir

df=pd.read_csv(base_dir+'data/CO2_ts.csv',names=['t','y'])
df=df.groupby('t', as_index=False).mean()
#df.drop_duplicates(inplace=True)
t=np.array(df['t']).astype(float)
y=np.array(df['y'])
'''
idx=np.where(t>100)
y=y[idx]
t=t[idx]
'''
for i in xrange(1,len(y)-1):
    if t[i+1]-t[i]==0:
		print t[i]
yerr=np.zeros(len(y))
for i in xrange(len(y)):
	yerr[i]+=0.1

pBest=3;qBest=0
nwalkers_=500
nsteps_=200
lc=kali.lc.externalLC('gw','r',tIn=t,yIn=y,yerrIn=yerr,maskIn=np.ones(len(yerr)))

newTask = kali.carma.CARMATask(pBest, qBest, nwalkers=nwalkers_, nsteps=nsteps_)
newTask.fit(lc)
medRho=newTask.medRho
uRho=newTask.upperRho
lRho=newTask.lowerRho
sig2=newTask.Sigma()
print "medRho is", medRho
print "upperRho is", uRho
print "lowerRho is", lRho
print "sigsqr is", sig2

bestLabelList = list()
for i in range(pBest):
    bestLabelList.append("$a_{%d}$"%(i + 1))
for i in range(qBest + 1):
    bestLabelList.append("$b_{%d}$"%(i))
bicc=newTask.bicc
bestFigTitle = 'Best Model for Global CO2 in Atmosphere: CARMA(%d,%d); BIC: %+4.3e'%(
            pBest, qBest, bicc)

mcmcviz.vizTriangle(pBest, qBest, newTask.Chain,labelList=bestLabelList, figTitle=bestFigTitle)
plt.savefig(base_dir+'CARMA('+str(pBest)+','+str(qBest)+ ') fit for gw_no_current.png')
plt.close()

newTask.plotpsd(doShow = False)
plt.savefig(base_dir+'CARMA('+str(pBest)+','+str(qBest)+ ') fit for gw_no_current_psd.png')
plt.close()
#mcmcviz.vizWalkers()
plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[0,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[0,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'AR0.png')
plt.close()

plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[1,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[1,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'AR1.png')
plt.close()

plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[2,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[2,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'AR2.png')
plt.close()
plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[3,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[3,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'AR3.png')
plt.close()
'''
plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[4,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[4,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'MA0.png')
plt.close()

plt.figure()
for i in xrange(nwalkers_):
    plt.plot(newTask.timescaleChain[5,i,:], c = '#000000', alpha = 0.25)
plt.plot(np.median(newTask.timescaleChain[5,:,],axis = 0), c = '#00ffff')
plt.xlim(0,100)
plt.savefig(base_dir+'MA1.png')
plt.close()
'''
