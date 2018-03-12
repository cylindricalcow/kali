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

#Initialization
base_dir='/home/jtrichm2/long_time_series/'
plot_dir=base_dir
data_dir=base_dir
outdir=base_dir

save_dir=base_dir
outfile_pi=base_dir

def fit(t,y,yerr):#,new_dir): 
    #Read in light curve of s82_id
    
    lc=kali.lc.externalLC('gw','r',tIn=t,yIn=y,yerrIn=yerr,maskIn=np.ones(len(yerr)))
    ps=[]
    qs=[]
    aics=[]
    timeLCARMA=0
    startLCARMA = time.time()
    taskDict = dict()
    for p in xrange(1,7):
        for q in xrange(6):
            if p>q:
                try:
                    newTask = kali.carma.CARMATask(p, q, nwalkers=100, nsteps=500)
                    newTask.fit(lc)
                    
                    aics.append(newTask.aicc)
                    ps.append(p)
                    qs.append(q)
                    taskDict['%d %d'%(p, q)] = newTask
                    print "Done with ", p,q
                except Exception:
                    pass
    stopLCARMA = time.time()  
    timeLCARMA = stopLCARMA - startLCARMA
    print timeLCARMA        
    idx=np.argmin(aics)
    pBest=ps[idx]
    qBest=qs[idx]
    df=pd.DataFrame()
    df['p']=ps
    df['q']=qs
    df['aic']=aics
    df.to_csv('gw_aicc.csv')
    print 'Best model for Global CO2 in atmosphere is C-ARMA(%d,%d)'%(pBest, qBest)  
    bestTask = taskDict['%d %d'%(pBest, qBest)]
    medRho=bestTask.medRho
    uRho=bestTask.upperRho
    lRho=bestTask.lowerRho
    sig2=bestTask.Sigma()
    print "medRho is", medRho
    print "upperRho is", uRho
    print "lowerRho is", lRho
    print "sigsqr is", sig2
    return pBest,qBest
df=pd.read_csv(base_dir+'data/CO2_ts.csv',names=['t','y'])
#df.drop_duplicates(inplace=True)
t=np.array(df['t'])
y=np.array(df['y'])
idx=np.where(t>100)
y=y[idx]
t=t[idx]
y=y/10

yerr=np.zeros(len(y))
for i in xrange(len(y)):
	yerr[i]+=0.1
pBest,qBest=fit(t,y,yerr)



   
