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
    
    lc=kali.lc.externalLC('tabby','r',tIn=t,yIn=y,yerrIn=yerr,maskIn=np.ones(len(yerr)))
    ps=[]
    qs=[]
    bics=[]
    timeLCARMA=0
    startLCARMA = time.time()
    taskDict = dict()
    for p in xrange(1,7):
        for q in range(6):
            if p>q:
                try:
                    newTask = kali.carma.CARMATask(p, q, nwalkers=50, nsteps=500)
                    newTask.fit(lc)
                    bics.append(newTask.bicc)
                    ps.append(p)
                    qs.append(q)
                    taskDict['%d %d'%(p, q)] = newTask
                    print "Done with ", p,q
                except Exception:
                    pass
    stopLCARMA = time.time()  
    timeLCARMA = stopLCARMA - startLCARMA
    print timeLCARMA        
    idx=np.argmin(bics)
    pBest=ps[idx]
    qBest=qs[idx]
    df=pd.DataFrame()
    df['p']=ps
    df['q']=qs
    df['bic']=bics
    df.to_csv('tabby_bic_rband.csv')
    print 'Best model for Tabby Star is C-ARMA(%d,%d)'%(pBest, qBest)  
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
df=pd.read_csv(base_dir+'tabby_r_current.csv')
df.drop_duplicates(inplace=True)
df=df.groupby('JD', as_index=False).mean()
t=np.array(df['JD'])
y=np.array(df['mag'])
yerr=np.array(df['err'])

pBest,qBest=fit(t,y,yerr)



   
