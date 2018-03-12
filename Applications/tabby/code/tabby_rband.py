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
    
    lc=kali.lc.externalLC('tabby','g',tIn=t,yIn=y,yerrIn=yerr,maskIn=np.ones(len(yerr)))
    ps=[]
    qs=[]
    bics=[]
    timeLCARMA=0
 
    startLCARMA = time.time()
    bestTask = kali.carma.CARMATask(4, 0, nwalkers=50, nsteps=500) 
    bestTask.fit(lc)
    stopLCARMA = time.time()  
    timeLCARMA = stopLCARMA - startLCARMA
    print timeLCARMA      
    medRho=bestTask.medRho
    uRho=bestTask.upperRho
    lRho=bestTask.lowerRho
    sig2=bestTask.Sigma()
    print "medRho is", medRho
    print "upperRho is", uRho
    print "lowerRho is", lRho
    print "sigsqr is", sig2
    print bestTask.logLikelihood(lc)

    psd=bestTask.psd()
    pfile = open(save_dir + 'tabby_rband.pickle', 'wb')
    cPickle.dump(psd, save_dir + 'tabby_rband.pickle', 'wb')
    pfile.close()

    bestTask.plotpsd(doShow = False)
    plt.savefig(save_dir+'tabby_rband_psd.png')
    plt.close()
   
df=pd.read_csv(base_dir+'tabby_r_current.csv')
df.drop_duplicates(inplace=True)
df=df.groupby('JD', as_index=False).mean()
t=np.array(df['JD'])
y=np.array(df['mag'])
yerr=np.array(df['err'])
fit(t,y,yerr)



   
