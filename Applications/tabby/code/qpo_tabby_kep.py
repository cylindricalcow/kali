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
    bestTask = kali.carma.CARMATask(1, 0, nwalkers=50, nsteps=500) #2,1
    bestTask.fit(lc)
    stopLCARMA = time.time()  
    timeLCARMA = stopLCARMA - startLCARMA
    print timeLCARMA   
    medTau=bestTask.medTau
    uTau=bestTask.upperTau
    lTau=bestTask.lowerTau
    sig2=bestTask.Sigma()
    print "medTau is", medTau
    print "upperTau is", uTau
    print "lowerTau is", lTau
    print "sigsqr is", sig2 
    '''    
    medRho=bestTask.medRho
    uRho=bestTask.upperRho
    lRho=bestTask.lowerRho
    sig2=bestTask.Sigma()
    print "medRho is", medRho
    print "upperRho is", uRho
    print "lowerRho is", lRho
    print "sigsqr is", sig2
    '''
    print bestTask.logLikelihood(lc)
df=pd.read_csv(base_dir+'taby.csv')
t=np.array(df['JD'])
y=np.array(df['mag'])
yerr=np.array(df['err'])
fit(t,y,yerr)



   
