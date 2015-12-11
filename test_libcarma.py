import math as math
import cmath as cmath
import numpy as np
import random as random
import cffi as cffi
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import gridspec, cm
import matplotlib.cm as colormap
import matplotlib.mlab as mlab
import time
import pdb

from _libcarma import ffi
from mpl_settings import *

goldenRatio=1.61803398875
fhgt=10.0
fwid=fhgt*goldenRatio
dpi = 300

AnnotateXXLarge = 72
AnnotateXLarge = 48
AnnotateLarge = 32
AnnotateMedium = 28
AnnotateSmall = 24
AnnotateXSmall = 20
AnnotateXXSmall = 16

LegendLarge = 24
LegendMedium = 20
LegendSmall = 16

LabelXLarge = 32
LabelLarge = 28
LabelMedium = 24
LabelSmall = 20
LabelXSmall = 16

AxisXXLarge = 32
AxisXLarge = 28
AxisLarge = 24
AxisMedium = 20
AxisSmall = 16
AxisXSmall = 12
AxisXXSmall = 8

normalFontSize=32
smallFontSize=24
footnoteFontSize=20
scriptFontSize=16
tinyFontSize=12

LabelSize = LabelXLarge
AxisSize = AxisLarge
AnnotateSize = AnnotateXLarge
AnnotateSizeAlt = AnnotateMedium
AnnotateSizeAltAlt = AnnotateLarge
LegendSize = LegendMedium

gs = gridspec.GridSpec(1000, 1000) 

set_plot_params(fontfamily='serif',fontstyle='normal',fontvariant='normal',fontweight='normal',fontstretch='normal',fontsize=AxisMedium,useTex='True')

ffiObj = cffi.FFI()
C = ffi.dlopen("./libcarma.so")

dt = 0.01
p = 2
q = 1
Theta = [0.75, 0.01, 7.0e-9, 1.2e-9]
numBurn = 1000000
numCadences = 60000
noiseSigma = 1.0e-18
startCadence = 0
burnSeed = 1311890535
distSeed = 2603023340
noiseSeed = 2410288857
cadence = np.array(numCadences*[0])
mask = np.array(numCadences*[0.0])
t = np.array(numCadences*[0.0])
y = np.array(numCadences*[0.0])
yerr = np.array(numCadences*[0.0])

dt_cffi = dt
p_cffi = p
q_cffi = q
Theta_cffi = ffiObj.new("double[%d]"%(len(Theta)))
for i in xrange(len(Theta)):
	Theta_cffi[i] = Theta[i]
numBurn_cffi = numBurn
numCadences_cffi = numCadences
noiseSigma_cffi = noiseSigma
startCadence_cffi = startCadence
burnSeed_cffi = burnSeed
distSeed_cffi = distSeed
noiseSeed_cffi = noiseSeed
cadence_cffi = ffiObj.new("int[%d]"%(numCadences))
mask_cffi = ffiObj.new("double[%d]"%(numCadences))
t_cffi = ffiObj.new("double[%d]"%(numCadences))
y_cffi = ffiObj.new("double[%d]"%(numCadences))
yerr_cffi = ffiObj.new("double[%d]"%(numCadences))

for i in xrange(numCadences):
	cadence_cffi[i] = i
	mask_cffi[i] = 1.0
	t_cffi[i] = dt_cffi*i
	y_cffi[i] = 0.0
	yerr_cffi[i] = 0.0

makeLCStart = time.time()

YesOrNo = C.cffi_makeMockLC(dt_cffi, p_cffi, q_cffi, Theta_cffi, numBurn_cffi, numCadences_cffi, noiseSigma_cffi, startCadence_cffi, burnSeed_cffi, distSeed_cffi, noiseSeed_cffi, cadence_cffi, mask_cffi, t_cffi, y_cffi, yerr_cffi)

makeLCStop = time.time()

print "Time to make LC: %f (s)"%(makeLCStop - makeLCStart)

computeLnLikeStart = time.time()

LnLike = C.cffi_computeLnLike(dt_cffi, p_cffi, q_cffi, Theta_cffi, numCadences_cffi,cadence_cffi, mask_cffi, t_cffi, y_cffi, yerr_cffi)

computeLnLikeStop = time.time()

print "LnLike: %+16.17e"%(LnLike)

print "Time to compute LnLike: %f (s)"%(computeLnLikeStop - computeLnLikeStart)

for i in xrange(numCadences):
	cadence[i] = cadence_cffi[i]
	mask[i] = mask_cffi[i]
	t[i] = t_cffi[i]
	y[i] = y_cffi[i]
	yerr[i] = yerr_cffi[i]

plt.figure(1,figsize=(fwid,fhgt))
plt.errorbar(t,y,yerr,fmt='.',capsize=0,color='#d95f02',markeredgecolor='none',zorder=10)
yMax=np.max(y[np.nonzero(y[:])])
yMin=np.min(y[np.nonzero(y[:])])
plt.ylabel(r'$F$ (arb units)')
plt.xlabel(r'$t$ (d)')
plt.xlim(t[0],t[-1])
plt.ylim(yMin,yMax)
plt.show()

pdb.set_trace()