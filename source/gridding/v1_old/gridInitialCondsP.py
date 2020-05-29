
""" gridERAIdays.py
	
	Script to grid the ERA-I snowfall data
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: ERA-I snowfall data
	Output: Gridded ERA-I snowfall data

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		10/01/2018: Version 1
"""

import matplotlib
matplotlib.use("AGG")

from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
from matplotlib import rc
from glob import glob
from netCDF4 import Dataset
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import commonFuncs as cF
import os

m = Basemap(projection='npstere',boundinglat=56,lon_0=-45, resolution='l', round=False)

ancDataPath='/data/users/aapetty/Analysis/NESOSIMdev/AncData/'
figpath='/data/users/aapetty/Figures/NESOSIMdev/InitialConds/'
outPath = '/data/users/aapetty/Forcings/temp2m/ERAI/'
concPath = '/data/users/aapetty/Forcings/'
icPath = '/data/users/aapetty/Forcings/InitialConds/'

dx=100000.
dxStr=str(int(dx/1000))+'km'
print(dxStr)

lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

region_mask, xptsI, yptsI = cF.get_region_mask(ancDataPath, m, xypts_return=1)
region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')

varStr='t2m'
#region_maskG=load(outPath+'regionMaskG'+dxStr)


t2mdurGAll=[]
for y in range(1979, 1991+1, 1):
	if (y==1987):
 		continue
	t2mdurGT=load(outPath+'duration/'+varStr+'duration'+dxStr+'-'+str(y)+'v11')
	t2mdurGAll.append(t2mdurGT)

t2mdurGclim=ma.mean(t2mdurGAll, axis=0)

#cF.plotSnow(m, xptsG, yptsG, t2mdurGclim, out=figpath+'/'+varStr+'/'+varStr+'-clim_duration', units_lab=r'max cumulative days > 0', minval=0, maxval=80, base_mask=0, norm=0, cmap_1=cm.viridis)


w99=cF.getWarren(lonG, latG, 7)

#cF.plotSnow(m, xptsG, yptsG, w99, out=figpath+'/'+varStr+'/W99clim', units_lab=r'snowDepth', minval=0, maxval=10, base_mask=0, norm=0, cmap_1=cm.viridis)


for y in range(2018, 2018+1, 1):
	if (y==1987):
 		continue

	t2mdurGT=load(outPath+'duration/'+varStr+'duration'+dxStr+'-'+str(y)+'v11')
	W99yrT=w99*(t2mdurGclim/t2mdurGT)
	W99yrT[where(latG<70)]=0
	W99yrT[where(region_maskG>8.2)]=0
	W99yrT[where(region_maskG<=7.8)]=0
	W99yrT[where(W99yrT<0)]=0
	W99yrT[where(W99yrT>10)]=10

	day=226
	dayStr=str(day) #226 is middle of August

	iceConcDayG=load(concPath+'/IceConc/CDR/'+str(y)+'/iceConcG_CDR'+dxStr+'-'+str(y)+'_d'+dayStr+'v11')
	W99yrT[where(iceConcDayG<0.15)]=0


	cF.plotSnow(m, xptsG, yptsG, W99yrT, out=figpath+'/W99scaled-IC'+str(y)+'v11', units_lab=r'snowDepth', minval=0, maxval=12, base_mask=0, norm=0, cmap_1=cm.viridis)

	W99yrT.dump(icPath+'/ICsnow'+str(y)+'-'+dxStr+'-v11')



	

