""" gridOSISAFdrift.py
	
	Script to grid the OSISAF derived Arctic ice drifts
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: OSISAF ice drifts
	Output: Gridded OSISAF ice drifts

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		10/01/2018: Version 1
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import utils as cF
import os
import pyproj
import cartopy.crs as ccrs
import xarray as xr
from config import nsidc_raw_path
from config import forcing_save_path
from config import figure_path

xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=25000.)

fileT=nsidc_raw_path+'/weekly/icemotion_weekly_nh_25km_20180101_20181231_v4.1.nc'
nsidc_drifts_weekly= xr.open_dataset(fileT)

# combine the weekly data across years so when we interpolate it fills in the end

nsidc_drifts_monthly=nsidc_drifts_weekly.resample(time='MS').mean()

nsidc_drifts_daily=nsidc_drifts_weekly.resample(time="1D").interpolate("linear")


xtest = nsidc_drifts_daily.u[0].values
ytest = nsidc_drifts_daily.v[0].values

# convert the drift vectors to zonal/meridional (luckily this dosn't need the original map projection, at least i don't think it does)
lontest=nsidc_drifts_daily.longitude[0].values
lattest=nsidc_drifts_daily.latitude[0].values
alpha = lontest*np.pi/180.
uvelT = ytest*np.sin(alpha) + xtest*np.cos(alpha)
vvelT = ytest*np.cos(alpha) - xtest*np.sin(alpha) 

xdrift, ydrift=proj.rotate_vector(uvelT, vvelT, lontest, lattest)
	# x=0
	# for file in files:

	# 	uvelT, vvelT=getFowlerDrift(file,lonsF)

	# 	xvel,yvel = m.rotate_vector(uvelT,vvelT,lonsF,latsF)

	# 	xvel[where(ma.getmask(xvel))]=np.nan
	# 	yvel[where(ma.getmask(yvel))]=np.nan
	# 	driftFmon[x, 0]=xvel
	# 	driftFmon[x, 1]=yvel
	# 	x+=1

	# if (mean==1):
	# 	driftFmon=ma.mean(driftFmon, axis=0)

	# return xptsF, yptsF, driftFmon, lonsF, latsF