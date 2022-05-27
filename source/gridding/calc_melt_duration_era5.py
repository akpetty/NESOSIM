""" calc_melt_duration_era5.py
	
	Calculate the melt season duration for estimating summer snow conditions.
	Model written by Alek Petty (06/01/2020)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: ERA5 2m air temperatures
	Output: Gridded ERA5 melt duration

	Python dependencies:
		See below for the relevant module imports
		Also some function in utils.py

	Update history:
		05/01/2020: Version 1 (adapted from the earlier ERA-I script)
						- utilize xarray's resample and reduce functions to optimize the code.
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
import itertools
import xarray as xr

from config import reanalysis_raw_path
from config import forcing_save_path
from config import figure_path

def get_ERA5_temps(proj, data_pathT, yearT):
	
	tempdata=xr.open_mfdataset(data_pathT+'/ERA5_temp6hour'+'_'+str(yearT)+'*cds.nc')
		
	numDaysYearT=np.size(tempdata['time'][:])/4
	print (numDaysYearT)

	lon = tempdata['longitude'][:]
	lowerlat=20
	tempdata=tempdata.where(tempdata.latitude>lowerlat, drop=True)
	
	lon = tempdata['longitude'][:]
	lat = tempdata['latitude'][:]
	tempdata_daily=tempdata.resample(time='1d').mean()-273.15
	xpts, ypts=proj(*np.meshgrid(lon, lat))

	return xpts, ypts, tempdata_daily

def get_cumulative_hotdays(temp_year_gridcell):
	
	a=np.where(temp_year_gridcell>0)[0]
	#print(a)
	
	if (np.size(a)>1):
		# If there is more than one day when temps are above freezing
		# This function is a bit abstract but it finds the maximum number of cumulative days of above freezing temps within the year (using the fact that the difference between these locations should be 1).
		t2mdur=max(len(list(v)) for g,v in itertools.groupby(np.diff(a)))			
	else: 
		# If there is only one day or less of above freezing temps then set the duration to one or zero
		t2mdur=np.size(a)
	
	return t2mdur

def main(yearT, extraStr='v11_1', dx=100000, data_path=reanalysis_raw_path+'ERA5/', out_path=forcing_save_path+'Temp/ERA5/', fig_path=figure_path+'Temp/ERA5/', anc_data_path='../../anc_data/'):
	

	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	print(xptsG)
	print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	print(dxStr)

	if not os.path.exists(out_path):
		os.makedirs(out_path)

	if not os.path.exists(fig_path):
		os.makedirs(fig_path)

	region_mask, xptsI, yptsI, _, _ = cF.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

	varStr='t2m'

	xptsM, yptsM, temp2mYear =get_ERA5_temps(proj, data_path, yearT)

	hotdays_duration=xr.apply_ufunc(get_cumulative_hotdays, temp2mYear['t2m'].compute(), input_core_dims=[["time"]], vectorize=True)

	print(xptsM.flatten().shape)
	print(yptsM.flatten().shape)
	print(temp2mYear)
	print(hotdays_duration.values.flatten().shape)

	t2mdurG = griddata((xptsM.flatten(), yptsM.flatten()), hotdays_duration.values.flatten(), (xptsG, yptsG), method='linear')
	#t2mdurG[where(region_maskG>10)]=0.
	#t2mdur=t2mdur.astype('f2')

	cF.plot_gridded_cartopy(lonG, latG, t2mdurG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=fig_path+'/duration'+str(yearT)+extraStr, date_string=str(yearT), extra=extraStr, varStr='hot days', units_lab=r'>0', minval=0, maxval=100, cmap_1=plt.cm.viridis)
		
	#monthStr='%02d' %(month+1)
	t2mdurG.dump(out_path+'/ERA5_duration'+dxStr+'-'+str(yearT)+extraStr)

#-- run main program
if __name__ == '__main__':
	for y in range(2021, 2021+1, 1):
		print(y)
		main(y)


	

