
""" gridInitialConditions.py
	
	Script to create initial (August) snow depths from reanalysis temperatures/warren scaling
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: Cumulative hot (> 0 C) days
	Output: Initial (August) snow depth

	Python dependencies:
		See below for the relevant module imports
		Also some function in utils.py

	Update history:
		10/01/2018: Version 1
		06/01/2020: Version 2
		 - New projection/pyproj over basemap.
		 - Limit snow to 10 cm maximum depth.
		 - Added guassian smoothing (2 sigma).
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
from scipy.ndimage.filters import gaussian_filter

from config import reanalysis_raw_path
from config import forcing_save_path
from config import figure_path


def main(yearT, extraStr='v11_1', dx=100000, data_path=reanalysis_raw_path+'ERA5/', out_path=forcing_save_path+'InitialConditions/ERA5/', fig_path=figure_path+'InitialConditions/ERA5/', anc_data_path='../../anc_data/'):
	

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

	reanalysis='ERA5'
	varStr='t2m'

	iceconc_path=forcing_save_path+'/IceConc/CDR/'
	temp_path=forcing_save_path+'/Temp/ERA5/'


	t2mdurGAll=[]
	for y in range(1980, 1991+1, 1):
		if (y==1987):
	 		continue
		t2mdurGT=np.load(temp_path+'duration'+dxStr+'-'+str(y)+extraStr, allow_pickle=True)
		t2mdurGAll.append(t2mdurGT)

	t2mdurGclim=ma.mean(t2mdurGAll, axis=0)
	print(t2mdurGclim.shape)
	w99=cF.getWarren(lonG, latG, 7)


	for yearT in range(2021, 2021+1, 1):
		if (yearT==1987):
	 		continue

		t2mdurGT=np.load(temp_path+'duration'+dxStr+'-'+str(yearT)+extraStr, allow_pickle=True)
		print(t2mdurGT.shape)
		W99yrT=w99*(t2mdurGclim/t2mdurGT)
		W99yrT[np.where(latG<70)]=0
		W99yrT[np.where(region_maskG>8.2)]=0
		W99yrT[np.where(region_maskG<=7.8)]=0
		W99yrT[np.where(W99yrT<0)]=0
		W99yrT[np.where(W99yrT>10)]=10

		day=226
		dayStr=str(day) #226 is middle of August

		iceConcDayG=np.load(iceconc_path+str(yearT)+'/iceConcG_CDR'+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr, allow_pickle=True)
		W99yrT[np.where(iceConcDayG<0.15)]=0

		W99yrT = gaussian_filter(W99yrT, sigma=1)

		# Convert to meters
		W99yrT=W99yrT/100.

		cF.plot_gridded_cartopy(lonG, latG, W99yrT, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=fig_path+'/initial_conditions'+str(yearT)+dxStr+extraStr, date_string=str(yearT), extra=extraStr, varStr='Snow depth ', units_lab=r'm', minval=0, maxval=0.12, cmap_1=plt.cm.viridis)
			
		W99yrT.dump(out_path+'ICsnow'+dxStr+'-'+str(yearT)+extraStr)



	

