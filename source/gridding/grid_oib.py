
""" grid_OIB.py
	
	Script to grid the OIB snow depths to the NESOSIM grid
	Model code written by Alek Petty (05/01/2020)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: OIB quicklook data
	Output: Gridded OIB snow depths

	Python dependencies:
		See below for the relevant module imports
		Also reads in some functions in utils.py

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

from config import forcing_save_path, figure_path, oib_data_path, anc_data_path



dx=100000
year=2019
extraStr='v11'
xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
print(xptsG)
print(yptsG)

dxStr=str(int(dx/1000))+'km'
print(dxStr)

#region_mask, xptsI, yptsI = cF.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
#region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

xptsDays, yptsDays,oibdates, snowDays= cF.read_icebridge_snowdepths(proj, oib_data_path, year)

for x in range(len(oibdates)):
	# Loop through dates of each flight. I want to keep separate to compre to the daily NESOSIM data.
	oib_dayG=cF.bin_oib(xptsDays[x], yptsDays[x], xptsG, yptsG, snowDays[x])

	cF.plot_gridded_cartopy(lonG, latG, oib_dayG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/OIB/'+oibdates[x]+dxStr+extraStr, date_string=oibdates[x], month_string='', varStr='OIB snow depth ', units_lab=r'm', minval=0, maxval=0.6, cmap_1=plt.cm.viridis)
		
	oib_dayG.dump(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+extraStr)

arr = np.hstack(xptsDays)

oib_daysG=bin_oib(np.hstack(xptsDays), np.hstack(yptsDays), xptsG, yptsG, np.hstack(snowDays))

cF.plot_gridded_cartopy(lonG, latG, oib_daysG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/OIB/'+str(year)+dxStr+extraStr, date_string=str(year), month_string='', varStr='OIB snow depth ', units_lab=r'm', minval=0, maxval=0.6, cmap_1=plt.cm.viridis)
		


