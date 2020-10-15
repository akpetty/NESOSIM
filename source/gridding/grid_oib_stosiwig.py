
""" grid_oib_stosiwig.py
	
	Script to grid the STOSIWIG raw OIB snow picks to the NESOSIM grid
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



def main(year, dx=100000, snowType='SRLD', extraStr='v11'):
	
	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	print(xptsG)
	print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	print(dxStr)

	#region_mask, xptsI, yptsI = cF.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
	#region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

	#xptsDays, yptsDays, oibdates, snowDays= cF.read_icebridge_snowdepths(proj, oib_data_path, year)
	xptsDays, yptsDays, _, _,snowDays, oibdates= cF.getSTOSIWIGyear_proj(proj, oib_data_path+'/stosiwig/', snowType, year)


	for x in range(len(oibdates)):
		# Loop through dates of each flight. I want to keep separate to compre to the daily NESOSIM data.
		oib_dayG=bin_oib(dx, xptsDays[x], yptsDays[x], xptsG, yptsG, snowDays[x])
		cF.plot_gridded_cartopy(lonG, latG, oib_dayG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/OIB/'+oibdates[x]+dxStr+snowType+extraStr, date_string=oibdates[x], month_string='', varStr='OIB snow depth ', units_lab=r'm', minval=0, maxval=0.6, cmap_1=plt.cm.viridis)
			
		oib_dayG.dump(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+snowType+extraStr)

	arr = np.hstack(xptsDays)

	oib_daysG=cF.bin_oib(dx, np.hstack(xptsDays), np.hstack(yptsDays), xptsG, yptsG, np.hstack(snowDays))

	cF.plot_gridded_cartopy(lonG, latG, oib_daysG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/OIB/'+str(year)+dxStr+snowType+extraStr, date_string=str(year), month_string='', varStr='OIB snow depth ', units_lab=r'm', minval=0, maxval=0.6, cmap_1=plt.cm.viridis)
			

#-- run main program
if __name__ == '__main__':
	for y in range(2009, 2016+1, 1):
		print (y)
		main(y)
