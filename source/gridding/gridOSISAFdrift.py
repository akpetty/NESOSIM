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
		05/01/2020: Version 2 - new domain and using cartopy/pyproj for gridding
		05/01/2022: Version 3 - added weights for interpolation to speed up processing
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

from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

from config import osisaf_raw_path
from config import forcing_save_path
from config import figure_path


def main(year, startMonth=0, endMonth=11, extraStr='v11_1', dx=100000, data_path=osisaf_raw_path, out_path=forcing_save_path, fig_path=figure_path+'IceDrift/OSISAF/', anc_data_path='../../AncData/'):

	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)

	#srcProj = cF.OSISAF() 
	#newProj = cF.P3413() 

	print(xptsG)
	print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	print(dxStr)

	yearT=year
	yearT2=year

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	else:
		monIndex = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	if not os.path.exists(out_path+'/'+str(year)):
		os.makedirs(out_path+'/'+str(year))

	if not os.path.exists(fig_path):
		os.makedirs(fig_path)

	calc_weights = 1 # start as one to calculate weightings then gets set as zero for future files

	for month in range(startMonth, endMonth+1):
		print(month)
		numDays=monIndex[month]

		for x in range(numDays):
			dayT=sum(monIndex[0:month])+x
			dayStr='%03d' %dayT
			print(dayStr)

			# reset year
			yearT=year
			yearT2=year

			if (x==0):
				mstr1='%02d' %(month)
				mstr2='%02d' %(month+1)
				xstr1='%02d' %monIndex[month-1]
				xstr2='%02d' %(x+2)
				if (month==0):
					yearT=year-1
					mstr1='%02d' %12

			elif (x==numDays-1):
				mstr1='%02d' %(month+1)
				mstr2='%02d' %(month+2)
				xstr1='%02d' %(x) 
				xstr2='%02d' %1
				if (month==11):
					yearT2=yearT2+1
					mstr2='%02d' %1
			else:
				mstr1='%02d' %(month+1)
				mstr2='%02d' %(month+1)
				xstr1='%02d' %(x) 
				xstr2='%02d' %(x+2) 

			print('Drift, mon:', month, 'day:', xstr1+'-'+xstr2)
			print('Drift, mstr1:', mstr1, 'mstr2:', mstr2)
			print(data_path+str(yearT2)+'/'+mstr2+'/ice_drift_nh_polstere-625_*'+str(yearT)+mstr1+xstr1+'1200-'+str(yearT2)+mstr2+xstr2+'1200.nc')
			try:
				fileT=glob(data_path+str(yearT2)+'/'+mstr2+'/ice_drift_nh_polstere-625_*'+str(yearT)+mstr1+xstr1+'1200-'+str(yearT2)+mstr2+xstr2+'1200.nc')[0]
			except:
				print('no file')
				continue

			print(fileT)

			if (np.size(glob(fileT))>0):
				ux, vy, latsO, lonsO, xptsO, yptsO = cF.get_osisaf_drifts_proj(proj, fileT)

				#Set masked values back to nan for gridding purposes
				ux[np.where(ma.getmask(ux))]=np.nan
				vy[np.where(ma.getmask(vy))]=np.nan
				drift_day_xy=np.stack((ux, vy))

				# if it's the first day, calculate weights
				if calc_weights == 1:
					# calculate Delaunay triangulation interpolation weightings for first file of the year
					print('calculating interpolation weightings')
					print(latsO)
					print(yptsO)
					ptM_arr = np.array([xptsO.flatten(),yptsO.flatten()]).T
					tri = Delaunay(ptM_arr) # delaunay triangulation
					calc_weights = 0


				#drift_xyG = cF.int_smooth_drifts_v2(xptsG, yptsG, xptsO, yptsO, latsO, drift_day_xy, sigma_factor=1)
				drift_xyG = cF.int_smooth_drifts_v3(tri, xptsG, yptsG, xptsO, yptsO, latsO, drift_day_xy, sigma_factor=1)

			#drift_day_xy[1] = vy 
			print(drift_xyG)	

			cF.plot_drift_cartopy(lonG , latG , xptsG, yptsG, drift_xyG[0], drift_xyG[1], np.sqrt(drift_xyG[0]**2+drift_xyG[1]**2) , out=fig_path+str(year)+'_d'+dayStr+dxStr+extraStr, units_lab='m/s', units_vec=r'm s$^{-1}$',
				minval=0, maxval=0.5, vector_val=0.1, date_string=str(yearT)+mstr1+xstr1+'-'+str(year)+mstr2+xstr2, month_string='', varStr='OSI SAF ice drift ',cbar_type='max', cmap_1=plt.cm.viridis)
				
			drift_xyG.dump(out_path+str(year)+'/OSISAF_driftG'+dxStr+'-'+str(year)+'_d'+dayStr+extraStr)

#-- run main program
if __name__ == '__main__':
	
	for year in range(2010, 2010+1, 1):
		print(year)
		main(year)
	


	#years=np.arange(2019, 2019+1, 1)
	#from itertools import repeat
	#import concurrent.futures
	#with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:

		#args=((campaign, beam) for beam in beams)
		#print(args)
		#esult=executor.map(main, years)






