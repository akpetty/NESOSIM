
""" gridCDR.py
	
	Script to grid the CDR sea ice concentration data
	Model written by Alek Petty (20/04/2019)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: CDR sea ice concentration data 
	Output: Gridded sea ice concentration data

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		20/04/2019: Version 1
		05/01/2020: Version 2 - new domain and using cartopy/pyproj for gridding
		05/01/2022: Version 3 - added weights for interpolation to speed up processing

	To do:
		Introduce a better pole hole interpolation method

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


from config import cdr_raw_path
from config import forcing_save_path
from config import figure_path


def main(year, startMonth=3, endMonth=3, extraStr='v11_1', dx=100000, data_path=cdr_raw_path, out_path=forcing_save_path, fig_path=figure_path+'IceConc/CDR/', anc_data_path='../../anc_data/'):
		
	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	print(xptsG)
	print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	print(dxStr)

	region_mask, xptsI, yptsI, lonsI, latsI = cF.get_region_mask_pyproj(anc_data_path, proj, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

	product='CDR'

	numDaysYr=cF.getLeapYr(year)
	if (numDaysYr>365):
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
		mstr='%02d' %(month+1)
		# Get pole hole	
		pmask=cF.get_pmask(year, month)

		numDays=monIndex[month]
		
		# should return array with nans not masked as needed for regridding.
		for x in range(numDays):
			dayT=sum(monIndex[0:month])+x
			daySumStr='%03d' %(dayT)
			dayMonStr='%02d' %(x+1)
			print('day month string', dayMonStr)
			
			try:
				# Try final first
				fileT=glob(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
				print(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')
				print('final data')
			except:
				try:
					# Try nrt
					fileT=glob(data_path+'/nrt/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
					print(data_path+'/nrt/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')
					print('nrt data')
				except:	
					try:
						dayMonStr='%03d' %(dayT-1)
						fileT=glob(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
						print(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')
						print('previous day final data')
					except:
						try:
							dayMonStr='%03d' %(dayT+1)
							fileT=glob(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
							print(data_path+'/final/'+str(year)+'/'+mstr+'/*'+str(year)+mstr+dayMonStr+'*.nc')
							print('following day final data')
						except:
							print('no conc')
							# previously was pass but this meant it would just use the previous loop data below
							continue

			iceConcDay = cF.getCDRconcproj(fileT, mask=0, maxConc=1, lowerConc=1)
			print(iceConcDay.shape)

			# if it's the first day, calculate weights
			if calc_weights == 1:

				# calculate Delaunay triangulation interpolation weightings for first file of the year
				print('calculating interpolation weightings')

				# Use region mask which is on the same grid!
				ptM_arr = np.array([xptsI.flatten(),yptsI.flatten()]).T
				tri = Delaunay(ptM_arr) # delaunay triangulation
				calc_weights = 0

			
			iceConcDay[np.where(region_mask>10)]=np.nan

			#iceConcDayG = griddata((xpts0.flatten(), ypts0.flatten()), iceConcDay.flatten(), (xptsG, yptsG), method='linear')
			# grid using linearNDInterpolator with triangulation calculated above 
			# (faster than griddata but produces identical output)
			interp = LinearNDInterpolator(tri,iceConcDay.flatten())
			iceConcDayG = interp((xptsG,yptsG))

			iceConcDayG[np.where(iceConcDayG<0.15)]=0
			iceConcDayG[np.where(iceConcDayG>1)]=1

			iceConcDayG[np.where(region_maskG>10)]=np.nan

			cF.plot_gridded_cartopy(lonG, latG, iceConcDayG, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=fig_path+'iceConcG_-'+str(year)+mstr+dayMonStr+dxStr+extraStr, 
				date_string=str(year), month_string=mstr+dayMonStr, extra=extraStr, varStr='CDR ice conc', units_lab='', minval=0, maxval=1, cmap_1=plt.cm.viridis)
		
			iceConcDayG.dump(out_path+str(year)+'/iceConcG_CDR'+dxStr+'-'+str(year)+'_d'+daySumStr+extraStr)


#-- run main program
if __name__ == '__main__':
	for y in range(2020, 2020+1, 1):
		print(y)
		main(y)


