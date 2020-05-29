
""" gridCDR.py
	
	Script to grid the CDR sea ice concentration data
	Model written by Alek Petty (20/04/2019)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: sea ice concentration data (NASA Team or Bootstrap)
	Output: Gridded sea ice concentration data

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		20/04/2019: Version 1
		05/01/2020: Version 2
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

#dataPath = '/data/users/aapetty/Data/'
#ancDataPath='/data/users/aapetty/Analysis/NESOSIMdev/AncData/'
#figpath='/data/users/aapetty/Figures/NESOSIMdev/IceConc/CDR/'
#outPath = '/data/users/aapetty/Forcings/IceConc/CDR/'


def main(year, startMonth=0, endMonth=11, extraStr='v2', dx=100000, data_path='.', out_path='.', fig_path='.', anc_data_path='../../AncData/'):
	
	m = Basemap(projection='npstere',boundinglat=56,lon_0=-45, resolution='l', round=False)
	dx=100000.
	dxStr=str(int(dx/1000))+'km'
	print dxStr
	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)
	region_mask, xptsI, yptsI = cF.get_region_mask_sect(ancDataPath, m, xypts_return=1)
	latsI, lonsI = cF.get_psnlatslons(ancDataPath)

	product='CDR'


	numDaysYr=cF.getLeapYr(year)
	if (numDaysYr>365):
		monIndex = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	else:
		monIndex = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


	if not os.path.exists(outPath+'/'+str(year)):
		os.makedirs(outPath+'/'+str(year))

	
	for month in xrange(startMonth, endMonth+1):
		print month
		mstr='%02d' %(month+1)
		# Get pole hole	
		pmask=cF.get_pmask(year, month)

		numDays=monIndex[month]
		
		# should return array with nans not masked as needed for regridding.
		for x in xrange(numDays):
			dayT=sum(monIndex[0:month])+x
			daySumStr='%03d' %(dayT)
			dayMonStr='%02d' %(x+1)
			print('day month string', dayMonStr)
			if (year>2017):
				try:
					fileT=glob(dataPath+'ICE_CONC/CDR/nrt/'+str(year)+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
				except:
					try:
						dayMonStr='%03d' %(dayT-1)
						fileT=glob(dataPath+'ICE_CONC/CDR/nrt/'+str(year)+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
					except:
						try:
							dayMonStr='%03d' %(dayT+1)
							fileT=glob(dataPath+'ICE_CONC/CDR/nrt/'+str(year)+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
						except:
							print('no conc')
							pass

				print(fileT)
			else:
				fileT=glob(dataPath+'ICE_CONC/CDR/final/'+str(year)+'/*'+str(year)+mstr+dayMonStr+'*.nc')[0]
				print(fileT)

			iceConcDay, lats0, lons0, xpts0, ypts0= cF.getCDRconc(m, fileT, mask=1, maxConc=1, lowerConc=1)
			print(iceConcDay.shape, lats0.shape)
			concHole=ma.mean(iceConcDay[(lats0>pmask-2.) & (lats0<pmask-1.5)])
			#print concHole
			iceConcDay = where((lats0 >=pmask-2.), concHole, iceConcDay)

			iceConcDay[where(region_mask>18)]=0
			#iceConcDay[where(region_mask>18)]=0
			#print iceConcDay
			#iceConcDay=ma.filled(np.nan)
			
			iceConcDayG = griddata((xpts0.flatten(), ypts0.flatten()), iceConcDay.flatten(), (xptsG, yptsG), method='linear')

			cF.plotSnow(m, xptsG, yptsG, iceConcDayG, out=figpath+'iceConcG_'+str(year)+'_d'+daySumStr+'v56', units_lab=r'm', minval=0, maxval=1, base_mask=0, norm=0, cmap_1=cm.viridis)
			
			iceConcDayG.dump(outPath+'/'+str(year)+'/iceConcG_'+product+dxStr+'-'+str(year)+'_d'+daySumStr+'v56')


#-- run main program
if __name__ == '__main__':
	for y in xrange(2019, 2019+1, 1):
		print y
		main(y)


