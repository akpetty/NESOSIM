""" gridERAIWinddays.py
	
	Script to grid the ERA-I wind data
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: ERA-I wind data
	Output: Gridded ERA-I wind data

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




#dataPath = '/data/users/aapetty/Data/'
#ancDataPath='/data/users/aapetty/Analysis/NESOSIMdev/AncData/'
#figpath='/data/users/aapetty/Figures/NESOSIMdev/ERAI/'
#outPath = '/data/users/aapetty/Forcings/Winds/ERAI/'



def main(year, startMonth=0, endMonth=11, extraStr='v2', dx=100000, data_path='.', out_path='.', fig_path='.', anc_data_path='../../AncData/'):
	

	reanalysis='ERAI'

	m = Basemap(projection='npstere',boundinglat=56,lon_0=-45, resolution='l', round=False)
	dx=100000.
	dxStr=str(int(dx/1000))+'km'
	print dxStr
	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)
	region_mask, xptsI, yptsI = cF.get_region_mask(ancDataPath, m, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')
	#region_maskG.dump(outPath+'regionMaskG'+dxStr)
	#region_maskG=load(outPath+'regionMaskG'+dxStr)

	yearT=year

	if not os.path.exists(figpath):
		os.makedirs(figpath)

	if not os.path.exists(outPath+str(year)):
		os.makedirs(outPath+str(year))

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
	else:
		monIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]


	# Month index starts at 0
	#startMonth=2
	#endMonth=11

	#varStr='sf'
	varStr='WindMag'

	WindMagDaysG=ma.masked_all((0, nx, ny))

	if not os.path.exists(outPath+'/'+str(year)):
		os.makedirs(outPath+'/'+str(year))
	if not os.path.exists(figpath+'/'+varStr+'/'):
		os.makedirs(figpath+'/'+varStr+'/')

	startDay=monIndex[startMonth]
	if (endMonth>11):
		endDay=monIndex[endMonth+1-12]+monIndex[-1]-1
	else:
		endDay=monIndex[endMonth+1]


	for day in xrange(startDay, endDay):
		dayT=day
		print 'Precip day:', day
		if (day>=numDays):
			dayT=day-numDays
			yearT=year+1
		dayStr='%03d' %dayT

		#print xptsM
		#in  kg/m2 per day
		xptsM, yptsM, lonsM, latsM, WindMag =cF.get_ERA_wind_days(m, dataPath, yearT, dayT, extra=extraStr)

		WindMagG = griddata((xptsM.flatten(), yptsM.flatten()), WindMag.flatten(), (xptsG, yptsG), method='linear')
		WindMagG[where(region_maskG>10)]=0

		# plot map as a check
		cF.plotSnow(m, xptsG, yptsG, WindMagG, out=figpath+'/'+varStr+'/'+varStr+'-'+str(yearT)+'_d'+str(dayT)+'v56', units_lab=r'm/s', minval=0, maxval=10, base_mask=0, norm=0, cmap_1=cm.viridis)


		WindMagG.dump(outPath+'/'+str(yearT)+'/'+reanalysis+varStr+dxStr+str(yearT)+'_d'+dayStr+'v56')


#-- run main program
if __name__ == '__main__':
	for y in xrange(2019, 2019+1, 1):
		print y
		main(y)










#lonS, latS=np.meshgrid(*(lon, lat))


#varMYr = varM[numdays]

	

