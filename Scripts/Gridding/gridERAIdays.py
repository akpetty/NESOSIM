
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

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False)
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )

reanalysis='ERAI'

rawDataPath = '/Volumes/PETTY_PASSPORT3/DATA/'
reanalysisDataPath = '/Volumes/PETTY_PASSPORT3/'
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/ERAI/'
outPath = '/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/Precip/ERAI/'

dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr

lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

region_mask, xptsI, yptsI = cF.get_region_mask(rawDataPath, m, xypts_return=1)

region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')
region_maskG.dump(outPath+'regionMaskG'+dxStr)

varStr='sf'
#varStr='tp'

if not os.path.exists(figpath+'/'+varStr+'/'):
	os.makedirs(figpath+'/'+varStr+'/')

def main(year):
	yearT=year

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
	else:
		monIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

	# FOR ALL MONTHS ADD 1, SO 9 iS OCTOBER
	startMonth=0
	endMonth=11


	if not os.path.exists(outPath+varStr+'/'+str(year)):
		os.makedirs(outPath+varStr+'/'+str(year))


	startDay=monIndex[startMonth]
	if (endMonth>11):
		endDay=monIndex[endMonth+1-12]+monIndex[-1]-1
	else:
		endDay=monIndex[endMonth+1]

	#glob(dataPath+'/REANALYSES/MERRA2-precip/PrecipTotal_'+str(year)+str(day))

	for day in xrange(startDay, endDay):
		dayT=day
		print 'Precip day:', day
		if (day>=numDays):
			dayT=day-numDays
			yearT=year+1
		dayStr='%03d' %dayT


		#in  kg/m2 per day
		xptsM, yptsM, lonsM, latsM, Precip =cF.get_ERA_precip_days(m, reanalysisDataPath, yearT, dayT, varStr=varStr)

		#xptsM, yptsM, lonsM, latsM, Precip =pF.get_ERA_precip_days2000to2017(m, reanalysisDataPath, yearT, dayT, varStr=varStr)


		PrecipG = griddata((xptsM.flatten(), yptsM.flatten()), Precip.flatten(), (xptsG, yptsG), method='linear')
		#PrecipG=PrecipG[0]
		PrecipG[where(region_maskG>10)]=0
		PrecipG=PrecipG.astype('f2')
		
		cF.plotSnow(m, xptsG, yptsG, PrecipG, out=figpath+'/'+varStr+'/'+varStr+'-'+str(yearT)+'_d'+str(dayT)+'T2', units_lab=r'kg/m2', minval=0, maxval=10, base_mask=0, norm=0, cmap_1=cm.viridis)

		#monthStr='%02d' %(month+1)
		PrecipG.dump(outPath+varStr+'/'+str(yearT)+'/'+reanalysis+varStr+dxStr+'-'+str(yearT)+'_d'+dayStr)

	#PrecipDaysG.dump(outPath+reanalysis+varStr+dxStr+'-'+str(year)+str(startMonth)+'-'+str(yearT)+str(endMonth))

#-- run main program
if __name__ == '__main__':
	for y in xrange(2017, 2017+1, 1):
		print y
		main(y)








#lonS, latS=np.meshgrid(*(lon, lat))


#varMYr = varM[numdays]

	

