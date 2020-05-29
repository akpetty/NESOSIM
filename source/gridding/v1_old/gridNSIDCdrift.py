""" gridNSIDCdrift.py
	
	Script to grid the OSISAF derived Arctic ice drifts
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: NSIDC (Polar Pathfinder) ice drifts
	Output: Gridded NSIDC ice drifts

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


dataPath = '/Volumes/PETTY_PASSPORT3/DATA/'
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/Drift/NSIDCv3/'
outPath = '/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/Drifts/NSIDCv3/'

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False)
dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr
lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

def main(year):
	yearT=year

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
	else:
		monIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]


	if not os.path.exists(outPath+str(year)):
		os.makedirs(outPath+str(year))

	startMonth=0
	endMonth=11

	for month in xrange(startMonth, endMonth+1):
		print month

		if (month>11):
			month=month-12
			yearT=year+1

		mstr = '%02d' %(month+1)
		dateStr=str(yearT)+mstr
		print dateStr

	
		xptsF, yptsF, driftFmon, lonsF, latsF = cF.getFowlerdriftMonthV3(dataPath, yearT, month, m)
			
		numDays=driftFmon.shape[0]
		#driftFmon=ma.masked_invalid(driftFmon)

		# should return array with nans not masked as needed for regridding.
		for x in xrange(numDays):
			xstr='%03d' %x
			print 'Drift, mon:', month, 'day:', xstr

			dayT=monIndex[month]+x
			dayStr='%03d' %dayT

			driftFG = cF.smoothDriftDaily(xptsG, yptsG, xptsF, yptsF, latsF, driftFmon[x], sigmaG=0.75)

			cF.plot_CSAT_DRIFT(m, xptsG , yptsG , driftFG[0], driftFG[1], sqrt(driftFG[0]**2+driftFG[1]**2), out=figpath+'NSIDCv3'+dateStr+'_d'+xstr, units_lab='m/s', units_vec=r'm s$^{-1}$',
				minval=0, maxval=0.5, res=2, vector_val=0.1, year_string=dateStr, month_string=xstr, extra='',cbar_type='max', cmap=plt.cm.viridis)

			driftFG=driftFG.astype('f2')

			driftFG.dump(outPath+str(yearT)+'/NSIDCv3DriftG'+dxStr+'-'+str(yearT)+'_d'+dayStr)



#-- run main program
if __name__ == '__main__':
	for y in xrange(2016, 2016+1, 1):
		print y
		main(y)


