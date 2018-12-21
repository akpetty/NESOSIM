""" gridKIMURAdrift.py
	
	Script to grid the KIMURA derived Arctic ice drifts
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: KIMURA ice drifts
	Output: Gridded KIMURA ice drifts

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
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/Drift/Kimura/'
outPath = '/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/Drifts/Kimura/'

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False  )
dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr
lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)


def main(year):

	numDaysYr=cF.getLeapYr(year)
	if (numDaysYr>365):
		monIndex = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	else:
		monIndex = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	if not os.path.exists(outPath+str(year)):
		os.makedirs(outPath+str(year))

	startMonth=0
	endMonth=11

	for month in xrange(startMonth, endMonth+1):
		print month
		numDays=monIndex[month]

		for x in xrange(numDays):
			dayT=sum(monIndex[0:month])+x
			dayStr='%03d' %dayT

			mstr = '%02d' %(month+1)
			dstr='%02d' %(x+1)

			print dayStr

			if (year>=2015):
				ystr = str(year)[2:]

				filesT=glob(dataPath+'/ICE_DRIFT/KIMURA/'+str(year)+'/'+ystr+mstr+dstr+'*')

				print filesT
				print size(filesT)

				if (size(filesT)>0):
					xptsK, yptsK, driftKday, lonsK, latsK = cF.getKimuradriftDayRaw(dataPath, filesT[0], m)

					driftCG = cF.smoothDriftDaily(xptsG, yptsG, xptsK, yptsK, latsK, driftKday, sigmaG=0.75)
					driftCG=driftCG.astype('f2')

				else:
					print 'NO DRIFT THIS DAY'
					# just set the daily drift to a masked array (no drifts available)
					driftCG=ma.masked_all((2,xptsG.shape[0], xptsG.shape[1]))

			else:
				ystr = str(year)
				filesT=glob(dataPath+'/CURRENTS/Kimura_drift//kimura_drift_NH_'+ystr+mstr+dstr+'.xy')

				print filesT
				print size(filesT)

				if (size(filesT)>0):
					xptsK, yptsK, driftKday, lonsK, latsK = cF.getKimuradriftDayC(dataPath, filesT[0], m)

					driftCG = cF.smoothDriftDaily(xptsG, yptsG, xptsK, yptsK, latsK, driftKday, sigmaG=0.75)
					driftCG=driftCG.astype('f2')

				else:
					print 'NO DRIFT THIS DAY'
					# just set the daily drift to a masked array (no drifts available)
					driftCG=ma.masked_all((2,xptsG.shape[0], xptsG.shape[1]))

			cF.plot_CSAT_DRIFT(m, xptsG , yptsG , driftCG[0], driftCG[1], sqrt(driftCG[0]**2+driftCG[1]**2) , out=figpath+'Kimura'+str(year)+'_d'+dayStr+'N', units_lab='m/s', units_vec=r'm s$^{-1}$',
				minval=0, maxval=0.5, base_mask=1, res=2, vector_val=0.1, year_string=ystr+mstr+dstr, month_string='', extra='',cbar_type='max', cmap=plt.cm.viridis)
				
			
			driftCG.dump(outPath+str(year)+'/KimuraDriftG'+dxStr+'-'+str(year)+'_d'+dayStr)

#-- run main program
if __name__ == '__main__':
	for year in xrange(2010, 2016+1, 1):
		print year
		main(year)
