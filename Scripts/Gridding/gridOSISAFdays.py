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


rcParams['ytick.major.size'] = 2
rcParams['axes.linewidth'] = .5
rcParams['lines.linewidth'] = .5
rcParams['patch.linewidth'] = .5
rcParams['ytick.labelsize']=8
rcParams['legend.fontsize']=8
rcParams['font.size']=8
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False  )
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )


dataPath = '/Volumes/PETTY_PASSPORT3/DATA/'
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/Drift/OSISAF/'
outPath = '/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/Drifts/OSISAF/'

osisafPath = dataPath+'/ICE_DRIFT/OSISAF/'

dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr

lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

def main(year):
	yearT=year
	yearT2=year

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		#monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
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
			print dayStr

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

			
			print 'Drift, mon:', month, 'day:', xstr1+'-'+xstr2

			fileT=osisafPath+str(yearT2)+'/'+mstr2+'/ice_drift_nh_polstere-625_multi-oi_'+str(yearT)+mstr1+xstr1+'1200-'+str(yearT2)+mstr2+xstr2+'1200.nc'

			print fileT
			print glob(fileT)

			if (size(glob(fileT))>0):
				ux, vy, mag, latsO, lonsO, xptsO, yptsO = cF.getOSISAFDrift(m, fileT)

				#if we want to set masked values back to nan for gridding purposes
				ux[where(ma.getmask(ux))]=np.nan
				vy[where(ma.getmask(vy))]=np.nan
				drift_day_xy=stack((ux, vy))
				#print drift_day_xy.shape

				driftCG = cF.smoothDriftDaily(xptsG, yptsG, xptsO, yptsO, latsO, drift_day_xy, sigmaG=0.75)
				driftCG=driftCG.astype('f2')

			else:
				# just set the daily drift to a masked array (no drifts available)
				driftCG=ma.masked_all((2,xptsG.shape[0], xptsG.shape[1]))
			#drift_day_xy[1] = vy 
			#sqrt(driftKGKx**2+driftKGKy**2)
			cF.plot_CSAT_DRIFT(m, xptsG , yptsG , driftCG[0], driftCG[1], sqrt(driftCG[0]**2+driftCG[1]**2) , out=figpath+'OSISAF'+str(year)+'_d'+dayStr, units_lab='m/s', units_vec=r'm s$^{-1}$',
				minval=0, maxval=0.5, base_mask=1, res=2, vector_val=0.1, year_string=str(yearT)+mstr1+xstr1+'-'+str(year)+mstr2+xstr2, month_string='', extra='',cbar_type='max', cmap=plt.cm.viridis)
				
			
			driftCG.dump(outPath+str(year)+'/OSISAFDriftG'+dxStr+'-'+str(year)+'_d'+dayStr)

#-- run main program
if __name__ == '__main__':
	for year in xrange(2017, 2017+1, 1):
		print year
		main(year)
