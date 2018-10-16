
""" gridIceConcDays.py
	
	Script to grid the sea ice concentration data
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: sea ice concentration data (NASA Team or Bootstrap)
	Output: Gridded sea ice concentration data

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

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False)
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )


fataPath = '/Volumes/PETTY_PASSPORT3/DATA/'
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/Diagnostic/IceConc/'
outPath = '/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/IceConc/'


dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr

lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

region_mask, xptsI, yptsI = cF.get_region_mask_sect(dataPath, m, xypts_return=1)

latsI, lonsI = cF.get_psnlatslons(dataPath)

def main(year):

	numDaysYr=cF.getLeapYr(year)
	if (numDaysYr>365):
		#monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
		monIndex = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	else:
		monIndex = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	alg=1

	if (alg==0):
		team_s = 'nt'
		vStr='v1.1'
	if (alg==1):
		team_s = 'bt'
		vStr='v03'

	if not os.path.exists(outPath+team_s+'/'+str(year)):
		os.makedirs(outPath+team_s+'/'+str(year))

	startMonth=0
	endMonth=11

	for month in xrange(startMonth, endMonth+1):
		print month

		# Get pole hole	
		pmask=cF.get_pmask(year, month)

		numDays=monIndex[month]
		# should return array with nans not masked as needed for regridding.
		for x in xrange(numDays):
			dayT=sum(monIndex[0:month])+x
			dayStr='%03d' %dayT

			if (year>2017):
				iceConcDay = cF.get_day_concSN_NRT(dataPath, year, month, x, alg=alg, vStr=vStr, pole='A',lowerConc=1, maxConc=1,mask=1)

			else:
				iceConcDay = cF.get_day_concSN_daily(dataPath, year, month, x, alg=alg, vStr=vStr, pole='A',lowerConc=1, maxConc=1,mask=1)

			concHole=ma.mean(iceConcDay[(latsI>pmask-0.5) & (latsI<pmask)])

			#iceConcMon[x] = where((latsI >=pmask-0.5), 1, iceConcMon[x])
			iceConcDay = where((latsI >=pmask-0.5), concHole, iceConcDay)
			#iceConcMon[x]=ma.filled(np.nan)

			iceConcDay[where(region_mask>18)]=0
			iceConcDayG = griddata((xptsI.flatten(), yptsI.flatten()), iceConcDay.flatten(), (xptsG, yptsG), method='linear')

			cF.plotSnow(m, xptsG, yptsG, iceConcDayG, out=figpath+'iceConcG_'+team_s+str(year)+'_d'+dayStr, units_lab=r'm', minval=-1, maxval=1, base_mask=0, norm=0, cmap_1=cm.viridis)
			
			iceConcDayG.dump(outPath+team_s+'/'+str(year)+'/iceConcG_'+team_s+dxStr+'-'+str(year)+'_d'+dayStr)


#-- run main program
if __name__ == '__main__':
	for y in xrange(2016, 2017+1, 1):
		print y
		main(y)


