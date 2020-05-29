
""" gridERAIsf.py
	
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
		01/05/2020: Version 2
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
#out_path = '/data/users/aapetty/Forcings/Precip/ERAI/'

def main(year, startMonth=0, endMonth=11, dx=100000, extraStr='v2', data_path='.', out_path='.', fig_path='.', anc_data_path='../../AncData/'):

	# Grid projection etc
	m = Basemap(projection='npstere',boundinglat=56,lon_0=-45, resolution='l', round=False)
	
	dxStr=str(int(dx/1000))+'km'
	#print dxStr
	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)
	region_mask, xptsI, yptsI = cF.get_region_mask(anc_data_path, m, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')

	varStr='sf'

	if not os.path.exists(figpath+'/'+varStr+'/'):
		os.makedirs(figpath+'/'+varStr+'/')


	yearT=year

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
	else:
		monIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

	# index starts at 0 (January)

	if not os.path.exists(outPath+varStr+'/'+str(year)):
		os.makedirs(outPath+varStr+'/'+str(year))


	startDay=monIndex[startMonth]

	endDay=monIndex[endMonth+1]

	for dayT in range(startDay, endDay):
		#dayT=day
		print 'Precip day:', dayT

		dayStr='%03d' %dayT

		# Units of kg/m2 per day
		xptsM, yptsM, lonsM, latsM, Precip =cF.get_ERA_precip_days(m, data_path, yearT, dayT, varStr=varStr, extra=extraStr)

		# Grid the data
		PrecipG = griddata((xptsM.flatten(), yptsM.flatten()), Precip.flatten(), (xptsG, yptsG), method='linear')
		PrecipG[where(region_maskG>10)]=0
		PrecipG=PrecipG.astype('f2')
		
		# Plot map as a test
		cF.plotSnow(m, xptsG, yptsG, PrecipG, out=figpath+'/'+varStr+'/'+varStr+'-'+str(yearT)+'_d'+str(dayT)+'v56', units_lab=r'kg/m2', minval=0, maxval=10, base_mask=0, norm=0, cmap_1=cm.viridis)

		PrecipG.dump(outPath+varStr+'/'+str(yearT)+'/'+'ERAI'+varStr+dxStr+'-'+str(yearT)+'_d'+dayStr+'v56')

# Can either run main program here or call main from the run script
if __name__ == '__main__':
	for y in range(2018, 2019+1, 1):
		print y
		main(y)








#lonS, latS=np.meshgrid(*(lon, lat))


#varMYr = varM[numdays]

	

