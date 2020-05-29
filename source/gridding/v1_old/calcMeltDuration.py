
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

m = Basemap(projection='npstere',boundinglat=56,lon_0=-45, resolution='l', round=False)
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )




#reanalysisDataPath = '/data/users/aapetty/Data/'
#ancDataPath='/data/users/aapetty/Analysis/NESOSIMdev/AncData/'
#figpath='/data/users/aapetty/Figures/NESOSIMdev/ERAI/'
#outPath = '/data/users/aapetty/Forcings/temp2m/ERAI/'


def main(year, startMonth=0, endMonth=11, extraStr='v2', dx=100000, data_path='.', out_path='.', fig_path='.', anc_data_path='../../AncData/'):
	
	reanalysis='ERAI'

	dxStr=str(int(dx/1000))+'km'
	print dxStr

	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

	region_mask, xptsI, yptsI = cF.get_region_mask(anc_data_path, m, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')
	#region_maskG=load(outPath+'regionMaskG'+dxStr)

	varStr='t2m'
	#varStr='tp'


	yearT=year

	xptsM, yptsM, lonsM, latsM, temp2mYear =cF.get_ERA_2mt_meltduration(m, data_path, yearT)
	t2mdur=np.zeros((temp2mYear.shape[1], temp2mYear.shape[2]))

	# temp2mYear is the gridded daily mean temperature data for all days within the year (day, x, y)
	# t2mdur should provide the gridded melt season duration (x, y)


	import itertools
	# loop over each grid-cell
	for i in xrange(t2mdur.shape[0]):
		for j in xrange(t2mdur.shape[1]):
			# Find all the locations where the temps are above freezing for a specific grid cell
			a=where(temp2mYear[:, i, j]>0)[0]
			if (size(a)>1):
				# If there is more than one day when temps are above freezing
				# This function is a bit abstract but it finds the maximum number of cumulative days of above freezing temps within the year (using the fact that the difference between these locations should be 1).
				t2mdur[i, j]=max(len(list(v)) for g,v in itertools.groupby(diff(a)))
				
			else: 
				# If there is only one day or less of above freezing temps then set the duration to one or zero
				t2mdur[i, j]=size(a)



	t2mdurG = griddata((xptsM.flatten(), yptsM.flatten()), t2mdur.flatten(), (xptsG, yptsG), method='linear')
	#PrecipG=PrecipG[0]
	t2mdurG[where(region_maskG>10)]=0.
	#t2mdur=t2mdur.astype('f2')

	cF.plotSnow(m, xptsG, yptsG, t2mdurG, out=fig_path+'/'+varStr+'/'+varStr+'-'+str(yearT)+'_durationv56', units_lab=r'max cumulative days > 0', minval=0, maxval=80, base_mask=0, norm=0, cmap_1=cm.viridis)

	#cF.plotSnow(m, xptsG, yptsG, region_maskG, out=figpath+'/'+varStr+'/region_mask', units_lab=r'region', minval=0, maxval=11, base_mask=0, norm=0, cmap_1=cm.viridis)

	#monthStr='%02d' %(month+1)
	t2mdurG.dump(outPath+'duration/'+varStr+'duration'+dxStr+'-'+str(yearT)+'v56')

#-- run main program
if __name__ == '__main__':
	for y in xrange(2015, 2018+1, 1):
		print y
		main(y)


	

