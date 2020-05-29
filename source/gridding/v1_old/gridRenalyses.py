
############################################################## 
# Date: 01/01/17
# Name: gridERAI.py
# Author: Alek Petty

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
sys.path.append('../../../common/')
sys.path.append('../')
import commonFuncs as cF
import pastaFunctions as pF
import os

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False)
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )

#reanalysis='NCEP_R1' # MERRA, MERRA_2, ASR, JRA55, CFSR, NCEP_R1, NCEP_R2, NCEP_CIRES

rawDataPath = '../../../../../DATA/'
#reanalysisDataPath = '/Volumes/PETTY_PASSPORT2/'
#reanalysisDataPath ='smb://gs615-oibserve.ndc.nasa.gov/icebridgedata2/data/Linette_otherdrive/'
reanalysisDataPath ='/Volumes/icebridgedata2/data/Linette_otherdrive/'


dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr

lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)

region_mask, xptsI, yptsI = cF.get_region_mask(rawDataPath, m, xypts_return=1)

region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='linear')
#region_maskG=load(outPath+'regionMaskG'+dxStr)

def main(year, reanalysis):
	print year, reanalysis
	yearT=year

	numDays=365
	monIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]


	outPath='../../DataOutput/Reanalyses/'+reanalysis+'/'
	figpath = '/Volumes/PETTY_PASSPORT2/pasta/Figures/Budgets/'+reanalysis+'/'
	
	# Index starts at 0 (January)
	startMonth=0
	endMonth=11

	varStr='sf'
	#varStr='tp'

	if (varStr=='sf'):
		if (reanalysis=='NCEP_R1' or reanalysis=='NCEP_R2' or reanalysis=='NCEP_CIRES'):
			precipStr='SnowfallRate'
		elif (reanalysis=='MERRA' or reanalysis=='MERRA2' or reanalysis=='MERRA_2'):
			precipStr='Snowfall' 
		elif (reanalysis=='JRA55'):
			precipStr='SNowfallRate' 
		else:
			precipStr='SNowfall'
	elif (varStr=='tp'):
		if (reanalysis=='MERRA' or reanalysis=='MERRA_2'):
			precipStr='PrecipTotal'
		elif (reanalysis=='CFSR'):
			precipStr='PrecipRate'
		elif (reanalysis=='JRA55'):
			precipStr='TotalPrecip'
	
	print outPath
	print precipStr

	#varStr='sftv2'

	if not os.path.exists(outPath+varStr+'/'+str(yearT)):
		os.makedirs(outPath+varStr+'/'+str(yearT))


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
		dayStr='%03d' %(dayT)


		#in  kg/m2 per day
		# NOTE THE PLUS 1 ON THE DAY AS MERRA/MERRA2/JRA START AT 001 not 000
		Precip =pF.get_GRIDDED_precip_days(m, reanalysisDataPath, precipStr, reanalysis, yearT, dayT+1, mask_val=-999, multi_var=0, offset=0)

	
		PrecipG = griddata((xptsI.flatten(), yptsI.flatten()), Precip.flatten(), (xptsG, yptsG), method='linear')
		#PrecipG=PrecipG[0]
		PrecipG[where(region_maskG>10)]=0
		PrecipG=PrecipG.astype('f2')
		
		#pF.plotSnow(m, xptsG, yptsG, PrecipG, out=figpath+'/'+varStr+'-'+str(yearT)+'_d'+str(dayT)+'T2', units_lab=r'kg/m2', minval=0, maxval=10, base_mask=0, norm=0, cmap_1=cm.viridis)

		PrecipG.dump(outPath+varStr+'/'+str(yearT)+'/'+reanalysis+varStr+dxStr+'-'+str(yearT)+'_d'+dayStr)

#-- run main program
if __name__ == '__main__':
	for y in xrange(2000, 2015+1, 1):
		for reanalysis in ['MERRA2']: # NCEP_R2, NCEP_R1, NCEP_CIRES
			
			main(y, reanalysis)








#lonS, latS=np.meshgrid(*(lon, lat))


#varMYr = varM[numdays]

	

