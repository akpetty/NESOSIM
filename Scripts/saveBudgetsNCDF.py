############################################################## 
# Date: 01/01/16
# Name: saveBudgetsNCDF.py
# Author: Alek Petty
# Description: Script to psave budgets as netCDF
# Input requirements: budget data

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
import netCDF4 as nc4
import os

def outData(lons, lats, snowVolT,snowDepthT, densityT, iceConcT, precipT, datesT, totalOutStr, folderPath):

	f = nc4.Dataset(folderPath+'/snowBudget_'+totalOutStr+'.nc','w', format='NETCDF4') #'w' stands for write
	#tempgrp = f.createGroup('DragData')
	print 'dimensions:', lons.shape[0], lons.shape[1], snowVolT.shape[0]
	f.createDimension('x', lons.shape[0])
	f.createDimension('y', lons.shape[1])
	f.createDimension('day', snowVolT.shape[0])

	longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
	latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  
	snowDepth = f.createVariable('snowDepth', 'f4', ('day', 'x', 'y'))
	snowVol = f.createVariable('snowVol', 'f4', ('day', 'x', 'y'))
	density = f.createVariable('density', 'f4', ('day', 'x', 'y'))
	precip = f.createVariable('Precip', 'f4', ('day', 'x', 'y'))
	iceConc = f.createVariable('iceConc', 'f4', ('day', 'x', 'y'))
	day = f.createVariable('day', 'i4', ('day'))

	#dates = f.createDimension('dates', None)

	#startYr = f.createVariable('startYear', 'i2')
	#date_range = f.createVariable('year', 'str')

	longitude.units = 'degrees East'
	latitude.units = 'degrees North'
	#Cdafr.units = ''
	#Cda.units = ''
	snowVol.description = 'Daily snow volume per unit grid cell (m)'
	snowDepth.description = 'Daily snow depth (m)'
	density.description = 'Bulk snow density (kg/m3)'
	precip.description = 'Precipitation, normally the explicit snowfall component of precip given in the filename (kg/m2)'
	iceConc.description = 'Ice concentration, product given in the filename (nt = NASA Team, bt = Bootstrap)'
	#dates.description = 'Date'
	
	longitude[:] = np.around(lons, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats, decimals=4)
	snowVol[:] = np.around(snowVolT, decimals=4)
	snowDepth[:] = np.around(snowDepthT, decimals=4)
	density[:] = np.around(densityT, decimals=4)
	iceConc[:] = np.around(iceConcT, decimals=4)
	precip[:] = np.around(precipT, decimals=4)
	day[:]=datesT
	#date_range[:]=date1+'-'+date2
	#get time in days since Jan 01,01
	from datetime import datetime
	today = datetime.today()
	#time_num = today.toordinal()
	#time[0] = time_num

	#Add global attributes
	f.author = "Alek Petty"
	f.contact = " alek.a.petty@nasa.gov, www.alekpetty.com, @alekpetty"
	f.description = "Daily NESOSI snow budget data"
	f.history = "Created " + today.strftime("%d/%m/%y")
	f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()


m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False  )
#m = Basemap(projection='npstere',boundinglat=30.52,lon_0=0, resolution='l'  )


outPath='..//DataOutput/'
dataPath = '/Volumes/PETTY_PASSPORT2/pasta/Data/'

def main(year1, month1, day1, yearIC1=0, reanalysisP='ERAI', varStr='sf', driftP='NSIDCv3', 
	team_s='nt', densityTypeT='variable', outStr='', IC=0, windPackFactorT=0.1, windPackThreshT=5., 
	leadLossFactorT=0.1, dynamicsInc=1, leadlossInc=1, windpackInc=1, saveData=1, plotBudgets=1):

	# GRID INFO
	dx=100000.
	dxStr=str(int(dx/1000))+'km'
	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)


	year2=year1+1
	month2=4 # 4=May
	day2=0
	extraStr=''
	layers='2layer'

	
	outStrings=['snowDepthTotal','snowDepthTotalConc', 'density', 'iceConc', 'Precip']

	startDay, numDays, numDaysYear1, dateOut=pF.getDays(year1, month1, day1, year2, month2, day2)
	print numDays
	daysT=np.arange(startDay, startDay+numDays, 1)

	saveStr= driftP+'_'+extraStr+reanalysisP+'_'+varStr+'_SIC'+team_s+'_Rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+outStr+'-'+dateOut
	saveStrNoDate=driftP+'_'+extraStr+reanalysisP+'_'+varStr+'_SIC'+team_s+'_Rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+outStr
	

	print saveStr

	folderPath=outPath+'SnowModel/'+layers+'/'+saveStrNoDate+'/'
	print 'Folderpath:', folderPath
	snowBudgetT=pF.get_budgets_time_multi(outStrings, folderPath, 'xarray'+saveStr, startDay=0, endDay=numDays, mafill=1)

	#print snowBudgetT[0].shape
	#print snowBudgetT[1].shape

	import datetime
	dates=[]
	for x in xrange(1, numDays+1):
		print x
		date = datetime.datetime(year1, month1+1, day1+1) + datetime.timedelta(x) #This assumes that the year is 2007
		print int(date.strftime('%Y%m%d'))
		dates.append(int(date.strftime('%Y%m%d')))

	saveFolderPath=outPath+'netCDF/'+saveStrNoDate+'/'
	if not os.path.exists(saveFolderPath):
		os.makedirs(saveFolderPath)
	#print snowBudgetT[0].shape, size(dates)
	outData(lonG, latG, snowBudgetT[0],snowBudgetT[1], snowBudgetT[2], snowBudgetT[3], snowBudgetT[4], dates, saveStr, saveFolderPath)


	return snowBudgetT


#years = np.concatenate([np.arange(1980, 1990+1), np.arange(2004, 2014+1)])

years=np.arange(2000, 2004)
if __name__ == '__main__':
	for y in years:
		if (y==1987):
			continue
		print y
		#'ERAI', MERRA', 'MERRA_2', 'JRA55', 'NCEPR2', 'NCEPR1''NCEP_R1', 'NCEP_R2','NCEP_CIRES'
		reanalyses=['MEDIAN'] 
 		for reanalysis in reanalyses:
			main(y, 7, 14, reanalysisP=reanalysis, varStr='sf', driftP='NSIDCv3', team_s='bt', densityTypeT='variable', 
				outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
				dynamicsInc=1, leadlossInc=1, windpackInc=1)


# if __name__ == '__main__':
# 	for y in xrange(2000, 2014+1, 1):
# 		if (y==1987):
# 			continue
# 		print y
# 		reanalyses=['MERRA2']
#  		for reanalysis in reanalyses:
# 			snowT=main(y, 7, 14, reanalysisP=reanalysis, varStr='sf', driftP='NSIDCv3', team_s='bt', densityTypeT='variable', 
# 					outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
# 					dynamicsInc=1, leadlossInc=1, windpackInc=1)
		

# f = Dataset('../../Data/snowBudget_'+reanalysisP+driftP+'_'+str(dates[0])+'-'+str(dates[-1])+'.nc','r') #'w' stands for write
	
# datesT = f.variables['day'][:]
# latT = f.variables['latitude'][:]
# snowDepthT = f.variables['snowDepth'][:]
