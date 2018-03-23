""" NESOSIM.py
	
	The NASA Euelerian Snow on Sea Ice Model (NESOSIM). 
	Model written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov) or refer to the GitHub site (https://github.com/akpetty/NESOSIM)

	Run this python script with the run.py script in this same directory. 

	Input:
		Gridded/daily data of snowfall, ice drift, ice concentration, wind speeds

	Output:
		Gridded/daily data of the snow depth/density and snow budget terms.
		The DataOutput/MODELRUN/budgets/ netcdf files are all the snow budget terms needed for the analysis scripts/
		The DataOutput/MODELRUN/final/ netcdf files are the finalized netcdf files of the key variables, including metadata.

	Python dependencies:
		See below for the relevant module imports. Of note:
		xarray/pandas
		netCDF4
		matplotlib
		basemap

		More information on installation is given in the README file.

	Update history:
		03/01/2018: Version 1

    
"""

from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from pylab import *
import numpy.ma as ma
import xarray as xr
import pandas as pd
import os
from glob import glob
from scipy.ndimage.filters import gaussian_filter
import netCDF4 as nc4
#from matplotlib import rc
#from netCDF4 import Dataset
import commonFuncs as cF



def OutputSnowModelRaw(savePath, saveStr, snowDepths, density, \
	precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, snowDiv, snowWind, snowWindPack):
	""" Output snow model data using xarray

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved

	Output:
		xarray data as basic netCDF files   
    
    """

	precipData = xr.DataArray(precipDays, dims=('time', 'x', 'y'))
	snowDepthsData = xr.DataArray(snowDepths, dims=('time', 'lyrs',  'x', 'y'))
	snowAccData = xr.DataArray(snowAcc, dims=('time', 'x', 'y'))
	snowDivData = xr.DataArray(snowDiv, dims=('time', 'x', 'y'))
	snowAdvData = xr.DataArray(snowAdv, dims=('time', 'x', 'y'))
	snowWindData = xr.DataArray(snowWind, dims=('time', 'x', 'y'))
	snowWindPackData=xr.DataArray(snowWindPack, dims=('time', 'x', 'y'))
	snowOceanData = xr.DataArray(snowOcean, dims=('time', 'x', 'y'))
	#snowRidgeData = xr.DataArray(snowRidge, dims=('time', 'x', 'y'))
	#snowDcationData = xr.DataArray(snowDcation, dims=('time', 'x', 'y'))
	densityData = xr.DataArray(density, dims=('time', 'x', 'y'))
	iceConcData = xr.DataArray(iceConcDays, dims=('time', 'x', 'y'))
	windData = xr.DataArray(windDays, dims=('time', 'x', 'y'))

	dataSet = xr.Dataset({'Precip': precipData, 'snowDepth': snowDepthsData, 'snowAcc': snowAccData, 'snowDiv': \
		snowDivData,'snowAdv': snowAdvData, 'snowWind': snowWindData,'snowWindPack': snowWindPackData,'snowOcean': \
		snowOceanData, 'density': densityData, 'iceConc': iceConcData, 'winds': windData})


	dataSet.to_netcdf(savePath+'/budgets/'+saveStr) 

def OutputSnowModelFinal(savePath, saveStr, lons, lats, snowVolT,snowDepthT, densityT, iceConcT, precipT, datesT):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  

	Output:

    
    """

     
	f = nc4.Dataset(savePath+'/final/'+saveStr+'.nc','w', format='NETCDF4') #'w' stands for write
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

def calcLeadLoss(snowDepthT, WindDayT, iceConcDaysT):
	""" Snow loss to leads due to wind forcing

	Use a variable leadlossfactor parameter. This is relatively unconstrained!

	Args:
		snowDepthT (var): Daily gridded snowdepth 
		WindDayT (var): Daily gridded wind magnitude
		iceConcDaysT (var): Daily gridded ice concentration

	returns:
		sowWindPackLossT (var): Snow lost from fresh snow layer

	"""
	#
	 
	snowWindT = -(leadLossFactor*snowDepthT*WindDayT*(1-iceConcDaysT)) #*iceConcDaysG[x]
	return snowWindT

def calcWindPacking(windDayT, snowDepthT0):
	""" Snow pack densification through wind packing

	Calculated using the amount of snow packed down and the 
	difference in density between the fresh snow density and the old snow density

	Args:
		snowDepthT (var): Daily gridded snowdepth 
		WindDayT (var): Daily gridded wind magnitude
		iceConcDaysT (var): Daily gridded ice concentration

	returns:
		snowWindPackLossT (var): Snow lost from fresh snow layer
		snowWindPackGainT (var): Snow gained to old snow layer
		snowWindPackNetT (var): Net snow gain/loss

	"""


	windT= where(windDayT>windPackThresh, 1, 0)
	
	# snow loss from fresh layer through wind packing to old layer
	snowWindPackLossT=-windPackFactor*windT*snowDepthT0 #*iceConcDaysG[x]

	# snow gain to old layer through wind packing from fresh layer
	snowWindPackGainT=windPackFactor*windT*snowDepthT0*(snowDensityFresh/snowDensityOld) #*iceConcDaysG[x]
	
	snowWindPackNetT=snowWindPackLossT+snowWindPackGainT#*iceConcDaysG[x]
	return snowWindPackLossT, snowWindPackGainT, snowWindPackNetT

def calcRidgeLoss(driftGday, snowDepthT, iceConcDaysT):
	""" Snow loss through ridging

	NB NOT USED IN CURRENT BUDGET CALCULATIONS

	Args:
		snowDepthT (var): Daily gridded snowdepth 
		driftGday (var): Daily gridded ice drift
		iceConcDaysT (var): Daily gridded ice concentration

	returns:
		snowRidgeT (var): Snow lost through ridging

	""" 

	# Snow loss from ridging as a percentage of current snow depth
	snowRidgeFactor=0.5
	coefb=0.05
	maxConc=0.95

	# Calculate convergence/divergence without snow
	dvelxdx = np.gradient(driftGday[0], dx, axis=(1)) #convert from m/s to m per day, #1 here is the columns, so in the x direction
	dvelydy = np.gradient(driftGday[1], dx, axis=(0)) #0 here is the rows, so in the y direction
    
	dvelxdx[np.isnan(dvelxdx.filled(0.))]=0.
	dvelydy[np.isnan(dvelydy.filled(0.))]=0.

	DynT= -(dvelxdx + dvelydy)
	
	# See where we have convergence and high ice concentration.
	# Could change to just losing some of the new snow fraction
	convergeT= where((DynT>0)&(iceConcDaysT>=maxConc), 1, 0)
	snowRidgeT= -convergeT*snowRidgeFactor*snowDepthT*((1-iceConcDaysT)/coefb) #*iceConcDaysG[x]
	return snowRidgeT

def calcDynamics(driftGday, snowDepthsT, dx):
	""" Snow loss/gain from ice dynamics

	Args:
		driftGday (var): Daily gridded ice drift
		snowDepthT (var): Daily gridded snowdepth 
		dx (var): grid spacing

	returns:
		snowAdvAllT (var): Snow change through advection (gain is positive)
		snowDivAllT (var): Snow change through convergence/divergence (convergence is positive)
		
	"""

	dhsvelxdxDiv = snowDepthsT*np.gradient(driftGday[0]*60.*60.*24., dx, axis=(1)) #convert from m/s to m per day, #1 here is the columns, so in the x direction
	dhsvelydyDiv = snowDepthsT*np.gradient(driftGday[1]*60.*60.*24., dx, axis=(0)) #0 here is the rows, so in the y direction

	# fill masked and nans with 0  
	dhsvelxdxDiv[np.isnan(dhsvelxdxDiv)]=0.
	dhsvelxdxDiv=ma.filled(dhsvelxdxDiv, 0.)
	dhsvelxdxDiv[~np.isfinite(dhsvelxdxDiv)]=0.
	dhsvelydyDiv[np.isnan(dhsvelydyDiv)]=0.
	dhsvelydyDiv=ma.filled(dhsvelydyDiv, 0.)
	dhsvelydyDiv[~np.isfinite(dhsvelydyDiv)]=0.

	#print 'dh:', np.amax(dhsvelydyDiv)

	snowDivAllT= -(dhsvelxdxDiv + dhsvelydyDiv)
	snowDivAllT=ma.filled(snowDivAllT, 0.)

	# Snow change from ice dynamics (divergence/convergence). Convergence is positive
	dhsvelxdxAdv = driftGday[0]*60.*60.*24.*np.gradient(snowDepthsT, dx, axis=(2)) #convert from m/s to m per day, #1 here is the columns, so in the x direction
	dhsvelydyAdv = driftGday[1]*60.*60.*24.*np.gradient(snowDepthsT, dx, axis=(1))  #0 here is the rows, so in the y direction

	#print 'dh1:', np.amax(dhsvelxdxAdv)
	# fill masked and nans with 0  
	dhsvelxdxAdv[np.isnan(dhsvelxdxAdv)]=0.
	dhsvelxdxAdv[~np.isfinite(dhsvelxdxAdv)]=0.
	dhsvelxdxAdv=ma.filled(dhsvelxdxAdv, 0.)
	dhsvelydyAdv[np.isnan(dhsvelydyAdv)]=0.
	dhsvelydyAdv=ma.filled(dhsvelydyAdv, 0.)
	dhsvelydyAdv[~np.isfinite(dhsvelydyAdv)]=0.
	#print 'dh2:', np.amax(dhsvelxdxAdv)


	snowAdvAllT= -(dhsvelxdxAdv + dhsvelydyAdv)
	snowAdvAllT=ma.filled(snowAdvAllT, 0.)


	# Set limits on how much snow can be lost? 
	# May be redundant
	mask0=where(-snowAdvAllT[0]>snowDepthsT[0])
	snowAdvAllT[0][mask0]=-snowDepthsT[0][mask0]

	mask1=where(-snowAdvAllT[1]>snowDepthsT[1])
	snowAdvAllT[1][mask1]=-snowDepthsT[1][mask1]

	mask2=where(-snowDivAllT[0]>snowDepthsT[0])
	snowDivAllT[0][mask2]=-snowDepthsT[0][mask2]

	mask3=where(-snowDivAllT[1]>snowDepthsT[1])
	snowDivAllT[1][mask3]=-snowDepthsT[1][mask3]
	

	return snowAdvAllT, snowDivAllT

def calcBudget(snowDepths, iceConcDayT, precipDayT, driftGdayT, windDayT, 
	density, precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, 
	snowDiv, snowWind, snowWindPackLoss, snowWindPackGain, snowWindPack, region_maskG, dx, x, dayT,
	densityType='variable', dynamicsInc=1, leadlossInc=1, windpackInc=1):
	""" Main snow budget calculation function

	"""

	# Convert current fields to new temp arrays
	precipDays[x]=precipDayT
	iceConcDays[x]=iceConcDayT
	windDays[x]=windDayT

	print 'Density:', densityType
	
	if (densityType=='clim'):
		# returns a fixed density value assigned to all grid cells based on climatology. 
		# Applies the same value to both snow layers.
		snowDensityNew=densityClim(dayT)
	else:
		# Two layers so a new snow density and an evolving old snow density
		snowDensityNew=snowDensityFresh

		# Convert precip to m
	precipDayDelta=precipDayT/snowDensityNew
	precipDayDelta[where(region_maskG>10)]=np.nan

	# Snow accumulated onto the ice
	snowAccDelta= (precipDayDelta * iceConcDayT)
	snowAcc[x+1] = snowAcc[x] + snowAccDelta

	# Snow deposited into the ocean
	snowOceanDelta= -(precipDayDelta * (1-iceConcDayT))
	snowOcean[x+1] = snowOcean[x] + snowOceanDelta

	# Snow change from dynamics
	if (dynamicsInc==1):
		snowAdvDelta, snowDivDelta = calcDynamics(driftGdayT, snowDepths[x], dx)
	else:
		snowAdvDelta=np.zeros((iceConcDayT.shape))
		snowDivDelta=np.zeros((iceConcDayT.shape))
	
	snowAdv[x+1] = snowAdv[x] + snowAdvDelta[0]+ snowAdvDelta[1]
	snowDiv[x+1] = snowDiv[x] + snowDivDelta[0]+ snowDivDelta[1]
	
	# Snow loss from ridging
	#snowRidgeT= calcRidgeLoss(driftGdays[x], snowDepthsT[0], iceConcDaysG[x])	
	#snowRidge[x+1] = snowRidge[x]+ snowRidgeT

	# Lead loss term
	if (leadlossInc==1):
		snowWindDelta= calcLeadLoss(snowDepths[x, 0], windDayT, iceConcDayT)
	else:
		snowWindDelta=np.zeros((iceConcDayT.shape))

	snowWind[x+1]=snowWind[x] + snowWindDelta

	# Wind packing
	if (windpackInc==1):
		snowWindPackLossDelta, snowWindPackGainDelta, snowWindPackNetDelta=calcWindPacking(windDayT, snowDepths[x, 0])
	else:
		snowWindPackLossDelta =np.zeros((iceConcDayT.shape))
		snowWindPackGainDelta=np.zeros((iceConcDayT.shape))
		snowWindPackNetDelta=np.zeros((iceConcDayT.shape))

	snowWindPackLoss[x+1]=snowWindPackLoss[x]+snowWindPackLossDelta
	snowWindPackGain[x+1]=snowWindPackGain[x]+snowWindPackGainDelta
	snowWindPack[x+1]=snowWindPack[x]+snowWindPackNetDelta



	snowDepths[x+1, 0]=snowDepths[x, 0]+snowAccDelta  +snowWindPackLossDelta +snowWindDelta +snowAdvDelta[0]+snowDivDelta[0] #+snowRidgeT
	# Old snow layer
	snowDepths[x+1, 1]=snowDepths[x, 1] +snowWindPackGainDelta + snowAdvDelta[1] + snowDivDelta[1] #+ snowDcationT

	# Set negative (false) snow values to zero
	snowDepths[x+1, 0][where(snowDepths[x+1, 0]<0.)]=0.
	snowDepths[x+1, 1][where(snowDepths[x+1, 1]<0.)]=0.
	snowDepths[x+1][where(np.isnan(snowDepths[x+1]))]=0.
	
	print 'sd1:', np.amax(snowDepths[x+1, 0])
	print 'sd11:', np.amax(snowDepths[x+1, 1])
	#snowDepths[x+1].filled(0.)

	snowDepths[x+1, 0] = gaussian_filter(snowDepths[x+1, 0], sigma=0.3)
	snowDepths[x+1, 1] = gaussian_filter(snowDepths[x+1, 1], sigma=0.3)

	print 'sd2:',np.amax(snowDepths[x+1])

	# Do this again after smoothing
	snowDepths[x+1, 0][where(snowDepths[x+1, 0]<0.)]=0.
	snowDepths[x+1, 1][where(snowDepths[x+1, 1]<0.)]=0.
	snowDepths[x+1][where(np.isnan(snowDepths[x+1]))]=0.
	
	print 'sd3:',np.amax(snowDepths[x+1])

	#snowDepths[x+1].filled(0.)
	# Set snow depths over land/coasts to zero
	snowDepths[x+1, 0][where(region_maskG>10)]=0.
	snowDepths[x+1, 1][where(region_maskG>10)]=0.

	if (densityType=='clim'):
		# returns a fixed density value assigned to all grid cells based on climatology. 
		# Applies the same value to both snow layers.
		
		density[x+1]=snowDensityNew
		density[x+1][where(region_maskG>10)]=np.nan
		density[x+1][where(iceConcDayT<minConc)]=np.nan
		density[x+1][where((snowDepths[x+1][0]+snowDepths[x+1][1])<minSnowD)]=np.nan
	else:
		# Two layers so a new snow density and an evolving old snow density
		
		density[x+1]=densityCalc(snowDepths[x+1], iceConcDayT, region_maskG)


def genEmptyArrays(numDaysT, nxT, nyT):
	""" Declare empty arrays for various budget terms"""
	
	precipDays=np.zeros((numDaysT, nxT, nyT)) 
	iceConcDays=np.zeros((numDaysT, nxT, nyT)) 
	windDays=np.zeros((numDaysT, nxT, nyT)) 

	snowDepths=np.zeros((numDaysT, 2, nxT, nyT))
	density=np.zeros((numDaysT, nxT, nyT))

	snowDiv=np.zeros((numDaysT, nxT, nyT))
	snowAdv=np.zeros((numDaysT, nxT, nyT))
	snowAcc=np.zeros((numDaysT, nxT, nyT))
	snowOcean=np.zeros((numDaysT, nxT, nyT))
	snowWindPack=np.zeros((numDaysT, nxT, nyT))
	snowWindPackLoss=np.zeros((numDaysT, nxT, nyT))
	snowWindPackGain=np.zeros((numDaysT, nxT, nyT))
	snowWind=np.zeros((numDaysT, nxT, nyT))
	#snowRidge=np.zeros((numDays, nx, ny))
	#snowDcation=np.zeros((numDays, nx, ny))

	#densDyn=np.zeros((numDays, nx, ny))
	#densDcation=np.zeros((numDays, nx, ny))
	
	return precipDays, iceConcDays, windDays, snowDepths, density, snowDiv, snowAdv, snowAcc, snowOcean, snowWindPack, \
	snowWindPackLoss, snowWindPackGain, snowWind

def plotEndBudgets(m, xptsG, yptsG, precipDaysT, windDaysT, snowDepthsT, snowOceanT, snowAccT, snowDivT, \
	snowAdvT, snowWindT, snowWindPackT, snowWindPackLossT, snowWindPackGainT, densityT, dateStr, totalOutStr='test'):
	""" Plot snow budget terms """

	cF.plotSnow(m, xptsG, yptsG, precipDaysT, date_string=dateStr, out=figpath+'/precip'+totalOutStr, units_lab=r'kg/m2', minval=0., maxval=1., base_mask=1, norm=0, cmap_1=cm.cubehelix_r)
	cF.plotSnow(m, xptsG, yptsG, windDaysT, date_string=dateStr, out=figpath+'/wind'+totalOutStr, units_lab=r'm/s', minval=0., maxval=10, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)

	cF.plotSnow(m, xptsG, yptsG, snowDepthsT[0], date_string=dateStr, out=figpath+'/snowNew_'+totalOutStr, units_lab=r'm', minval=0., maxval=0.5, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)
	cF.plotSnow(m, xptsG, yptsG, snowDepthsT[1], date_string=dateStr, out=figpath+'/snowOld_'+totalOutStr, units_lab=r'm', minval=0., maxval=0.5, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)
	cF.plotSnow(m, xptsG, yptsG, (snowDepthsT[0]+snowDepthsT[1]), date_string=dateStr, out=figpath+'/snowTot_'+totalOutStr, units_lab=r'm', minval=0., maxval=0.5, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)

	cF.plotSnow(m, xptsG, yptsG, snowOceanT, date_string=dateStr, out=figpath+'/snowOcean_'+totalOutStr, units_lab=r'm', minval=-0.6, maxval=0, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)
	cF.plotSnow(m, xptsG, yptsG, snowAccT, date_string=dateStr, out=figpath+'/snowAcc_'+totalOutStr, units_lab=r'm', minval=0., maxval=0.6, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)

	cF.plotSnow(m, xptsG, yptsG, snowDivT, date_string=dateStr, out=figpath+'/snowDiv_'+totalOutStr, units_lab=r'm', minval=-0.3, maxval=0.3, base_mask=1, norm=0, cmap_1=cm.RdBu)
	cF.plotSnow(m, xptsG, yptsG, snowAdvT, date_string=dateStr, out=figpath+'/snowAdv_'+totalOutStr, units_lab=r'm', minval=-0.3, maxval=0.3, base_mask=1, norm=0, cmap_1=cm.RdBu)
	cF.plotSnow(m, xptsG, yptsG, snowWindT, date_string=dateStr, out=figpath+'/snowWind_'+totalOutStr, units_lab=r'm', minval=-0.3, maxval=0.3, base_mask=1, norm=0, cmap_1=cm.RdBu)
	cF.plotSnow(m, xptsG, yptsG, snowWindPackT, date_string=dateStr, out=figpath+'/snowWindPackNet_'+totalOutStr, units_lab=r'm', minval=-0.3, maxval=0., base_mask=1, norm=0, cmap_1=cm.cubehelix)
	cF.plotSnow(m, xptsG, yptsG, snowWindPackLossT, date_string=dateStr, out=figpath+'/snowWindPackLoss_'+totalOutStr, units_lab=r'm', minval=-0.3, maxval=0., base_mask=1, norm=0, cmap_1=cm.cubehelix)
	cF.plotSnow(m, xptsG, yptsG, snowWindPackGainT, date_string=dateStr, out=figpath+'/snowWindPackGain_'+totalOutStr, units_lab=r'm', minval=0, maxval=0.3, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)

	cF.plotSnow(m, xptsG, yptsG, densityT, date_string=dateStr, out=figpath+'/dens'+totalOutStr, units_lab=r'kg/m3', minval=300, maxval=360, base_mask=1, norm=0, cmap_1=cm.cubehelix_r)



def loadData(yearT, dayT, driftP, reanalysisP, reanalysisWind,varStr, windStr, dxStr, team_s, yearT2=0):
	"""Load daily model forcings"""

	dayStr='%03d' %dayT

	if (yearT2<1900):
		# If we dont feed in a different forcing year then just set this to the main forcing year
		yearT2=yearT


	print yearT, yearT2
	precipDayG=load(dataPath+'Reanalysis/'+reanalysisP+'/'+varStr+'/'+str(yearT)+'/'+reanalysisP+varStr+dxStr+'-'+str(yearT)+'_d'+dayStr)
	windDayG=load(dataPath+'Reanalysis/'+reanalysisWind+'/'+windStr+'/'+str(yearT)+'/'+reanalysisWind+windStr+dxStr+str(yearT)+'_d'+dayStr)
	
	iceConcDayG=load(dataPath+'IceConc/'+team_s+'/'+str(yearT2)+'/iceConcG_'+team_s+dxStr+'-'+str(yearT2)+'_d'+dayStr)
	
	#print 'Path:', glob(outPath+'/Drifts/'+driftP+'/'+str(yearT)+'/'+driftP+'DriftG'+dxStr+'-'+str(yearT)+'_d'+dayStr)
	#print 'Num drift files:', size(glob(outPath+'/Drifts/'+driftP+'/'+str(yearT)+'/'+driftP+'DriftG'+dxStr+'-'+str(yearT)+'_d'+dayStr))
	if (size(glob(dataPath+'Drifts/'+driftP+'/'+str(yearT2)+'/'+driftP+'DriftG'+dxStr+'-'+str(yearT2)+'_d'+dayStr))>0):
		driftGdayG=load(dataPath+'Drifts/'+driftP+'/'+str(yearT2)+'/'+driftP+'DriftG'+dxStr+'-'+str(yearT2)+'_d'+dayStr)
	else:
		# if no drifts exist for that day then just set drifts to masked array (i.e. no drift).
		driftGdayG=ma.masked_all((2, iceConcDayG.shape[0], iceConcDayG.shape[1]))
	

	return iceConcDayG, precipDayG, driftGdayG, windDayG

def densityCalc(snowDepthsT, iceConcDayT, region_maskT):
	"""Assign initial snow densities based on snow depths"""

	# Assign initial density
	densityT=((snowDepthsT[0]*snowDensityFresh) + (snowDepthsT[1]*snowDensityOld))/(snowDepthsT[0]+snowDepthsT[1]) #+ densDcationT
	
	densityT[where(densityT>snowDensityOld)]=snowDensityOld
	densityT[where(densityT<snowDensityFresh)]=snowDensityFresh

	densityT[where(region_maskT>10)]=np.nan
	densityT[where(iceConcDayT<minConc)]=np.nan
	densityT[where((snowDepthsT[0]+snowDepthsT[1])<minSnowD)]=np.nan

	return densityT

def densityClim(dayT):
	"""Assign initial snow density based on daily climatology"""

	densityClim=pd.read_csv(dataPath+'/Daily_Density.csv', header=0, names=['Day', 'Density'])
	#find density of given day and multiply by 1000 to express in Kg
	densityClimDay=1000*densityClim['Density'].iloc[dayT-1]

	return densityClimDay


def main(year1, month1, day1, yearIC1=0, reanalysisP='ERAI', varStr='sf', driftP='NSIDCv3', 
	team_s='nt', densityTypeT='variable', outStr='', IC=0, windPackFactorT=0.1, windPackThreshT=5., 
	leadLossFactorT=0.1, dynamicsInc=1, leadlossInc=1, windpackInc=1, saveData=1, plotBudgets=1, saveFolder=''):
	""" Primary model function 

	Add in more of a description here

	Args:
		various

	"""

	# Map projection
	m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False  )
	
	global dataPath
	global outPath
	
	outPath='../DataOutput/'

	# Or try ../TestData/
	dataPath='../Data/'

	# Grid info
	dx=100000.
	dxStr=str(int(dx/1000))+'km'
	lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)
	region_maskG=load(dataPath+'/Extra/regionMaskG'+dxStr)

	# Products used
	extraStr=''
	reanalysisWind='ERAI' # NCEPR2, 
	windStr='WindMag'
	#reanalysisP='ERAI'
	#varStr='sf' # tp is precip, sf is snowfall for the reanalyses that provide that, e.g. ERA-I
	

	# Assign density of the two snow layers
	global snowDensityFresh, snowDensityOld, minSnowD, minConc, leadLossFactor, windPackThresh, windPackFactor
	
	snowDensityFresh=200. # density of fresh snow layer
	snowDensityOld=350. # density of old snow layer

	minSnowD=0.02 # minimum snow depth for density estimate
	minConc=0.15 # mask budget values with a conc below this

	leadLossFactor=leadLossFactorT #0.025 # snow loss to leads coefficient
	windPackThresh=windPackThreshT # 5. #minimum winds needed for wind packing
	windPackFactor=windPackFactorT # 0.05 #fraction of snow packed into old snow layer
	
	# Date info
	#year1=2013
	#month1=7 # 8=September
	#day1=14

	year2=year1+1
	month2=4 # 4=May
	day2=0


	yearIC2=yearIC1+1

	year=year1
	yearIC=yearIC1

	# Get time period info
	startDay, numDays, numDaysYear1, dateOut= cF.getDays(year1, month1, day1, year2, month2, day2)
	print startDay, numDays, numDaysYear1, dateOut

	# make this into a small function
	import datetime
	dates=[]
	for x in xrange(1, numDays+1):
		print x
		date = datetime.datetime(year1, month1+1, day1+1) + datetime.timedelta(x) #This assumes that the year is 2007
		print int(date.strftime('%Y%m%d'))
		dates.append(int(date.strftime('%Y%m%d')))


	saveStr= driftP+'_'+extraStr+reanalysisP+'_'+varStr+'_SIC'+team_s+'_Rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+outStr+'-'+dateOut
	saveStrNoDate=driftP+'_'+extraStr+reanalysisP+'_'+varStr+'_SIC'+team_s+'_Rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+outStr
	print saveStr
	 #'../../DataOutput/'

	savePath=outPath+saveFolder+'/'+saveStrNoDate
	# Declare empty arrays for compiling budgets
	if not os.path.exists(savePath+'/budgets/'):
		os.makedirs(savePath+'/budgets/')
	if not os.path.exists(savePath+'/final/'):
		os.makedirs(savePath+'/final/')

	global figpath
	figpath='../Figures/Diagnostic/'+saveStrNoDate+'/'
	if not os.path.exists(figpath):
		os.makedirs(figpath)

	precipDays, iceConcDays, windDays, snowDepths, density, snowDiv, snowAdv, snowAcc, snowOcean, snowWindPack, snowWindPackLoss, snowWindPackGain, snowWind= genEmptyArrays(numDays, nx, ny)

	if (IC>0):
		if (IC==1):
			#ICSnowDepth = load(outPath+'AugSnow'+dxStr)
			if (yearIC1>1900):
				# if we read in a valid ice conc year
				ICSnowDepth = load(dataPath+'InitialConds/ICsnow'+str(yearIC1)+'-'+dxStr)
			else:
				ICSnowDepth = load(dataPath+'InitialConds/ICsnow'+str(year1)+'-'+dxStr)
		elif (IC==2):
			ICSnowDepth = load(dataPath+'InitialConds/AugSnow'+dxStr)

		iceConcDayG, precipDayG, driftGdayG, windDayG=loadData(year1, startDay, driftP, reanalysisP, reanalysisWind,varStr, windStr, dxStr, team_s, yearT2=yearIC1)
		ICSnowDepth[where(iceConcDayG<minConc)]=0

		#pF.plotSnow(m, xptsG, yptsG, ICSnowDepth, date_string='T', out=figpath+'/Snow/2layer/snowIC', units_lab=r'm', minval=-0, maxval=.1, base_mask=0, norm=0, cmap_1=cm.viridis)

		# Split the initial snow depth over both layers
		snowDepths[0, 0]=ICSnowDepth*0.5
		snowDepths[0, 1]=ICSnowDepth*0.5


	#pF.plotSnow(m, xptsG, yptsG, densityT, date_string=str(startDay-1), out=figpath+'/Snow/2layer/densityD'+driftP+extraStr+reanalysisP+varStr+'_sy'+str(year1)+'d'+str(startDay)+outStr+'T0', units_lab=r'kg/m3', minval=180, maxval=360, base_mask=0, norm=0, cmap_1=cm.viridis)

	# Loop over days 
	for x in xrange(numDays-1):	
		day = x+startDay
		print 'day:', day
		if (day>=numDaysYear1):
			# If day goes beyond the number of days in initial year, jump to the next year
			day=day-numDaysYear1
			year=year2
			yearIC=yearIC2

		# Load daily data 
		iceConcDayG, precipDayG, driftGdayG, windDayG=loadData(year, day, driftP, reanalysisP, reanalysisWind,varStr, windStr, dxStr, team_s, yearT2=yearIC)
		
		# Calculate snow budgets
		calcBudget(snowDepths, iceConcDayG, precipDayG, driftGdayG, windDayG, 
			density, precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, 
			snowDiv, snowWind, snowWindPackLoss, snowWindPackGain, snowWindPack, region_maskG, dx, x, day,
			densityType=densityTypeT, dynamicsInc=dynamicsInc, leadlossInc=leadlossInc, windpackInc=windpackInc)
		
	# Load last data 
	iceConcDayG, precipDayG, _, windDayG=loadData(year, day+1, driftP, reanalysisP, reanalysisWind,varStr, windStr, dxStr, team_s, yearT2=yearIC)
	precipDays[x+1]=precipDayG
	iceConcDays[x+1]=iceConcDayG
	windDays[x+1]=windDayG
	
	outStrings=['snowDepthTotal','snowDepthTotalConc', 'density', 'iceConc', 'Precip']

	if (saveData==1):
		# Output snow budget terms to netcdf datafiles
		OutputSnowModelRaw(savePath, saveStr, snowDepths, density, precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, snowDiv, snowWind, snowWindPack)
		OutputSnowModelFinal(savePath, saveStr, lonG, latG, snowDepths[:, 0]+snowDepths[:, 1], (snowDepths[:, 0]+snowDepths[:, 1])/iceConcDays, density, iceConcDays, precipDays, dates)

	if (plotBudgets==1):
		# Plot final snow budget terms 
		plotEndBudgets(m, xptsG, yptsG, precipDayG, windDayG, snowDepths[x+1], snowOcean[x+1], snowAcc[x+1], snowDiv[x+1], \
		snowAdv[x+1], snowWind[x+1], snowWindPack[x+1], snowWindPackLoss[x+1], snowWindPackGain[x+1], density[x+1], dates[-1], totalOutStr=saveStr)

