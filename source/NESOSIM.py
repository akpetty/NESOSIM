""" NESOSIM.py
	
	The NASA Euelerian Snow on Sea Ice Model (NESOSIM) v1.1 
	Model written by Alek Petty
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
		cartopy

		More information on installation is given in the README file.

	Update history:
		1st March 2018: Version 0.1
		1st October 2018: Version 1.0 (updated through review process)
		1st May 2020: Version 1.1 (updated for ICESat-2 processing, new domain using cartopy/pyproj, bug fixes)
		1st October 2020: Version 1.1 (bug fixes, replace masked array with nan throughout, smoothed dynamics terms)
		15th October 2020: Version 1.1 (new atmosphere-wind loss term)

"""

import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import xarray as xr
import pandas as pd
import os
from glob import glob
from scipy.ndimage.filters import gaussian_filter
import netCDF4 as nc4
import utils as cF
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import datetime
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

def OutputSnowModelRaw(savePath, saveStr, snowDepths, density, \
	precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, snowDiv, snowLead, snowAtm, snowWindPack):
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
	snowLeadData = xr.DataArray(snowLead, dims=('time', 'x', 'y'))
	snowAtmData = xr.DataArray(snowAtm, dims=('time', 'x', 'y'))
	snowWindPackData=xr.DataArray(snowWindPack, dims=('time', 'x', 'y'))
	snowOceanData = xr.DataArray(snowOcean, dims=('time', 'x', 'y'))
	#snowRidgeData = xr.DataArray(snowRidge, dims=('time', 'x', 'y'))
	#snowDcationData = xr.DataArray(snowDcation, dims=('time', 'x', 'y'))
	densityData = xr.DataArray(density, dims=('time', 'x', 'y'))
	iceConcData = xr.DataArray(iceConcDays, dims=('time', 'x', 'y'))
	windData = xr.DataArray(windDays, dims=('time', 'x', 'y'))

	dataSet = xr.Dataset({'Precip': precipData, 'snowDepth': snowDepthsData, 'snowAcc': snowAccData, 'snowDiv': \
		snowDivData,'snowAdv': snowAdvData, 'snowLead': snowLeadData,'snowAtm': snowAtmData,'snowWindPack': snowWindPackData,'snowOcean': \
		snowOceanData, 'density': densityData, 'iceConc': iceConcData, 'winds': windData})

	print ('saving to:', savePath+'/budgets/'+saveStr)

	dataSet.to_netcdf(savePath+'/budgets/'+saveStr+'.nc') 


def OutputSnowModelFinal(savePath, saveStr, lons, lats, snowVolT,snowDepthT, densityT, iceConcT, precipT, windsT, tempT, datesT):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  
    
    """

     
	f = nc4.Dataset(savePath+'/final/'+saveStr+'.nc','w', format='NETCDF4') #'w' stands for write
	#tempgrp = f.createGroup('DragData')
	print ('dimensions:', lons.shape[0], lons.shape[1], snowVolT.shape[0])
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
	winds = f.createVariable('windMag', 'f4', ('day', 'x', 'y'))
	temps = f.createVariable('airTemp2m', 'f4', ('day', 'x', 'y'))
	day = f.createVariable('day', 'i4', ('day'))

	longitude.units = 'degrees East'
	latitude.units = 'degrees North'

	snowVol.description = 'Daily snow volume per unit grid cell (m)'
	snowDepth.description = 'Daily snow depth (effetive over the ice fraction) (m)'
	density.description = 'Bulk snow density (kg/m3)'
	precip.description = 'Precipitation, normally the explicit snowfall component of precip given in the filename (kg/m2)'
	iceConc.description = 'Ice concentration, product given in the filename (nt = NASA Team, bt = Bootstrap, cd = climate data record)'
	winds.description = 'Wind speed (m)'
	temps.description = '2m air temperature (C)'
	
	longitude[:] = np.around(lons, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats, decimals=4)
	snowVol[:] = np.around(snowVolT, decimals=4)
	snowDepth[:] = np.around(snowDepthT, decimals=4)
	density[:] = np.around(densityT, decimals=4)
	iceConc[:] = np.around(iceConcT, decimals=4)
	precip[:] = np.around(precipT, decimals=4)
	winds[:] = np.around(windsT, decimals=4)
	temps[:] = np.around(tempT, decimals=4)
	day[:]=datesT

	today = datetime.datetime.today()

	#Add global attributes
	f.author = "Alek Petty"
	f.contact = " alek.a.petty@nasa.gov, www.alekpetty.com, @alekpetty"
	f.description = "Daily NESOSIM snow on sea ice estimates (v1.1)"
	f.history = "Created " + today.strftime("%d/%m/%y")
	f.data_range = "Date range of the snow budgets: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def calcLeadLoss(snowDepthT, windDayT, iceConcDaysT):
	""" Snow loss to leads due to winds

	Use a variable leadlossfactor parameter. This is relatively unconstrained!

	Args:
		snowDepthT (var): Daily gridded snowdepth 
		WindDayT (var): Daily gridded wind magnitude
		iceConcDaysT (var): Daily gridded ice concentration

	returns:
		snowLeadT (var): Snow lost from fresh snow layer

	Updates:
		v1.0 (during review process) added wind packing threshold

	"""

	windT= np.where(windDayT>windPackThresh, 1, 0)
	 
	snowLeadT = -(windT*leadLossFactor*deltaT*snowDepthT*windDayT*(1-iceConcDaysT)) #*iceConcDaysG[x]
	return snowLeadT

def calcAtmLoss(snowDepthT, windDayT, atm_loss_factor=0.1):
	""" Snow lost to the atmosphere due to winds

	Using a multiple of the variable windlossfactor parameter. This is relatively unconstrained!

	Args:
		snowDepthT (var): Daily gridded snowdepth 
		WindDayT (var): Daily gridded wind magnitude
		atm_loss_factor (var): extra atm loss factor, replacing the concentration dependence of the lead loss factor

	returns:
		snowAtmLossT (var): Snow lost from fresh snow layer

	Updates:
		v1.1: testing new wind-atmosphere snow loss term

	"""

	windT= np.where(windDayT>windPackThresh, 1, 0)
	
	# half this randomly!
	snowAtmLossT = -(windT*leadLossFactor*deltaT*snowDepthT*windDayT*atm_loss_factor) #*iceConcDaysG[x]
	return snowAtmLossT

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


	windT= np.where(windDayT>windPackThresh, 1, 0)
	
	# snow loss from fresh layer through wind packing to old layer
	snowWindPackLossT=-windPackFactor*deltaT*windT*snowDepthT0 #*iceConcDaysG[x]

	# snow gain to old layer through wind packing from fresh layer
	snowWindPackGainT=windPackFactor*deltaT*windT*snowDepthT0*(snowDensityFresh/snowDensityOld) #*iceConcDaysG[x]
	
	snowWindPackNetT=snowWindPackLossT+snowWindPackGainT#*iceConcDaysG[x]
	return snowWindPackLossT, snowWindPackGainT, snowWindPackNetT

def fillMaskAndNaNWithZero(arr):
	""" Helper function: Fill masked and nan values in an array 
	with 0, in place

	Args:
		arr (var): A numpy ndarray
	returns:
		None (performs operation in place)

	"""
	arr[np.isnan(arr)] = 0.
	arr = ma.filled(arr, 0.)
	arr[~np.isfinite(arr)] = 0.

def fill_nan_no_negative(arr, region_maskG, negative_to_zero=True):
	""" Helper function: Fill masked with nan and don't allow negative snow depths if selected

	Args:
		arr (var): A numpy ndarray x, y
	returns:
		None (performs operation in place)

	"""

	

	# Set infinte snow depths to nan
	arr[~np.isfinite(arr)]=np.nan
	arr[~np.isfinite(arr)]=np.nan
	
	# Set snow depths over land/coasts to nan (changed from zero)
	arr[np.where(region_maskG>10)]=np.nan
	arr[np.where(region_maskG>10)]=np.nan

	# Set snow depths over lakes (lake sic included in CDR) to nan (changed from zero)
	arr[np.where(region_maskG<1)]=np.nan

	# Set negative snow to zero.
	if (negative_to_zero):
		arr[np.where(arr<0.)]=0.
	#arr[x+1][np.where(np.isnan(arr[x+1]))]=0.


def smooth_snow(arr, stddev_val=1,x_size_val=3, y_size_val=3):
	""" Smooths array using a Gaussian kernal (from astropy)

	Args:
		arr (var): A numpy ndarray x, y
		stddev_val (constant): gaussian width of filter
		x_size_val (constant): grid-cell width in x
		y_size_val (constant): grid-cell width in y
		
	returns:
		None (performs operation in place)
	"""
	
	# x_stddev applies to both x and y if y isn't specified.
	kernel = Gaussian2DKernel(x_stddev=stddev_val, x_size=x_size_val, y_size=x_size_val)
	arr = convolve(arr, kernel)
	
	return arr

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

	#------------  Snow change from ice divergence/convergence. Convergence is positive
	dhsvelxdxDiv = snowDepthsT*np.gradient(driftGday[0]*deltaT, dx, axis=(1)) #convert from m/s to m per day, #1 here is the columns, so in the x direction
	dhsvelydyDiv = snowDepthsT*np.gradient(driftGday[1]*deltaT, dx, axis=(0)) #0 here is the rows, so in the y direction
	# Snow divergence
	snowDivAllT= -(dhsvelxdxDiv + dhsvelydyDiv)

	#------------ Snow change from ice advection. Advetion away from a grid-cell is negative
	dhsvelxdxAdv = driftGday[0]*deltaT*np.gradient(snowDepthsT, dx, axis=(2)) #convert from m/s to m per day, #1 here is the columns, so in the x direction
	dhsvelydyAdv = driftGday[1]*deltaT*np.gradient(snowDepthsT, dx, axis=(1))  #0 here is the rows, so in the y direction
	# Snow advection
	snowAdvAllT= -(dhsvelxdxAdv + dhsvelydyAdv)

	# Set bad values to zero
	fillMaskAndNaNWithZero(snowAdvAllT[0])
	fillMaskAndNaNWithZero(snowAdvAllT[1])
	fillMaskAndNaNWithZero(snowDivAllT[0])
	fillMaskAndNaNWithZero(snowDivAllT[1])
	

	return snowAdvAllT, snowDivAllT

def calcBudget(xptsG, yptsG, snowDepths, iceConcDayT, precipDayT, driftGdayT, windDayT, tempDayT, 
	density, precipDays, iceConcDays, windDays, tempDays, snowAcc, snowOcean, snowAdv, 
	snowDiv, snowLead, snowAtm, snowWindPackLoss, snowWindPackGain, snowWindPack, region_maskG, dx, x, dayT,
	densityType='variable', dynamicsInc=1, leadlossInc=1, windpackInc=1, atmlossInc=0):
	""" Snow budget calculations

	Args:
		xptsG (var, x/y): x coordinates of grid
		yptsG (var, x/y): y coordinates of grid
		snowDepths (var, day/x/y): daily snow depth grid
		iceConcDayT (var, x/y): ice concentration for that day
		precipDayT (var, x/y): precip for that day
		driftGdayT (var, x/y): ice drift for that day
		windDayT (var, x/y): wind speed magnitude for that day
		tempDayT (var, x/y): near surface air temperture for that day
		density (var, x/y/j): 2 layer snow density

	returns:
		Updated snow budget arrays
		
	"""

	precipDays[x]=precipDayT
	iceConcDays[x]=iceConcDayT
	windDays[x]=windDayT
	tempDays[x]=tempDayT

	if (densityType=='clim'):
		# returns a fixed density value assigned to all grid cells based on climatology. 
		# Applies the same value to both snow layers.
		snowDensityNew=cF.densityClim(dayT, ancDataPath)
	else:
		# Two layers so a new snow density and an evolving old snow density
		snowDensityNew=snowDensityFresh
		
	# Convert precip to m/day
	precipDayDelta=precipDayT/snowDensityNew

	# ------ Snow accumulation
	snowAccDelta= (precipDayDelta * iceConcDayT)
	snowAcc[x+1] = snowAcc[x] + snowAccDelta

	# ------ Ocean freshwater flux
	snowOceanDelta= -(precipDayDelta * (1-iceConcDayT))
	snowOcean[x+1] = snowOcean[x] + snowOceanDelta

	# ------ Snow/ice dynamics calulation

	if (dynamicsInc==1):
		snowAdvDelta, snowDivDelta = calcDynamics(driftGdayT, snowDepths[x], dx)

		# Smooth dynamics terms as generally quite noisy
		snowAdvDelta[0]=smooth_snow(snowAdvDelta[0])
		snowAdvDelta[1]=smooth_snow(snowAdvDelta[1])
		snowDivDelta[0]=smooth_snow(snowDivDelta[0])
		snowDivDelta[1]=smooth_snow(snowDivDelta[1])
		
		fill_nan_no_negative(snowAdvDelta[0], region_maskG, negative_to_zero=False)
		fill_nan_no_negative(snowAdvDelta[1], region_maskG, negative_to_zero=False)
		fill_nan_no_negative(snowDivDelta[0], region_maskG, negative_to_zero=False)
		fill_nan_no_negative(snowDivDelta[1], region_maskG, negative_to_zero=False)

	else:
		snowAdvDelta=np.zeros((iceConcDayT.shape))
		snowDivDelta=np.zeros((iceConcDayT.shape))
	
	snowAdv[x+1] = snowAdv[x] + snowAdvDelta[0]+ snowAdvDelta[1]
	snowDiv[x+1] = snowDiv[x] + snowDivDelta[0]+ snowDivDelta[1]

	# ------ Lead loss calulation

	if (leadlossInc==1):
		snowLeadDelta= calcLeadLoss(snowDepths[x, 0], windDayT, iceConcDayT)
	else:
		snowLeadDelta=np.zeros((iceConcDayT.shape))

	snowLead[x+1]=snowLead[x] + snowLeadDelta

	# ------ Wind loss calulation

	if (atmlossInc==1):
		snowAtmDelta= calcAtmLoss(snowDepths[x, 0], windDayT)
	else:
		snowAtmDelta=np.zeros((iceConcDayT.shape))

	snowAtm[x+1]=snowAtm[x] + snowAtmDelta

	#---------- Wind packing calulation

	if (windpackInc==1):
		snowWindPackLossDelta, snowWindPackGainDelta, snowWindPackNetDelta=calcWindPacking(windDayT, snowDepths[x, 0])
	else:
		snowWindPackLossDelta =np.zeros((iceConcDayT.shape))
		snowWindPackGainDelta=np.zeros((iceConcDayT.shape))
		snowWindPackNetDelta=np.zeros((iceConcDayT.shape))

	snowWindPackLoss[x+1]=snowWindPackLoss[x]+snowWindPackLossDelta
	snowWindPackGain[x+1]=snowWindPackGain[x]+snowWindPackGainDelta
	snowWindPack[x+1]=snowWindPack[x]+snowWindPackNetDelta

	#------------ Update snow depths

	# New (upper) layer
	snowDepths[x+1, 0]=snowDepths[x, 0]+snowAccDelta  +snowWindPackLossDelta + snowLeadDelta + snowAtmDelta +snowAdvDelta[0]+snowDivDelta[0] #+snowRidgeT
	# Old snow layer
	snowDepths[x+1, 1]=snowDepths[x, 1] +snowWindPackGainDelta + snowAdvDelta[1] + snowDivDelta[1] #+ snowDcationT

	# Fill negatives and set nans
	fill_nan_no_negative(snowDepths[x+1, 0], region_maskG)
	fill_nan_no_negative(snowDepths[x+1, 1], region_maskG)
	
	if (densityType=='clim'):
		# returns a fixed density value assigned to all grid cells based on climatology. 
		# Applies the same value to both snow layers.
		
		density[x+1]=snowDensityNew
		# mask over land, lakes and coast
		density[x+1][np.where(region_maskG>10)]=np.nan
		density[x+1][np.where(region_maskG<1)]=np.nan
		density[x+1][np.where(iceConcDayT<minConc)]=np.nan
		density[x+1][np.where((snowDepths[x+1][0]+snowDepths[x+1][1])<minSnowD)]=np.nan
	else:
		# Two layers so a new snow density and an evolving old snow density	
		density[x+1]=densityCalc(snowDepths[x+1], iceConcDayT, region_maskG)


def genEmptyArrays(numDaysT, nxT, nyT):
	""" 
	Declare empty arrays to store the various budget terms

	"""
	
	precipDays=np.zeros((numDaysT, nxT, nyT)) 
	iceConcDays=np.zeros((numDaysT, nxT, nyT)) 
	windDays=np.zeros((numDaysT, nxT, nyT)) 
	tempDays=np.zeros((numDaysT, nxT, nyT)) 

	snowDepths=np.zeros((numDaysT, 2, nxT, nyT))
	density=np.zeros((numDaysT, nxT, nyT))

	snowDiv=np.zeros((numDaysT, nxT, nyT))
	snowAdv=np.zeros((numDaysT, nxT, nyT))
	snowAcc=np.zeros((numDaysT, nxT, nyT))
	snowOcean=np.zeros((numDaysT, nxT, nyT))
	snowWindPack=np.zeros((numDaysT, nxT, nyT))
	snowWindPackLoss=np.zeros((numDaysT, nxT, nyT))
	snowWindPackGain=np.zeros((numDaysT, nxT, nyT))
	snowLead=np.zeros((numDaysT, nxT, nyT))
	snowAtm=np.zeros((numDaysT, nxT, nyT))
	

	return precipDays, iceConcDays, windDays, tempDays, snowDepths, density, snowDiv, snowAdv, snowAcc, snowOcean, snowWindPack, \
	snowWindPackLoss, snowWindPackGain, snowLead, snowAtm

def plot_budgets_cartopy(lonG, latG, precipDaysT, windDaysT, snowDepthsT, snowOceanT, snowAccT, snowDivT, \
	snowAdvT, snowLeadT, snowAtmT, snowWindPackT, snowWindPackLossT, snowWindPackGainT, densityT, dateStr, totalOutStr='test'):
	""" Plot snow budget terms """

	cF.plot_gridded_cartopy(lonG, latG, precipDaysT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/precip'+totalOutStr, units_lab='kg/m2', varStr='Snowfall', minval=0., maxval=1, cmap_1=cm.cubehelix_r)
	cF.plot_gridded_cartopy(lonG, latG, windDaysT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/wind'+totalOutStr, units_lab='m/s', varStr='Wind speed', minval=0., maxval=10, cmap_1=cm.cubehelix_r)
	cF.plot_gridded_cartopy(lonG, latG, snowDepthsT[0], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowNew_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	cF.plot_gridded_cartopy(lonG, latG, snowDepthsT[1], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowOld_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	cF.plot_gridded_cartopy(lonG, latG, snowDepthsT[0]+snowDepthsT[1], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowTot_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	
	cF.plot_gridded_cartopy(lonG, latG, snowDivT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDiv_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	cF.plot_gridded_cartopy(lonG, latG, snowAdvT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowAdv_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	cF.plot_gridded_cartopy(lonG, latG, snowDivT+snowAdvT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDyn_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	cF.plot_gridded_cartopy(lonG, latG, snowLeadT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowLeadLoss_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	cF.plot_gridded_cartopy(lonG, latG, snowAtmT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowAtmLoss_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	cF.plot_gridded_cartopy(lonG, latG, snowWindPackT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowWindPack_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)

	cF.plot_gridded_cartopy(lonG, latG, densityT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDensity_'+totalOutStr, units_lab='kg/m3', varStr='Snow density', minval=220, maxval=340, cmap_1=cm.viridis)


def loadData(yearT, dayT, precipVar, windVar, concVar, driftVar, dxStr, extraStr):
	""" Load daily forcings

	Temp added transpose to convert grid to proper row/column index. 

	"""
	dayStr='%03d' %dayT
	
	#------- Read in precipitation -----------
	try:
		precipDayG=np.load(forcingPath+'Precip/'+precipVar+'/'+str(yearT)+'/'+precipVar+'sf'+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr, allow_pickle=True).T
	except:
		if (dayStr=='365'):
			print('no leap year data, using data from the previous day')
			precipDayG=np.load(forcingPath+'Precip/'+precipVar+'/sf/'+str(yearT)+'/'+precipVar+'sf'+dxStr+'-'+str(yearT)+'_d'+'364', allow_pickle=True).T
		
		else:
			print('No precip data so exiting!')
			exit()
	
	#------- Read in wind magnitude -----------
	try:
		windDayG=np.load(forcingPath+'Winds/'+windVar+'/'+str(yearT)+'/'+windVar+'winds'+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr, allow_pickle=True).T
	
	except:
		if (dayStr=='365'):
			print('no leap year data, using data from the previous day')
			windDayG=np.load(forcingPath+'Winds/'+windVar+'/'+str(yearT)+'/'+windVar+'winds'+dxStr+'-'+str(yearT)+'_d'+'364'+extraStr, allow_pickle=True).T
		
		else:
			print('No wind data so exiting!')
			exit()

	#------- Read in ice concentration -----------
	try:
		iceConcDayG=np.load(forcingPath+'IceConc/'+concVar+'/'+str(yearT)+'/iceConcG_'+concVar+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr, allow_pickle=True).T
	
	except:
		if (dayStr=='365'):
			print('no leap year data, using data from the previous day')
			iceConcDayG=np.load(forcingPath+'IceConc/'+concVar+'/'+str(yearT)+'/iceConcG_'+concVar+dxStr+'-'+str(yearT)+'_d'+'364'+extraStr, allow_pickle=True).T
	
		else:
			print('No ice conc data so exiting!')
			exit()

	# fill with zero
	iceConcDayG[~np.isfinite(iceConcDayG)]=0.
	
	#------- Read in ice drifts -----------
	try:
		driftGdayG=np.load(forcingPath+'IceDrift/'+driftVar+'/'+str(yearT)+'/'+driftVar+'_driftG'+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr+'_n', allow_pickle=True)	
	except:
		
		# if no drifts exist for that day then just set drifts to nan array (i.e. no drift).
		print('No drift data')
		driftGdayG = np.empty((2, iceConcDayG.shape[0], iceConcDayG.shape[1]))
		driftGdayG[:] = np.nan

	driftGdayG = ma.filled(driftGdayG, np.nan)
	print(driftGdayG)

	#------- Read in temps (not currently used, placeholder) -----------
	try:
		tempDayG=np.load(forcingPath+'Temp/'+precipVar+'/t2m/'+str(yearT)+'/t2m'+dxStr+'-'+str(yearT)+'_d'+dayStr+extraStr, allow_pickle=True)
	except:
		# if no drifts exist for that day then just set drifts to masked array (i.e. no drift).
		#print('No temp data')
		tempDayG = np.empty((iceConcDayG.shape[0], iceConcDayG.shape[1]))
		tempDayG[:] = np.nan
		#tempDayG=ma.masked_all((iceConcDayG.shape[0], iceConcDayG.shape[1]))
	
	return iceConcDayG, precipDayG, driftGdayG, windDayG, tempDayG

def densityCalc(snowDepthsT, iceConcDayT, region_maskT):
	"""Assign initial density based on snow depths

	dropped ice conc mask as based on old time step
	"""		

	densityT=((snowDepthsT[0]*snowDensityFresh) + (snowDepthsT[1]*snowDensityOld))/(snowDepthsT[0]+snowDepthsT[1]) #+ densDcationT
	
	densityT[np.where(densityT>snowDensityOld)]=snowDensityOld
	densityT[np.where(densityT<snowDensityFresh)]=snowDensityFresh

	densityT[np.where(region_maskT<1)]=np.nan
	densityT[np.where(region_maskT>10)]=np.nan
	densityT[np.where((snowDepthsT[0]+snowDepthsT[1])<minSnowD)]=np.nan

	return densityT

def doyToMonth(day, year):
	""" given the day-of-year day and the year, return an integer
	corresponding to the month during which the day occurs
	"""
	date_fmt = np.datetime64('{}-01-01'.format(year)) + np.timedelta64(day-1,'D')
	return date_fmt.astype(object).month

def applyScaling(product,factor,scaling_type='mul'):
	"""Apply a scaling factor to a given product; the factor
	must either be a scalar or have the same dimensions as
	the product
	
	"""
	if scaling_type=='mul':
		# multiplicative scaling
		product_scaled = product*factor

	return product_scaled


def main(year1, month1, day1, year2, month2, day2, outPathT='.', forcingPathT='.', anc_data_pathT='../anc_data/', figPathT='../Figures/', 
	precipVar='ERA5', windVar='ERA5', driftVar='OSISAF', concVar='CDR', icVar='ERAI', densityTypeT='variable', 
	outStr='', extraStr='', IC=2, windPackFactorT=0.1, windPackThreshT=5., leadLossFactorT=0.1, dynamicsInc=1, leadlossInc=1, 
	windpackInc=1, atmlossInc=0, saveData=1, plotBudgets=1, plotdaily=1, saveFolder='', dx=50000,scaleCS=False):
	""" 

	Main model function

	Args:
		The various model configuration parameters

	"""

	#------- Create map projection
	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	nx=xptsG.shape[0]
	ny=xptsG.shape[1]

	dxStr=str(int(dx/1000))+'km'
	print(nx, ny, dxStr)

	# Assign some global parameters
	global dataPath, forcingPath, outPath, ancDataPath
	
	outPath=outPathT+dxStr+'/'
	forcingPath=forcingPathT+dxStr+'/'
	ancDataPath=anc_data_pathT
	print('OutPath:', outPath)
	print('forcingPath:', forcingPath)
	print('ancDataPath:', ancDataPath)

	# Assign density of the two snow layers
	global snowDensityFresh, snowDensityOld, minSnowD, minConc, leadLossFactor, windPackThresh, windPackFactor, deltaT
	snowDensityFresh=200. # density of fresh snow layer
	snowDensityOld=350. # density of old snow layer
	minSnowD=0.02 # minimum snow depth for a density estimate
	minConc=0.15 # mask budget values with a concentration below this value

	deltaT=60.*60.*24. # time interval (seconds in a day)

	region_mask, xptsI, yptsI = cF.get_region_mask_pyproj(anc_data_pathT, proj, xypts_return=1)
	region_maskG = griddata((xptsI.flatten(), yptsI.flatten()), region_mask.flatten(), (xptsG, yptsG), method='nearest')

	leadLossFactor=leadLossFactorT # Snow loss to leads coefficient
	windPackThresh=windPackThreshT # Minimum winds needed for wind packing
	windPackFactor=windPackFactorT # Fraction of snow packed into old snow layer

	# Current year
	yearCurrent=year1
	
	#--------- Get time period info
	startDay, numDays, numDaysYear1, dateOut= cF.getDays(year1, month1, day1, year2, month2, day2)
	print (startDay, numDays, numDaysYear1, dateOut)

	# make this into a small function
	dates=[]
	for x in range(0, numDays):
		#print x
		date = datetime.datetime(year1, month1+1, day1+1) + datetime.timedelta(x)
		#print (int(date.strftime('%Y%m%d')))
		dates.append(int(date.strftime('%Y%m%d')))
	#print(dates)
	
	CSstr = ''
	if scaleCS:
		# load scaling factors; assumes scaling factors are in same directory as NESOSIM.py
		monthlyScalingFactors = xr.open_dataset('{}scale_coeffs_{}.nc'.format(ancDataPath, precipVar))['scale_factors']
		CSstr = 'CSscaled'

	#------ create output strings and file paths -----------
	saveStr= precipVar+CSstr+'sf'+windVar+'winds'+driftVar+'drifts'+concVar+'sic'+'rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_AL'+str(atmlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+extraStr+outStr+'-'+dateOut
	saveStrNoDate=precipVar+CSstr+'sf'+windVar+'winds'+driftVar+'drifts'+concVar+'sic'+'rho'+densityTypeT+'_IC'+str(IC)+'_DYN'+str(dynamicsInc)+'_WP'+str(windpackInc)+'_LL'+str(leadlossInc)+'_AL'+str(atmlossInc)+'_WPF'+str(windPackFactorT)+'_WPT'+str(windPackThreshT)+'_LLF'+str(leadLossFactorT)+'-'+dxStr+extraStr+outStr
	
	print ('Saving to:', saveStr)
	 #'../../DataOutput/'

	savePath=outPath+saveFolder+'/'+saveStrNoDate
	# Declare empty arrays for compiling budgets
	if not os.path.exists(savePath+'/budgets/'):
		os.makedirs(savePath+'/budgets/')
	if not os.path.exists(savePath+'/final/'):
		os.makedirs(savePath+'/final/')

	global figpath
	figpath=figPathT+'/Diagnostic/'+dxStr+'/'+saveStrNoDate+'/'
	if not os.path.exists(figpath):
		os.makedirs(figpath)
	if not os.path.exists(figpath+'/daily_snow_depths/'):
		os.makedirs(figpath+'/daily_snow_depths/')

	precipDays, iceConcDays, windDays, tempDays, snowDepths, density, snowDiv, snowAdv, snowAcc, snowOcean, snowWindPack, snowWindPackLoss, snowWindPackGain, snowLead, snowAtm = genEmptyArrays(numDays, nx, ny)

	print('IC:', IC)
	if (IC>0):
		if (IC==1):
			# August Warren climatology snow depths
			ICSnowDepth = np.load(forcingPath+'InitialConditions/AugSnow'+dxStr, allow_pickle=True)
			print('Initialize with August Warren climatology')
		
		elif (IC==2):
			# Alek v2 (capped at 10 m) ICs based on MW method
			try:
				# Convert to meters!!
				ICSnowDepth = np.load(forcingPath+'InitialConditions/'+icVar+'/ICsnow'+str(year1)+'-'+dxStr+extraStr, allow_pickle=True)
				print('Initialize with new v1.1 scaled initial conditions')
				print(np.amax(ICSnowDepth))
			except:
				print('No initial conditions file available')

		iceConcDayG, precipDayG, driftGdayG, windDayG, tempDayG =loadData(year1, startDay, precipVar, windVar, concVar, driftVar, dxStr, extraStr)
		ICSnowDepth[np.where(iceConcDayG<minConc)]=0

		# Split the initial snow depth over both layers
		snowDepths[0, 0]=ICSnowDepth*0.5
		snowDepths[0, 1]=ICSnowDepth*0.5


	#pF.plotSnow(m, xptsG, yptsG, densityT, date_string=str(startDay-1), out=figpath+'/Snow/2layer/densityD'+driftP+extraStr+reanalysisP+varStr+'_sy'+str(year1)+'d'+str(startDay)+outStr+'T0', units_lab=r'kg/m3', minval=180, maxval=360, base_mask=0, norm=0, cmap_1=cm.viridis)

	# Loop over days 
	for x in range(numDays-1):	
		day = x+startDay
		

		if (day>=numDaysYear1):
			# If day goes beyond the number of days in initial year, jump to the next year
			day=day-numDaysYear1
			yearCurrent=year2
		
		print ('day:', day)
		
		# Load daily data 
		iceConcDayG, precipDayG, driftGdayG, windDayG, tempDayG =loadData(yearCurrent, day, precipVar, windVar, concVar, driftVar, dxStr, extraStr)
		
		# apply CloudSat scaling if used
		if scaleCS:
			currentMonth = doyToMonth(day, yearCurrent) # get current month
			scalingFactor = monthlyScalingFactors.loc[currentMonth,:,:] # get scaling factor for current month
			# apply scaling to current day's precipitation
			precipDayG = applyScaling(precipDayG, scalingFactor,scaling_type='mul').values

		# Calculate snow budgets
		calcBudget(xptsG, yptsG, snowDepths, iceConcDayG, precipDayG, driftGdayG, windDayG, tempDayG,
			density, precipDays, iceConcDays, windDays, tempDays, snowAcc, snowOcean, snowAdv, 
			snowDiv, snowLead, snowAtm, snowWindPackLoss, snowWindPackGain, snowWindPack, region_maskG, dx, x, day,
			densityType=densityTypeT, dynamicsInc=dynamicsInc, leadlossInc=leadlossInc, windpackInc=windpackInc, atmlossInc=atmlossInc)
		
		if (plotdaily==1):
			cF.plot_gridded_cartopy(lonG, latG, snowDepths[x+1, 0]+snowDepths[x+1, 1], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'daily_snow_depths/snowTot_'+saveStrNoDate+str(x), units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	
	# Load last data 
	iceConcDayG, precipDayG, _, windDayG, tempDayG =loadData(yearCurrent, day+1, precipVar, windVar, concVar, driftVar, dxStr, extraStr)
	precipDays[x+1]=precipDayG
	iceConcDays[x+1]=iceConcDayG
	windDays[x+1]=windDayG
	tempDays[x+1]=tempDayG
	
	outStrings=['snowDepthTotal','snowDepthTotalConc', 'density', 'iceConc', 'Precip']

	if (saveData==1):
		# Output snow budget terms to netcdf datafiles
		OutputSnowModelRaw(savePath, saveStr, snowDepths, density, precipDays, iceConcDays, windDays, snowAcc, snowOcean, snowAdv, snowDiv, snowLead, snowAtm, snowWindPack)
		OutputSnowModelFinal(savePath, saveStr, lonG, latG, snowDepths[:, 0]+snowDepths[:, 1], (snowDepths[:, 0]+snowDepths[:, 1])/iceConcDays, density, iceConcDays, precipDays, windDays, tempDays, dates)

	if (plotBudgets==1):
		# Plot final snow budget terms 
		plot_budgets_cartopy(lonG, latG, precipDayG, windDayG, snowDepths[x+1], snowOcean[x+1], snowAcc[x+1], snowDiv[x+1], \
		snowAdv[x+1], snowLead[x+1], snowAtm[x+1], snowWindPack[x+1], snowWindPackLoss[x+1], snowWindPackGain[x+1], density[x+1], dates[-1], totalOutStr=saveStr)

