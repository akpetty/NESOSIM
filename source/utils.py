""" utils.py
	
	Common functions used by the NESOSIM.py script 
	Original code written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)


	Python dependencies:
		See below for the relevant module imports. Of note:
		matplotlib
		basemap

	Update history:
		03/01/2018: Version 1
		05/10/2020: Version 1.1: Converted to Python 3
								Changed name to utils.py

"""

from glob import glob
import matplotlib.pyplot as plt
import matplotlib.colorbar as mcbar
from scipy.interpolate import griddata
import xarray as xr
from scipy import stats
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import pyproj
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
import datetime
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
import pandas as pd
import netCDF4 as nc4
import matplotlib.cm as cm
import pandas as pd

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


def OutputSnowModelFinal(savePath, saveStr, lons, lats, xpts, ypts, snowVolT,snowDepthT, densityT, iceConcT, precipT, windsT, tempT, datesT, ice_conc_mask=0.5):
	""" Read in xrarray data and save as netCDF 

	Args:
		savePath (str): Path the the xarray data wil be saved to
		reanalysisP (str): Reanalysis snowfall forcing used for this model run
		saveStr (str): output string for saved filed
		Remaining arguments* (vars): Model variables being saved  
    
    """

	f = nc4.Dataset(savePath+'/final/'+saveStr+'.nc','w', format='NETCDF4') 

	projection = f.createVariable('projection', 'i4')
	projection.long_name = "WGS 84 / NSIDC Sea Ice Polar Stereographic North (3413)" ;
	projection.spatial_ref = "PROJCS[\"WGS 84 / NSIDC Sea Ice Polar Stereographic North\",GEOGCS[\"WGS 84\"[\"DATUM[\"WGS_1984\"[\"SPHEROID[\"WGS 84\",6378137,298.257223563[\"AUTHORITY[\"EPSG\",\"7030\"]][\"AUTHORITY[\"EPSG\",\"6326\"]][\"PRIMEM[\"Greenwich\",0[\"AUTHORITY[\"EPSG\",\"8901\"]][\"UNIT[\"degree\",0.0174532925199433[\"AUTHORITY[\"EPSG\",\"9122\"]][\"AUTHORITY[\"EPSG\",\"4326\"]][\"PROJECTION[\"Polar_Stereographic\"][\"PARAMETER[\"latitude_of_origin\",70][\"PARAMETER[\"central_meridian\",-45][\"PARAMETER[\"scale_factor\",1][\"PARAMETER[\"false_easting\",0][\"PARAMETER[\"false_northing\",0][\"UNIT[\"metre\",1[\"AUTHORITY[\"EPSG\",\"9001\"]][\"AXIS[\"X\",EAST][\"AXIS[\"Y\",NORTH][\"AUTHORITY[\"EPSG\",\"3413\"]]";
	projection.grid_mapping_name = "polar_stereographic"
	projection.proj4text = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs";
	
	print ('dimensions:', lons.shape[0], lons.shape[1], snowVolT.shape[0])
	f.createDimension('x', lons.shape[0])
	f.createDimension('y', lons.shape[1])
	f.createDimension('day', snowVolT.shape[0])

	longitude = f.createVariable('longitude', 'f4', ('x', 'y'))
	latitude = f.createVariable('latitude', 'f4', ('x', 'y'))  
	xgrid = f.createVariable('xgrid', 'f4', ('x', 'y'))
	ygrid = f.createVariable('ygrid', 'f4', ('x', 'y'))  

	snowDepth = f.createVariable('snow_depth', 'f4', ('day', 'x', 'y'))
	snowVol = f.createVariable('snow_volume', 'f4', ('day', 'x', 'y'))
	density = f.createVariable('snow_density', 'f4', ('day', 'x', 'y'))
	precip = f.createVariable('precipitation', 'f4', ('day', 'x', 'y'))
	iceConc = f.createVariable('ice_concentration', 'f4', ('day', 'x', 'y'))
	winds = f.createVariable('wind_speed', 'f4', ('day', 'x', 'y'))
	#temps = f.createVariable('airTemp2m', 'f4', ('day', 'x', 'y'))
	day = f.createVariable('day', 'i4', ('day'))

	day.units = 'yyyymmdd'
	day.long_name='calendar date'

	longitude.units = 'degrees East'
	longitude.long_name='longitude'

	xgrid.units = 'meters'
	xgrid.long_name='projection grid x values'
	xgrid.description = "center values of projection grid in x direction" 

	ygrid.units = 'meters'
	ygrid.long_name='projection grid y values'
	ygrid.description = "center values of projection grid in y direction" 

	latitude.units = 'degrees North'
	latitude.long_name='latitude'

	snowVol.description = 'Daily snow volume per unit grid cell'
	snowVol.units = 'meters'
	snowVol.long_name = 'snow volume'

	snowDepth.description = 'Daily snow depth (effective over the ice fraction)'
	snowDepth.units = 'meters'
	snowDepth.long_name = 'snow depth'

	density.description = 'Bulk snow density'
	density.units = 'kilograms per meters cubed'
	density.long_name = 'snow density'

	precip.description = 'Precipitation, generally the explicit snowfall component of total precipitaiton provided by the chosen reanalysis'
	precip.units = 'kilograms of snow per meters squared'
	precip.long_name = 'Snowfall'

	iceConc.description = 'Sea ice concentration derived from passive mcirowave concentration'
	iceConc.units = 'unitless (between 0 and 1)'
	iceConc.long_name = 'ice concentration'

	winds.description = 'Wind speed magnitude, calculated as the root mean square of the u/v wind vectors'
	winds.units = 'meters'
	winds.long_name = 'wind speed'
	
	if ice_conc_mask>0:
		snowVolT[np.where(iceConcT<ice_conc_mask)]=np.nan 
		snowDepthT[np.where(iceConcT<ice_conc_mask)]=np.nan 
		densityT[np.where(iceConcT<ice_conc_mask)]=np.nan 

		# Also mask ice concentration less than 0.15 due to knwon uncertainties with the passive microwave data
		iceConcT[np.where(iceConcT<0.15)]=np.nan 


	longitude[:] = np.around(lons, decimals=4) #The "[:]" at the end of the variable instance is necessary
	latitude[:] = np.around(lats, decimals=4)
	xgrid[:] = np.around(xpts, decimals=4) #The "[:]" at the end of the variable instance is necessary
	ygrid[:] = np.around(ypts, decimals=4)

	snowVol[:] = np.around(snowVolT, decimals=4)
	snowDepth[:] = np.around(snowDepthT, decimals=4)
	density[:] = np.around(densityT, decimals=4)
	iceConc[:] = np.around(iceConcT, decimals=4)
	precip[:] = np.around(precipT, decimals=4)
	winds[:] = np.around(windsT, decimals=4)
	#temps[:] = np.around(tempT, decimals=4)
	day[:]=datesT

	today = datetime.datetime.today()

	#Add global attributes
	
	f.reference = "github.com/akpetty/NESOSIM"
	f.contact = "alek.a.petty@nasa.gov"
	f.description = "Daily snow on sea ice (depth and density) from the NASA Eulerian Snow on Sea Ice Model (NESOSIM) version 1.1"
	f.history = "Created " + today.strftime("%d/%m/%y")
	f.data_range = "Date range: "+str(datesT[0])+'-'+str(datesT[-1])

	f.close()

def plot_budgets_cartopy(lonG, latG, precipDaysT, windDaysT, snowDepthsT, snowOceanT, snowAccT, snowDivT, \
	snowAdvT, snowLeadT, snowAtmT, snowWindPackT, snowWindPackLossT, snowWindPackGainT, densityT, dateStr, figpath, totalOutStr='test'):
	""" Plot snow budget terms """

	plot_gridded_cartopy(lonG, latG, precipDaysT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/precip'+totalOutStr, units_lab='kg/m2', varStr='Snowfall', minval=0., maxval=1, cmap_1=cm.cubehelix_r)
	plot_gridded_cartopy(lonG, latG, windDaysT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/wind'+totalOutStr, units_lab='m/s', varStr='Wind speed', minval=0., maxval=10, cmap_1=cm.cubehelix_r)
	plot_gridded_cartopy(lonG, latG, snowDepthsT[0], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowNew_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	plot_gridded_cartopy(lonG, latG, snowDepthsT[1], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowOld_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	plot_gridded_cartopy(lonG, latG, snowDepthsT[0]+snowDepthsT[1], proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowTot_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=0., maxval=0.6, cmap_1=cm.cubehelix_r)
	
	plot_gridded_cartopy(lonG, latG, snowDivT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDiv_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	plot_gridded_cartopy(lonG, latG, snowAdvT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowAdv_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	plot_gridded_cartopy(lonG, latG, snowDivT+snowAdvT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDyn_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	plot_gridded_cartopy(lonG, latG, snowLeadT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowLeadLoss_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	plot_gridded_cartopy(lonG, latG, snowAtmT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowAtmLoss_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)
	plot_gridded_cartopy(lonG, latG, snowWindPackT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowWindPack_'+totalOutStr, units_lab='m', varStr='Snow depth', minval=-0.3, maxval=0.3, cmap_1=cm.RdBu)

	plot_gridded_cartopy(lonG, latG, densityT, proj=ccrs.NorthPolarStereo(central_longitude=-45), date_string='', out=figpath+'/snowDensity_'+totalOutStr, units_lab='kg/m3', varStr='Snow density', minval=220, maxval=340, cmap_1=cm.viridis)


def bin_oib(xptsOIB, yptsOIB, xptsG, yptsG, oibVar):
	""" Bin data using numpy histogram"""

	xbins = yptsG[:, 0]+(dx/2) # these are the columns which are actually the y values
	ybins = xptsG[0]+(dx/2) 
	xbins=np.append(xbins, xbins[-1]+dx)
	ybins=np.append(ybins, ybins[-1]+dx)

	denominator, xedges, yedges = np.histogram2d(xptsOIB, yptsOIB,bins=(xbins, ybins))
	nominator, _, _ = np.histogram2d(xptsOIB, yptsOIB,bins=(xbins, ybins), weights=oibVar)
	oibG = nominator / denominator
	
	# transpose to go from columns/rows to x/y
	oibG=oibG.T
	return oibG

def getLeapYr(year):
	leapYrs=[1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]
	if year in leapYrs:
		numDays=366
	else:
		numDays=365
	return numDays


def getDays(year1, month1, day1, year2, month2, day2):
	"""Get days in model time period
	"""
	numDaysYear1=getLeapYr(year1)

	dT01 = datetime.datetime(year1, 1, 1)
	d1 = datetime.datetime(year1, month1+1, day1+1)
	d2 = datetime.datetime(year2, month2+1, day2+1)
	startDayT=(d1 - dT01).days
	numDaysT=(d2 - d1).days+1

	fmt = '%d%m%Y'
	date1Str=d1.strftime(fmt)
	date2Str=d2.strftime(fmt)

	return startDayT, numDaysT, numDaysYear1, date1Str+'-'+date2Str


def int_smooth_drifts_v2(xptsG, yptsG, xptsF, yptsF, latsF, driftFmon, sigma_factor=1, x_size_val=3, truncate=1):
	# use info from https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
	# to better rapply to nan regions
	nx=xptsG.shape[0]
	ny=xptsG.shape[1]

	driftFG=ma.masked_all((2, nx, ny))

	# Bodge way of making a mask
	badData = np.zeros((xptsF.shape))
	# & (latsF>88)) if we want to mask the pole hole
	badData[np.where(np.isnan(driftFmon[1])& (latsF>90))]=1

	driftFx = driftFmon[0][np.where(badData<0.5)]
	driftFy = driftFmon[1][np.where(badData<0.5)]
	xptsFM = xptsF[np.where(badData<0.5)]
	yptsFM = yptsF[np.where(badData<0.5)]

	driftFGx = griddata((xptsFM, yptsFM), driftFx, (xptsG, yptsG), method='linear')
	driftFGy = griddata((xptsFM, yptsFM), driftFy, (xptsG, yptsG), method='linear')

	#driftFGx[np.isnan(driftFGx)]=0
	#driftFGy[np.isnan(driftFGy)]=0

	kernel = Gaussian2DKernel(x_stddev=sigma_factor, y_stddev=sigma_factor,x_size=x_size_val, y_size=x_size_val)

	driftFGxg = convolve(driftFGx, kernel)
	driftFGyg = convolve(driftFGy, kernel)

	# Mask based on original gridded data
	driftFG[0]=ma.masked_where(np.isnan(driftFGx), driftFGxg)
	driftFG[1]=ma.masked_where(np.isnan(driftFGy), driftFGyg)
	

	return driftFG

def getOIBbudgetDays(datesOIBT, yearT, monthT, dayT):
	"""Get number of days from start month of model run for the OIB dates"""

	d0 = datetime.datetime(yearT, monthT+1, dayT+1)
	fmt = '%Y%m%d'
	OIBdatetimes = [datetime.datetime.strptime(s, fmt) for s in datesOIBT]
	OIBnumDays= [(d1 - d0).days for d1 in OIBdatetimes]
	return OIBnumDays

def getOSISAFDrift(m, fileT):
	"""
	Calculate the OSI-SAF vectors on our map projection
	With help from Thomas Lavergne.

	"""

	f = Dataset(fileT, 'r')
	print(fileT)

	# read lat/lon (start points) and lat1/lon1 (end points)
	lon = (f.variables['lon'][:])
	lon1 = (f.variables['lon1'][0])
	lat = (f.variables['lat'][:])
	lat1 = (f.variables['lat1'][0])
	f.close()

	# transform to map project coordinates (basemap's axes, assume they have unit of m)
	x0, y0=m(lon, lat)
	x1, y1=m(lon1, lat1)

	xpts=(x0+x1)/2.
	ypts=(y0+y1)/2.

	# normalize drift components to m/s (x1-x0 is already m, so we just divide by 2-days worth of seconds)
	xt=(x1-x0)/(60*60*24*2.)
	yt=(y1-y0)/(60*60*24*2.)

	# TL: no need to rotate : the xt, and yt are already in the basemap's projection

	# compute magnitude (speed scalar)
	mag=np.sqrt(xt**2+yt**2)
	#print mag.mean(), mag.min(), mag.max()

	return xt, yt, mag, lat, lon, xpts, ypts

import cartopy.crs as ccrs


class P3413(ccrs.Projection):
	# 3413 projection in CRS format
    
    def __init__(self):

        # see: http://www.spatialreference.org/ref/epsg/3408/
        proj4_params = {'proj': 'stere',
            'lat_0': 90.,
            'lon_0': -45,
            'x_0': 0,
            'y_0': 0,
            'units': 'm',
            'no_defs': ''}

        super(P3413, self).__init__(proj4_params)

    @property
    def boundary(self):
        coords = ((self.x_limits[0], self.y_limits[0]),(self.x_limits[1], self.y_limits[0]),
                  (self.x_limits[1], self.y_limits[1]),(self.x_limits[0], self.y_limits[1]),
                  (self.x_limits[0], self.y_limits[0]))

        return ccrs.sgeom.Polygon(coords).exterior

    @property
    def threshold(self):
        return 1e5

    @property
    def x_limits(self):
        return (-2353926.81, 2345724.36)
        #return (-9000000, 9000000)

    @property
    def y_limits(self):
        return (-382558.89, 383896.60)
        #return (-9000000, 9000000)

    def describe(self):
        for k, v in vars(self).items():
            print (f'{k}: {v}')


class EASE_North(ccrs.Projection):
	# Original EASE North projection class
    '''Projection class for NSIDC EASE grid north'''
    
    def __init__(self):

        # see: http://www.spatialreference.org/ref/epsg/3408/
        proj4_params = {'proj': 'laea',
            'lat_0': 90.,
            'lon_0': 0,
            'x_0': 0,
            'y_0': 0,
            'a': 6371228,
            'b': 6371228,
            'units': 'm',
            'no_defs': ''}

        super(EASE_North, self).__init__(proj4_params)

    @property
    def boundary(self):
        coords = ((self.x_limits[0], self.y_limits[0]),(self.x_limits[1], self.y_limits[0]),
                  (self.x_limits[1], self.y_limits[1]),(self.x_limits[0], self.y_limits[1]),
                  (self.x_limits[0], self.y_limits[0]))

        return ccrs.sgeom.Polygon(coords).exterior

    @property
    def threshold(self):
        return 1e5

    @property
    def x_limits(self):
        return (-9030575.88125, 9030575.88125)
        #return (-9000000, 9000000)

    @property
    def y_limits(self):
        return (-9030575.88125, 9030575.88125)
        #return (-9000000, 9000000)

    def describe(self):
        for k, v in vars(self).items():
            print (f'{k}: {v}')

            


def get_osisaf_drifts_proj(proj, fileT):
	"""
	Calculate the OSI-SAF vectors on our given map projection
	With help from Thomas Lavergne!

	v2: based on getOSISAFdrift but using pyproj instead of basemap 

	"""
	f = Dataset(fileT, 'r')
	print(fileT)

	# read lat/lon (start points) and lat1/lon1 (end points)
	lon = (f.variables['lon'][:])
	lon1 = (f.variables['lon1'][0])
	lat = (f.variables['lat'][:])
	lat1 = (f.variables['lat1'][0])
	f.close()

	lond = lon1-lon
	latd = lat1-lat
	
	# transform to map project coordinates (basemap's axes, assume they have unit of m)
	x0, y0=proj(lon, lat)
	x1, y1=proj(lon1, lat1)

	xpts=(x0+x1)/2.
	ypts=(y0+y1)/2.

	# normalize drift components to m/s (x1-x0 is already m, so we just divide by 2-days worth of seconds)
	xt=(x1-x0)/(60*60*24*2.)
	yt=(y1-y0)/(60*60*24*2.)

	# TL: no need to rotate : the xt, and yt are already in the basemap's projection

	#mag_scale = np.hypot(u_map, v_map) / np.hypot(u, v)

	# compute magnitude (speed scalar)
	#mag_scale = np.hypot(xt, yt) / np.hypot(lond/(60*60*24*2.), latd/(60*60*24*2.))
	#xt=xt/mag_scale
	#yt=yt/mag_scale

	return xt, yt, lat, lon, xpts, ypts

def getFowlerdriftMonthV3(rawdatapath, year, month, m, mean=0):

	fowlerPath = rawdatapath+'/ICE_DRIFT/FOWLER/V3/DAILY/'

	latsF, lonsF, xptsF, yptsF = getFowlerLonLat(m, rawdatapath)


	timeIndex = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
	if year in [1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]:
		print ('LEAP YEAR')
		timeIndex = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
	
	files = glob(fowlerPath+str(year)+'/*.bin')

	#files = glob(fowlerPath+str(year)+'/*.bin')
	files= files[timeIndex[month]:timeIndex[month+1]]
	print ('Num of days:', np.size(files))
	driftFmon=ma.masked_all((np.size(files), 2, lonsF.shape[0], lonsF.shape[1]))

	
	x=0
	for file in files:

		uvelT, vvelT=getFowlerDrift(file,lonsF)

		xvel,yvel = m.rotate_vector(uvelT,vvelT,lonsF,latsF)

		xvel[np.where(ma.getmask(xvel))]=np.nan
		yvel[np.where(ma.getmask(yvel))]=np.nan
		driftFmon[x, 0]=xvel
		driftFmon[x, 1]=yvel
		x+=1

	if (mean==1):
		driftFmon=ma.mean(driftFmon, axis=0)

	return xptsF, yptsF, driftFmon, lonsF, latsF

def getFowlerLonLat(mF, rawdatapath):
	# need to use Fowler map to ensure drifts are orientated correctly
	fowlerPath = rawdatapath+'/ICE_DRIFT/FOWLER/'
	lonlatF = loadtxt(fowlerPath+'/north_x_y_lat_lon.txt')
	lonsF = np.np.reshape(lonlatF[:, 3], (361, 361))
	latsF = np.np.reshape(lonlatF[:, 2], (361, 361))
	xptsF, yptsF = mF(lonsF, latsF)

	return latsF, lonsF, xptsF, yptsF

def get_day_concSN_NRT(datapath, year, month, day, alg=0, pole='A', vStr='v1.1', mask=1, maxConc=0, lowerConc=0, monthMean=0):
	if (alg==0):
		team = 'NASA_TEAM'
		team_s = 'nt'
		header = 300
		datatype='uint8'
		scale_factor=250.
	if (alg==1):
		team = 'BOOTSTRAP'
		team_s = 'bt'
		header = 0
		datatype='<i2'
		scale_factor=1000.
	
	if (pole=='A'):
		poleStr='ARCTIC'
		rows=448
		cols=304
	if (pole=='AA'):
		poleStr='ANTARCTIC'
		rows=332
		cols=316

	month_str = '%02d' % (month+1)
	day_str = '%02d' % (day+1)
	year_str=str(year)
	
	files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+day_str+'*')
	if (np.size(files)>0):
		print('Same day conc file exists:')

	if (np.size(files)==0):
		# first try day before
		day_str = '%02d' % (day)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+day_str+'*')
		if (np.size(files)>0):
			print('Using day before file:')

	# If still nothing try day after
	if (np.size(files)==0):
		# try day after
		day_str = '%02d' % (day+2)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+day_str+'*')
		if (np.size(files)>0):
			print('Using day after file:')

	
	fd = open(files[0], 'r')
	data = np.fromfile(file=fd, dtype=datatype)
	data = data[header:]
	#FIRST 300 FILES ARE HEADER INFO
	ice_conc = np.reshape(data, [rows, cols])
	#divide by 250 to express in concentration
	ice_conc = ice_conc/scale_factor
	#GREATER THAN 250 is mask/land etc

	if (mask==1):
		ice_conc = ma.masked_where(ice_conc>1., ice_conc)
	
	if (maxConc==1):
		ice_conc = ma.where(ice_conc>1.,0, ice_conc)

	if (lowerConc==1):
		ice_conc = ma.where(ice_conc<0.15,0, ice_conc)

	if (monthMean==1):
		ice_conc=ma.mean(ice_conc, axis=0)

	return ice_conc

def get_month_concSN_daily(datapath, year, month, alg=0, pole='A', vStr='v1.1', mask=1, maxConc=0, lowerConc=0, monthMean=0):
	if (alg==0):
		team = 'NASA_TEAM'
		team_s = 'nt'
		header = 300
		datatype='uint8'
		scale_factor=250.
	if (alg==1):
		team = 'BOOTSTRAP'
		team_s = 'bt'
		header = 0
		datatype='<i2'
		scale_factor=1000.
	
	if (pole=='A'):
		poleStr='ARCTIC'
		rows=448
		cols=304
	if (pole=='AA'):
		poleStr='ANTARCTIC'
		rows=332
		cols=316

	month_str = '%02d' % (month+1)
	year_str=str(year)
	files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/daily/'++str(year)+'/'+team_s+'_'+str(year)+month_str+'*'+vStr+'*')
	

	print('Num conc files:', np.size(files), 'in month:'+month_str)
	ice_conc = ma.masked_all((np.size(files), rows, cols))
	
	for x in range(np.size(files)):
		fd = open(files[x], 'r')
		data = np.fromfile(file=fd, dtype=datatype)
		data = data[header:]
		#FIRST 300 FILES ARE HEADER INFO
		ice_conc[x] = np.reshape(data, [rows, cols])
		
	#divide by 250 to express in concentration
	ice_conc = ice_conc/scale_factor
	#GREATER THAN 250 is mask/land etc

	if (mask==1):
		ice_conc = ma.masked_where(ice_conc>1., ice_conc)
	
	if (maxConc==1):
		ice_conc = ma.where(ice_conc>1.,0, ice_conc)

	if (lowerConc==1):
		ice_conc = ma.where(ice_conc<0.15,0, ice_conc)

	if (monthMean==1):
		ice_conc=ma.mean(ice_conc, axis=0)

	return ice_conc

def get_ERA_precip_days(m, dataPath, yearT, dayT, varStr='tp'):
	
	f1 = Dataset(dataPath+'REANALYSES/ERAI/ERAI_'+varStr+'_'+str(yearT)+'.nc', 'r')

	# Units given in m of freshwater in a 12 hour period. 
	# So to convert to kg/m2/s multiply by den
	#var=var*1000./(60.*60.*12.)

	lon = f1.variables['longitude'][:]


	lat = f1.variables['latitude'][:]
	#print lat[1]-lat[0]

	# just get data above 45N
	lowerLatidx=int(90-45./(lat[0]-lat[1]))
	#print lowerLatidx
	lat=lat[0:lowerLatidx]
	xpts, ypts=m(*np.meshgrid(lon, lat))

	numday=dayT

	#in units of m of water so times by 1000, the density of water, to express this as kg/m2
	# data is every 12-hours, so need to multiply numdays by 2, then also sum over the first two time intervals
	varT=f1.variables[varStr][numday*2:(numday*2)+2, 0:lowerLatidx, :].astype(np.float16)*1000.
	var=np.sum(varT, axis=0)

	return xpts, ypts, lon, lat, var


def get_ERA_wind_days(m, dataPath, yearT, dayT):

	print(dataPath+'REANALYSES/ERAI/ERAI_'+'winds'+'_'+str(yearT)+'.nc')
	f1 = Dataset(dataPath+'REANALYSES/ERAI/ERAI_'+'winds'+'_'+str(yearT)+'.nc', 'r')


	lon = f1.variables['longitude'][:]

	#0:60 as just want to extrat the higher latitiudes. Repeat this in the variable extraction too.
	lat = f1.variables['latitude'][0:60]
	xpts, ypts=m(*np.meshgrid(lon, lat))


	#numYears = yearT-2009
	numday= dayT
	print(numday)

	# data is every 6-hours, so need to multiply numdays by 4, then also sum over the first four time intervals
	# only pick out the first 60 rows as want Arctic data only
	u10=f1.variables['u10'][numday*4:(numday*4)+4, 0:60, :].astype(np.float16)
	v10=f1.variables['v10'][numday*4:(numday*4)+4, 0:60, :].astype(np.float16)
	mag=np.mean(np.sqrt((u10**2)+(v10**2)), axis=0)

	return xpts, ypts, lon, lat, mag

def get_ERA5_precip_days_pyproj(proj, era5_data_path, yearStr, monStr, numday, lowerlatlim=0, varStr='sf'):

	print(yearStr, monStr, numday)

	f1 = Dataset(era5_data_path+'/ERA5_'+varStr+'_'+yearStr+monStr+'cds.nc', 'r')

	# Units given in m of freshwater in the previous 1 hour period. 
	# So to convert to kg/m2/s multiply by den
	#var=var*1000./(60.*60.*12.)

	lon = f1.variables['longitude'][:]
	
	lat = f1.variables['latitude'][:]
	
	#print lat[1]-lat[0]
	lowerLatidx=int((90-lowerlatlim)/(lat[0]-lat[1]))
	#print lowerLatidx
	lat=lat[0:lowerLatidx]

	# Had to switch lat and lon around when using proj!
	xpts, ypts=proj(*np.meshgrid(lon, lat))

	# in units of m of water so times by 1000, the density of water, to express this as kg/m2
	# data is every 1 hour (starting at 1 am on the first day of each month), so sum over the first 24 time intervals 
	#(index 0=1am to 23=midnight)

	varT=f1.variables[varStr][(numday*24):(numday*24)+24, 0:lowerLatidx, :].astype(np.float16)*1000.
	
	var=np.sum(varT, axis=0)

	return xpts, ypts, lon, lat, var

def get_ERA5_wind_days_pyproj(proj, era5_data_path, yearStr, monStr, numday, freq=6, lowerlatlim=0):

	print(yearStr, monStr, numday)

	f1 = Dataset(era5_data_path+'/ERA5_winds_'+yearStr+monStr+'cds.nc', 'r')

	lon = f1.variables['longitude'][:]
	
	lat = f1.variables['latitude'][:]
	
	#print lat[1]-lat[0]
	lowerLatidx=int((90-lowerlatlim)/(lat[0]-lat[1]))
	#print lowerLatidx
	lat=lat[0:lowerLatidx]

	# Had to switch lat and lon around when using proj!
	xpts, ypts=proj(*np.meshgrid(lon, lat))

	# data is every 1 hour (starting at 1 am on the first day of each month), so sum over the first 24 time intervals 
	#(index 0=1am to 23=midnight)
	# data is every hour

	u10=f1.variables['u10'][(numday*24):(numday*24)+24:freq, 0:lowerLatidx, :].astype(np.float16)
	v10=f1.variables['v10'][(numday*24):(numday*24)+24:freq, 0:lowerLatidx, :].astype(np.float16)
	mag=np.mean(np.sqrt((u10**2)+(v10**2)), axis=0)


	return xpts, ypts, lon, lat, mag


def create_grid(epsg_string='3413', dxRes=50000, lllat=36, llon=-90, urlat=36, urlon=90):
	""" Use pyproj to generate a grid covering the given domain (defined by the projection and the corner lat/lons)"""

	crs = pyproj.CRS.from_string("epsg:"+epsg_string)
	p=pyproj.Proj(crs)
	llcrn=p(llon, lllat)
	urcrn=p(urlon, urlat)

	print(llcrn)
	print(urcrn)

	nx = int((urcrn[0]-llcrn[0])/dxRes)+1
	ny = int((urcrn[1]-llcrn[1])/dxRes)+1
	print(nx, ny)

	x = llcrn[0]+dxRes*np.indices((ny,nx),np.float32)[1] # 1=column indices
	y = llcrn[1]+dxRes*np.indices((ny,nx),np.float32)[0] # 0=row indices

	lons, lats = p(x, y, inverse=True)
	
	return x, y, lats, lons, p

def get_ERA5_meltduration(m, dataPath, yearT):


	glob(dataPath+'REANALYSES/ERA5/ERA5_temp6hour'+'_'+str(yearT)+'*cds.nc')

	numDaysYearT=np.size(f1.variables['time'][:])/4
	print(numDaysYearT)

	lon = f1.variables['longitude'][:]
	lat = f1.variables['latitude'][:]
	lowerLatidx=int(45./(lat[0]-lat[1]))
	print(lowerLatidx)
	lat=lat[0:lowerLatidx]

	xpts, ypts=m(*np.meshgrid(lon, lat))


	# data is every 6-hours, so need to multiply numdays by 4, then also sum over the first four time intervals
	# only pick out the first 60 rows as want Arctic data only
	temp2mAnnual=ma.masked_all((numDaysYearT, lowerLatidx, lon.shape[0]))
	for x in range(numDaysYearT):
		temp2mAnnual[x]=np.mean(f1.variables['t2m'][x*4:(x*4)+4, 0:lowerLatidx, :].astype(float16), axis=0)-273.15


	return xpts, ypts, lon, lat, temp2mAnnual

def getCDRconcproj(proj, fileT, mask=1, maxConc=0, lowerConc=0):
	"""
	Grab the CDR ice concentration

	"""

	f = Dataset(fileT, 'r')
	print(fileT)

	# read lat/lon (start points) and lat1/lon1 (end points)
	lon = (f.variables['longitude'][:])
	lat = (f.variables['latitude'][:])
	# Convert percent to conc!
	conc = (f.variables['seaice_conc_cdr'][0])
	f.close()

	if (mask==1):
		conc = ma.masked_where(conc>1., conc)
	
	if (maxConc==1):
		conc = ma.where(conc>1.,0, conc)

	if (lowerConc==1):
		conc = ma.where(conc<0.15,0, conc)

	# transform to map project coordinates (basemap's axes, assume they have unit of m)
	x0, y0=proj(lon, lat)
	
	return conc, lat, lon, x0, y0


def read_icebridge_snowdepths(proj, dataPath, year, mask=1):
	"""Script to read in (quicklook and final) icebridge sea ice data"""

	xpts_total=[] 
	ypts_total=[]
	snow_thickness_days=[]
	dates=[]
	if (year>2012):
		files = glob(dataPath+'/quicklook/*'+str(year)+'*/*.txt')
	else:
		files = glob(dataPath+'/final/*'+str(year)+'*/*.txt')
	
	for x in range(np.size(files)):
		data = np.genfromtxt(files[x], delimiter=',', skip_header=1, dtype=str)
		# data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
		lats = data[:, 0].astype(float)
		lons = data[:, 1].astype(float)
		snow_thickness = data[:, 7].astype(float)
		xpts,ypts = proj(lons, lats)

		date = str(year)+files[x].split('/')[-1].split(str(year))[-1][0:4]
		print(date)

		if (mask==1):
			good_data=np.where((snow_thickness>=0.)&(snow_thickness<=2.))
			snow_thickness = snow_thickness[good_data]
			xpts = xpts[good_data]
			ypts = ypts[good_data]

		xpts_total.append(xpts)
		ypts_total.append(ypts)
		snow_thickness_days.append(snow_thickness)
		dates.append(date)

	return xpts_total, ypts_total, dates, snow_thickness_days
	
def get_coastal_mask(path, proj):
	
	coast_file = data=xr.open_dataset(path+'distance-to-coast_2m.nc')
	lat_c = coast_file['lat'][:]
	lon_c = coast_file['lon'][:]
	z_c = coast_file['z'][:]

	#REMOVE SOME OF THE COAST DISTANCE DATA THAT ISN'T NEEDED
	lat_c=lat_c[4250:-1]
	lon_c=lon_c[::20]
	z_c=z_c[4250:-1, ::20]
	print(np.amin(lat_c), np.amax(lat_c))
	print(np.amin(lon_c), np.amax(lon_c))

	xpts, ypts = proj(lon_c, lat_c)

	return xpts, ypts, z_c

def plot_gridded_cartopy(lons, lats, var, proj=ccrs.NorthPolarStereo(central_longitude=-45), shading='flat', out='./figure', units_lab='units', varStr='',
 minval=1., maxval=1., date_string='', month_string='', extra='',cbar_type='both', cmap_1=plt.cm.viridis, norm=0):

	#proj = ccrs.epsg(epsg_string)
	#proj = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0, central_latitude=90, false_easting=0.0, false_northing=0.0, globe=None)
	
	# The projection keyword determines how the plot will look
	fig=plt.figure(figsize=(5, 6))
	ax = plt.axes(projection = proj)
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	cs=ax.pcolormesh(lons, lats, var, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)
	ax.coastlines(zorder=3)
	ax.gridlines(draw_labels=True,
              linewidth=0.22, color='gray', alpha=0.5, linestyle='--')
	
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	# for some reason this extent can freak out if you set 180 to 180
	ax.set_extent([-179, 179, 45, 90], ccrs.PlateCarree())
	cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
	cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
	cb.set_label(varStr+' ('+units_lab+')',size=8)
	ax.set_title(varStr+' '+date_string+month_string+extra)

	plt.tight_layout()
	#print 'Saving..'
	plt.savefig(out+'.png', dpi=200)
	plt.close(fig)

def get_uv_from_xy(xdrift, ydrift, lon):
	# convert the x/y drift vectors to zonal/meridional vectors
	alpha = lon*np.pi/180.
	uvelT = ydrift*np.sin(alpha) + xdrift*np.cos(alpha)
	vvelT = ydrift*np.cos(alpha) - xdrift*np.sin(alpha) 

	return uvelT, vvelT

def get_nsidc_driftv4(nsidc_raw_path, year, day, time_freq):

	fileT=nsidc_raw_path+time_freq+'/icemotion_'+time_freq+'_nh_25km_'+str(year)+'0101_'+str(year)+'1231_v4.1.nc'
	nsidc_drifts_daily= xr.open_dataset(fileT)

	# Read in day, convert cm per second to meters per second
	xeasedrift = nsidc_drifts_daily.u[day].values*0.01
	yeasedrift = nsidc_drifts_daily.v[day].values*0.01

	lont=nsidc_drifts_daily.longitude.values
	latt=nsidc_drifts_daily.latitude.values

	#uvelT, vvelT = get_uv_from_xy(xeasedrift, yeasedrift, lont)
	
	return xeasedrift, yeasedrift, lont, latt

def plot_drift_cartopy(lons, lats, xpts, ypts, var_x, var_y, var_mag, proj=ccrs.NorthPolarStereo(central_longitude=-45), shading='flat', out='./figure', units_lab='units', units_vec='', varStr='',
 minval=1., maxval=1., date_string='year', month_string='months', extra='', res=2, scale_vec=1, vector_val=1, cbar_type='both', cmap_1=plt.cm.viridis, norm=0):

	#proj = ccrs.epsg(epsg_string)
	#proj = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0, central_latitude=90, false_easting=0.0, false_northing=0.0, globe=None)
	
	# The projection keyword determines how the plot will look
	fig=plt.figure(figsize=(5, 6))
	ax = plt.axes(projection = proj)
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	cs = ax.pcolormesh(lons, lats, var_mag, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)
	
	Q = ax.quiver(xpts[::res, ::res], ypts[::res, ::res], var_x[::res, ::res], var_y[::res, ::res], units='inches',scale=scale_vec, zorder=5)
	
	ax.coastlines(zorder=3)
	ax.gridlines(draw_labels=True,
              linewidth=0.22, color='gray', alpha=0.5, linestyle='--')
	
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	# for some reason this extent can freak out if you set 180 to 180
	ax.set_extent([-179, 179, 50, 90], ccrs.PlateCarree())
	cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
	cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
	cb.set_label(varStr+units_lab,size=8)
	ax.set_title(' '+varStr+' '+date_string+month_string+extra)

	qk = plt.quiverkey(Q, 3389969, 3389969, vector_val, str(vector_val)+' '+units_vec, coordinates='data', zorder = 11)   


	plt.tight_layout()
	#print 'Saving..'
	plt.savefig(out+'.png', dpi=200)
	plt.close(fig)

def plot_drift_cartopy_uv(lons, lats, var_u, var_v, var_mag, proj=ccrs.NorthPolarStereo(central_longitude=-45), shading='flat', out='./figure', units_lab='units', units_vec='', varStr='',
 minval=1., maxval=1., date_string='year', month_string='months', extra='', res=2, scale_vec=1, vector_val=1, cbar_type='both', cmap_1=plt.cm.viridis, norm=0):

	#proj = ccrs.epsg(epsg_string)
	#proj = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0, central_latitude=90, false_easting=0.0, false_northing=0.0, globe=None)
	
	# The projection keyword determines how the plot will look
	fig=plt.figure(figsize=(5, 6))
	ax = plt.axes(projection = proj)
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	cs = ax.pcolormesh(lons, lats, var_mag, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

	# Following https://github.com/SciTools/cartopy/issues/1179 we need to rescale the 
	u_src_crs = var_u / np.cos(lats / 180 * np.pi)
	v_src_crs = var_v
	magnitude = ma.sqrt(var_u**2 + var_v**2)
	magn_src_crs = ma.sqrt(u_src_crs**2 + v_src_crs**2)

	var_u_scaled = u_src_crs * magnitude / magn_src_crs
	var_v_scaled = v_src_crs * magnitude / magn_src_crs

	Q = ax.quiver(lons[::res, ::res], lats[::res, ::res], var_u_scaled[::res, ::res], var_v_scaled[::res, ::res], transform=ccrs.PlateCarree(), units='inches',scale=scale_vec, zorder=5)
	
	ax.coastlines(zorder=3)
	ax.gridlines(draw_labels=True,
              linewidth=0.22, color='gray', alpha=0.5, linestyle='--')
	
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	# for some reason this extent can freak out if you set 180 to 180
	ax.set_extent([-179, 179, 50, 90], ccrs.PlateCarree())
	cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
	cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
	cb.set_label(varStr+units_lab,size=8)
	ax.set_title(' '+varStr+' '+date_string+month_string+extra)

	qk = plt.quiverkey(Q, 3389969, 3389969, vector_val, str(vector_val)+' '+units_vec, coordinates='data', zorder = 11)   


	plt.tight_layout()
	#print 'Saving..'
	plt.savefig(out+'.png', dpi=200)
	plt.close(fig)

def get_psnlatslons(data_path):
	""" Get NSIDC 25 km lat/lon grids"""

	mask_latf = open(data_path+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(data_path+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	return lats_mask, lons_mask

def get_pmask(year, month):
	#remove half a degree as gridding around the pole hole edge
	
	if (year<1987):
		pmask=84.4
	elif((year==1987)&(month<=6)):
		pmask=84.4
	elif ((year==1987)&(month>6)):
		pmask=86.7
	elif ((year>1987)&(year<2008)):
		pmask=87.2
	else:
		pmask=89.2
	
	return pmask

def getGrid(outPath, dx):
	"""Get model grid data"""

	dxStr=str(int(dx/1000))+'km'
	gridData=load(outPath+'gridData'+dxStr+'.txt')
	lonG=gridData[0]
	latG=gridData[1]
	xptsG=gridData[2]
	yptsG=gridData[3]

	return lonG, latG, xptsG, yptsG


def getSTOSIWIGday(m, dayFiles, delim, mask_hs=1):
	"""  Get all snow radar data from all files in one OIB campaign day
	""" 
	convFactor=0.7809 # Speed of light to speed of light through snow
	lats_total=[] 
	lons_total=[]
	xpts_total=[]
	ypts_total=[]
	snow_thickness_total=[]

	for fileT in dayFiles:
		#print 'File:', fileT
		print (fileT)
		data = np.genfromtxt(fileT, delimiter=delim, skip_header=0, dtype=float)
		# data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
		lats = data[:, 1].astype(float)
		lons = data[:, 2].astype(float)
		snowRange = np.round(data[:, 3].astype(float), decimals=3)
		snowDepth= snowRange*convFactor

		if (mask_hs==1):
			goodhs=np.where((snowDepth>=0.)&(snowDepth<=1.5))
			lats = np.array(lats[goodhs])
			lons = np.array(lons[goodhs])

			snowDepth = np.array(snowDepth[goodhs])

		lats_total.extend(lats)
		lons_total.extend(lons)

		snow_thickness_total.extend(snowDepth)

		xpts, ypts = m(lons, lats)
		xpts_total.extend(xpts)
		ypts_total.extend(ypts)

	return xpts_total, ypts_total, lats_total, lons_total, snow_thickness_total


def getSTOSIWIGyear_proj(proj, dataPath, snowTypeT, yearT):
	"""  Get all snow radar data from all days within a campaign year.

	 Calls getSTOSIWIGday
	""" 

	if (snowTypeT=='GSFC'):
		delim='\t'
		endStr='txt'
	elif (snowTypeT=='JPL'):
		delim=','
		endStr='JPL'
	elif (snowTypeT=='SRLD'):
		delim=','
		endStr='srld'

	print(snowTypeT, yearT)
	folders = glob(dataPath+'/'+snowTypeT+'/'+str(yearT)+'*')
	print ('folders', folders)
	datesY=[folder[-8:] for folder in folders]


	latsY=[] 
	lonsY=[]
	xptsY=[]
	yptsY=[]
	snowY=[]

	for folder in folders:
		dayFilesT=glob(folder+'/*.'+endStr)
		#for day in folder 

		xptsD, yptsD, latsD, lonsD, snowD = getSTOSIWIGday(proj, dayFilesT, delim)
		xptsY.append(xptsD)
		yptsY.append(yptsD)
		latsY.append(latsD)
		lonsY.append(lonsD)
		snowY.append(snowD)

	return xptsY, yptsY, latsY, lonsY, snowY, datesY

def getWarren(lonT, latT, monthT):
	H_0 = [28.01, 30.28, 33.89, 36.8, 36.93, 36.59, 11.02, 4.64, 15.81, 22.66, 25.57, 26.67]
	a = [.127, .1056, .5486, .4046, .0214, .7021, .3008, .31, .2119, .3594, .1496, -0.1876]
	b = [-1.1833, -0.5908, -0.1996, -0.4005, -1.1795, -1.4819, -1.2591, -0.635, -1.0292, -1.3483, -1.4643, -1.4229]
	c = [-0.1164, -0.0263, 0.0280, 0.0256, -0.1076, -0.1195, -0.0811, -0.0655, -0.0868, -0.1063, -0.1409, -0.1413]
	d = [-0.0051, -0.0049, 0.0216, 0.0024, -0.0244, -0.0009, -0.0043, 0.0059, -0.0177, 0.0051, -0.0079, -0.0316]
	e = [0.0243, 0.0044, -0.0176, -0.0641, -0.0142, -0.0603, -0.0959, -0.0005, -0.0723, -0.0577, -0.0258, -0.0029]

	#Convert lat and lon into degrees of arc, +x axis along 0 degrees longitude and +y axis along 90E longitude
	x = (90.0 - latT)*np.cos(lonT * np.pi/180.0)	

	# To convert lat lon values into proper x and y coordinates
	y = (90.0 - latT)*np.sin(lonT * np.pi/180.0) 

	Hsw = H_0[monthT] + a[monthT]*x + b[monthT]*y + c[monthT]*x*y + (d[monthT]*x*x) + (e[monthT]*y*y)
	#Hsw=Hsw/100.

	#Hsw[where(region_maskG<9.6)]=np.nan
	#Hsw[where(region_maskG==14)]=np.nan
	#Hsw[where(region_maskG>15.5)]=np.nan

	return np.array(Hsw)


def get_budgets2layers_day(outStrings, outPath, folderStr, dayT, totalOutStr, region='', converttocm=0):

	data=xr.open_dataset(outPath+folderStr+'/budgets/'+totalOutStr+'.nc') 

	iceConcDay=np.array(data['iceConc'][dayT])
	#iceConcDay=ma.masked_where(iceConcDay<0.15, iceConcDay)


	snowBudget=[]
	for outString in outStrings:
		if (outString=='density'):
			snowDataT=data['density'][dayT]
			#snowDataT= ma.masked_where(iceConcDays<0.15, snowDataT)
			snowDataT= ma.masked_where(iceConcDay<0.15, snowDataT)
			# Need at least 2 cm for a density to be used (stop weird behaviour at very low snow depth)
			snowDepthT=data['snowDepth'][dayT, 0] + data['snowDepth'][dayT, 1]
			snowDataT= ma.masked_where(snowDepthT<0.02, snowDataT)

		elif (outString=='snowDepthNew'):
			snowDataT=data['snowDepth'][dayT, 0]
			#snowDataT= ma.masked_where(iceConcDays<0.15, snowDataT)
			snowDataT= ma.masked_where(iceConcDay<0.15, snowDataT)
		elif (outString=='snowDepthOld'):
			snowDataT=data['snowDepth'][dayT, 1]
			snowDataT= ma.masked_where(iceConcDay<0.15, snowDataT)
		elif (outString=='snowDepthTotal'):
			snowDataT=data['snowDepth'][dayT, 0] + data['snowDepth'][dayT, 1]
			snowDataT= ma.masked_where(iceConcDay<0.15, snowDataT)
		elif (outString=='snowDepthTotalConc'):
			snowDataT=(data['snowDepth'][dayT, 0] + data['snowDepth'][dayT, 1])/iceConcDay
			snowDataT= ma.masked_where(iceConcDay<0.15, snowDataT)
		elif (outString=='Precip'):
			# Meters of water equivalent (summed over the time period)
			if (dayT>0):
				precipT=data[outString][0:dayT] #.fillna(0)
			else:
				precipT=data[outString]
			snowDataT = np.sum(precipT/200., axis=0)

		else:

			snowDataT = data[outString][dayT]

		#print 'region len', len(region)
		if (len(region)>0):
			print ('masking region')
			regionM=load(outPath+region)
			snowDataT=ma.masked_where(regionM<0.5, snowDataT)

		snowDataT= ma.masked_where(np.isnan(snowDataT), snowDataT)

		if (converttocm==1):
			snowDataT=snowDataT*100.

		snowBudget.append(snowDataT)

	if (np.size(outStrings)>1):
		return snowBudget
	else:
		print ('1 var')
		return snowDataT

def densityClim(dayT, dataPath):
	"""Assign initial snow density based on daily Warren climatology"""

	densityClim=pd.read_csv(dataPath+'/W99_density.csv', header=0, names=['Day', 'Density'])
	#find density of given day and multiply by 1000 to express in Kg
	densityClimDay=1000*densityClim['Density'].iloc[dayT-1]

	return densityClimDay

def get_region_mask_pyproj(anc_data_path, proj, xypts_return=0):
	""" Read in NSIDC Arctic Ocean mask and transform to prescribed projection
	"""

	header = 300
	datatype='uint8'
	file_mask = anc_data_path+'region_n.msk'

	#1 - Non regional ocean
	#2 - Sea of Okhotsk and Japan
	#3 - Bering 
	#4 - Hudson Bay
	#5 - Baffin Bay/Davis Strait/Labrador Sea
	#6 - Greenland
	#7 - Kara/Barents
	#8 - Arctic Ocean
	#9 - Canadian Archipelago
	#10 - Gulf of St Lawrence
	#11 - Land

	fd = open(file_mask, 'rb')
	region_mask = np.fromfile(file=fd, dtype=datatype)
	region_mask = np.reshape(region_mask[header:], [448, 304])

	if (xypts_return==1):
		mask_latf = open(anc_data_path+'/psn25lats_v3.dat', 'rb')
		mask_lonf = open(anc_data_path+'/psn25lons_v3.dat', 'rb')
		lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		# switched lat and lon from basemap
		xpts, ypts = proj(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask


