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
		05/10/2020: Version 2: Converted to Python 3
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



def int_smooth_drifts(xptsG, yptsG, xptsF, yptsF, latsF, driftFmon, sigma_factor=0.75):
	
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

	# COPY TO GENERATE AN ARRAY OF ZEROS INSTEAD OF NANS FOR FILTERING.
	driftFGxN=np.copy(driftFGx)
	driftFGyN=np.copy(driftFGy)

	driftFGx[np.isnan(driftFGx)]=0
	driftFGy[np.isnan(driftFGy)]=0

	driftFGxg = gaussian_filter(driftFGx, sigma=sigma_factor)
	driftFGyg = gaussian_filter(driftFGy, sigma=sigma_factor)

	driftFG[0]=ma.masked_where(np.isnan(driftFGxN), driftFGxg)
	driftFG[1]=ma.masked_where(np.isnan(driftFGyN), driftFGyg)
	

	return driftFG

def getOSISAFDrift(m, fileT):
	"""
	Calculate the OSI-SAF vectors on our map projection
	With help from Thomas Lavergne!

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

	# transform to map project coordinates (basemap's axes, assume they have unit of m)
	x0, y0=proj(lon, lat)
	x1, y1=proj(lon1, lat1)

	xpts=(x0+x1)/2.
	ypts=(y0+y1)/2.

	# normalize drift components to m/s (x1-x0 is already m, so we just divide by 2-days worth of seconds)
	xt=(x1-x0)/(60*60*24*2.)
	yt=(y1-y0)/(60*60*24*2.)

	# TL: no need to rotate : the xt, and yt are already in the basemap's projection

	# compute magnitude (speed scalar)
	mag=np.sqrt(xt**2+yt**2)

	return xt, yt, mag, lat, lon, xpts, ypts

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

	x = llcrn[0]+dxRes*np.indices((ny,nx),np.float32)[0] # 1=column indicdes
	y = llcrn[1]+dxRes*np.indices((ny,nx),np.float32)[1] # 0=row indices

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
	if (year>2013):
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

		date=files[x][-12:-4]

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
	ax.set_extent([-179, 179, 50, 90], ccrs.PlateCarree())
	cax,kw = mcbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
	cb=fig.colorbar(cs,cax=cax,extend='both',**kw)
	cb.set_label(varStr+' ('+units_lab+')',size=8)
	ax.set_title(varStr+' '+date_string+month_string+extra)

	plt.tight_layout()
	#print 'Saving..'
	plt.savefig(out+'.png', dpi=200)
	plt.close(fig)

def plot_drift_cartopy(lons, lats, xpts, ypts, var_u, var_v, var_mag, proj=ccrs.NorthPolarStereo(central_longitude=-45), shading='flat', out='./figure', units_lab='units', units_vec='', varStr='',
 minval=1., maxval=1., date_string='year', month_string='months', extra='', res=2, scale_vec=1, vector_val=1, cbar_type='both', cmap_1=plt.cm.viridis, norm=0):

	#proj = ccrs.epsg(epsg_string)
	#proj = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0, central_latitude=90, false_easting=0.0, false_northing=0.0, globe=None)
	
	# The projection keyword determines how the plot will look
	fig=plt.figure(figsize=(5, 6))
	ax = plt.axes(projection = proj)
	#ax.imshow(data, transform=ccrs.PlateCarree(), zorder=2)
	cs = ax.pcolormesh(lons, lats, var_mag, vmin=minval, vmax=maxval, transform=ccrs.PlateCarree(), zorder=2)

	Q = ax.quiver(xpts[::res, ::res], ypts[::res, ::res], var_u[::res, ::res], var_v[::res, ::res], units='inches',scale=scale_vec, zorder=5)
	
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

def get_region_maskAOsnow(datapath, mplot, xypts_return=0, latN=60):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'
	
	
	fd = open(file_mask, 'rb')
	region_mask = np.fromfile(file=fd, dtype=datatype)
	region_mask = np.reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskAO=np.zeros((lons_mask.shape))
	mask = np.where((region_mask<11) & (lats_mask>latN))
	region_maskAO[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskAO, xpts, ypts
	else:
		return region_maskAO

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
		data = genfromtxt(fileT, delimiter=delim, skip_header=0, dtype=float)
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


def getSTOSIWIGyear(m, dataPath, snowTypeT, yearT):
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
	folders = glob(dataPath+'/ICEBRIDGE/STOSIWIG/'+snowTypeT+'/'+str(yearT)+'*')
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

		xptsD, yptsD, latsD, lonsD, snowD = getSTOSIWIGday(m, dayFilesT, delim)
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

def bindataN(x, y, z, xG, yG, binsize=0.01, retbin=True, retloc=True):
	"""
	Place unevenly spaced 2D data on a grid by 2D binning (nearest
	neighbor interpolation).

	Parameters
	----------
	x : ndarray (1D)
		The idependent data x-axis of the grid.
	y : ndarray (1D)
		The idependent data y-axis of the grid.
	z : ndarray (1D)
		The dependent data in the form z = f(x,y).
	binsize : scalar, optional
		The full width and height of each bin on the grid.  If each
		bin is a cube, then this is the x and y dimension.  This is
		the step in both directions, x and y. Defaults to 0.01.
	retbin : boolean, optional
		Function returns `bins` variable (see below for description)
		if set to True.  Defaults to True.
	retloc : boolean, optional
		Function returns `wherebins` variable (see below for description)
		if set to True.  Defaults to True.

	Returns
	-------
	grid : ndarray (2D)
		The evenly gridded data.  The value of each cell is the median
		value of the contents of the bin.
	bins : ndarray (2D)
		A grid the same shape as `grid`, except the value of each cell
		is the number of points in that bin.  Returns only if
		`retbin` is set to True.
	wherebin : list (2D)
		A 2D list the same shape as `grid` and `bins` where each cell
		contains the indicies of `z` which contain the values stored
		in the particular bin.

	Revisions
	---------
	2010-07-11  ccampo  Initial version
	"""
	# get extrema values.
	xmin, xmax = xG.min(), xG.max()
	ymin, ymax = yG.min(), yG.max()

	# make coordinate arrays.
	xi      = xG[0]
	yi      = yG[:, 0] #np.arange(ymin, ymax+binsize, binsize)
	xi, yi = np.meshgrid(xi,yi)

	# make the grid.
	grid           = np.zeros(xi.shape, dtype=x.dtype)
	nrow, ncol = grid.shape
	if retbin: bins = np.copy(grid)

	# create list in same shape as grid to store indices
	if retloc:
		wherebin = np.copy(grid)
		wherebin = wherebin.tolist()

	# fill in the grid.
	for row in range(nrow):
		for col in range(ncol):
			xc = xi[row, col]    # x coordinate.
			yc = yi[row, col]    # y coordinate.

			# find the position that xc and yc correspond to.
			posx = np.abs(x - xc)
			posy = np.abs(y - yc)
			ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
			ind  = np.where(ibin == True)[0]

			# fill the bin.
			bin = z[ibin]
			if retloc: wherebin[row][col] = ind
			if retbin: bins[row, col] = bin.size
			if bin.size != 0:
				binval         = np.mean(bin)
				grid[row, col] = binval
			else:
				grid[row, col] = np.nan   # fill empty bins with nans.

	# return the grid
	if retbin:
		if retloc:
			return grid, bins, wherebin
		else:
			return grid, bins
	else:
		if retloc:
			return grid, wherebin
		else:
			return grid				

def correlateVars(var1, var2):
#correlate two variables
	trend, intercept, r_a, prob, stderr = stats.linregress(var1, var2)
	sig = 100.*(1.-prob)
	return trend, sig, r_a, intercept 

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

def get_region_mask(datapath, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'

	#8 - Arctic Ocean
	#9 - Canadian Archipelago
	#10 - Gulf of St Lawrence
	#11 - Land

	fd = open(file_mask, 'rb')
	region_mask = np.fromfile(file=fd, dtype=datatype)
	region_mask = np.reshape(region_mask[header:], [448, 304])

	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
		lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask

def densityClim(dayT):
	"""Assign initial snow density based on daily climatology"""

	densityClim=pd.read_csv(dataPath+'/Daily_Density.csv', header=0, names=['Day', 'Density'])
	#find density of given day and multiply by 1000 to express in Kg
	densityClimDay=1000*densityClim['Density'].iloc[dayT-1]

	return densityClimDay

def get_region_mask_pyproj(anc_data_path, proj, xypts_return=0):
	""" Read in NSIDC Arctic Ocean mask and transofrm to given projection
	"""

	header = 300
	datatype='uint8'
	file_mask = anc_data_path+'region_n.msk'

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

