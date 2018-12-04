""" commonFuncs.py
	
	Common functions used by the NESOSIM.py script 
	Model written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov) or refer to the GitHub site (ADD THIS)


	Python dependencies:
		See below for the relevant module imports. Of note:
		matplotlib
		basemap

	Update history:
		03/01/2018: Version 1

"""


from pylab import *
from glob import glob
from scipy.interpolate import griddata
import xarray as xr
from scipy import stats

def defGrid(m, dxRes=50000):
    nx = int((m.xmax-m.xmin)/dxRes)+1; ny = int((m.ymax-m.ymin)/dxRes)+1
    gridStr=str(int(dxRes/1000))+'km'
    lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)

    return lonsG, latsG, xptsG, yptsG, nx, ny


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

def getCSatDrift(file, m, numDaysLag=3, rotatetoXY=1):
	
	f = Dataset(file, 'r')
	u = f.variables['zonal_motion'][0]/(60.*60.*24.*numDaysLag)
	v = f.variables['meridional_motion'][0]/(60.*60.*24.*numDaysLag)
	lon = f.variables['longitude'][:]
	lat = f.variables['latitude'][:]
    
    #less than 0 are flags meaning the drift hasnt passed certain tests. 
    #q = f.variables['quality_flag'][0]
    #u = ma.masked_where((q<=0), u)
    #v = ma.masked_where((q<=0), v)

    #ROTATE VECOTRS TO X/Y GRID
	if (rotatetoXY==1):
		ux,vy = m.rotate_vector(u,v,lon,lat)
		return ux, vy, lon, lat
	else:
		return u, v, lon, lat

def getKimuradriftDayRaw(rawdatapath, fileT, m):

	lonlatK = loadtxt(rawdatapath+'/ICE_DRIFT/KIMURA/'+'latlon_amsr_ads145.txt', unpack=True)
	latsK=flipud(lonlatK[2].reshape(145, 145))
	lonsK=flipud(lonlatK[3].reshape(145, 145))
	xptsK, yptsK=m(lonsK, latsK)
	alphaK = lonsK*pi/180.

	# Comes in xy coordinates so need to rotate to UV
	     
	KFile = open(fileT, 'r')
	Kdrift = fromfile(file=KFile, dtype=float32)
	Kdrift=ma.masked_where(Kdrift>900, Kdrift/100.)
	xvel=-flipud(Kdrift[0::2].reshape(145, 145))
	yvel=-flipud(Kdrift[1::2].reshape(145, 145))

	uvelK = yvel*sin(alphaK) + xvel*cos(alphaK)
	vvelK = yvel*cos(alphaK) - xvel*sin(alphaK) 

	xvelG,yvelG = m.rotate_vector(uvelK,vvelK,lonsK,latsK)

	xvelG[where(ma.getmask(xvelG))]=np.nan
	yvelG[where(ma.getmask(yvelG))]=np.nan
	driftKday=stack((xvelG, yvelG))


	return xptsK, yptsK, driftKday, lonsK, latsK

def getFowlerDrift(file,lon):
	
	fd = open(file, 'rb')
	motionDat = fromfile(file=fd, dtype='<i2')
	motionDat = reshape(motionDat, [361, 361, 3])

	xt = motionDat[:, :, 0]/1000.
	yt = motionDat[:, :, 1]/1000.         
	q = motionDat[:, :, 2]/1000.

	mask = where((q<=0) | (q>1), 0, 1)

	xt = ma.masked_where(mask<0.5, xt)
	yt = ma.masked_where(mask<0.5, yt)


	alpha = lon*pi/180.
	uvelT = yt*sin(alpha) + xt*cos(alpha)
	vvelT = yt*cos(alpha) - xt*sin(alpha) 

	return uvelT, vvelT

def smoothDriftDaily(xptsG, yptsG, xptsF, yptsF, latsF, driftFmon, sigmaG=0.75):
	from scipy.ndimage.filters import gaussian_filter
	
	nx=xptsG.shape[0]
	ny=xptsG.shape[1]

	driftFG=ma.masked_all((2, nx, ny))
	#driftFGK=ma.masked_all((2, nx, ny))


	badData = np.zeros((xptsF.shape))
	badData[where(np.isnan(driftFmon[1]) & (latsF>88))]=1

	driftFx = driftFmon[0][where(badData<0.5)]
	driftFy = driftFmon[1][where(badData<0.5)]
	xptsFM = xptsF[where(badData<0.5)]

	yptsFM = yptsF[where(badData<0.5)]

	driftFGx = griddata((xptsFM, yptsFM), driftFx, (xptsG, yptsG), method='linear')
	driftFGy = griddata((xptsFM, yptsFM), driftFy, (xptsG, yptsG), method='linear')

	# COPY TO GENERATE AN ARRAY OF ZEROS INSTEAD OF NANS FOR FILTERING.
	driftFGxN=np.copy(driftFGx)
	driftFGyN=np.copy(driftFGy)

	driftFGx[np.isnan(driftFGx)]=0
	driftFGy[np.isnan(driftFGy)]=0

	driftFGxg = gaussian_filter(driftFGx, sigma=sigmaG)
	driftFGyg = gaussian_filter(driftFGy, sigma=sigmaG)

	driftFG[0]=ma.masked_where(np.isnan(driftFGxN), driftFGxg)
	driftFG[1]=ma.masked_where(np.isnan(driftFGyN), driftFGyg)
	return driftFG

def getOSISAFDrift(m, fileT):
	"""
	Calculate the OSI-SAF vectors on our map projection
	With help from Thomas Lavergne!

	"""

	f = Dataset(fileT, 'r')
	print fileT

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
	mag=sqrt(xt**2+yt**2)
        #print mag.mean(), mag.min(), mag.max()

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
	print ('Num of days:', size(files))
	driftFmon=ma.masked_all((size(files), 2, lonsF.shape[0], lonsF.shape[1]))

	
	x=0
	for file in files:
	     
		uvelT, vvelT=getFowlerDrift(file,lonsF)

		xvel,yvel = m.rotate_vector(uvelT,vvelT,lonsF,latsF)

		xvel[where(ma.getmask(xvel))]=np.nan
		yvel[where(ma.getmask(yvel))]=np.nan
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
    lonsF = np.reshape(lonlatF[:, 3], (361, 361))
    latsF = np.reshape(lonlatF[:, 2], (361, 361))
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
	if (size(files)>0):
		print 'Same day conc file exists:'

	if (size(files)==0):
		# first try day before
		day_str = '%02d' % (day)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+day_str+'*')
		if (size(files)>0):
			print 'Using day before file:'

	# If still nothing try day after
	if (size(files)==0):
		# try day after
		day_str = '%02d' % (day+2)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+day_str+'*')
		if (size(files)>0):
			print 'Using day after file:'

	
	fd = open(files[0], 'r')
	data = fromfile(file=fd, dtype=datatype)
	data = data[header:]
	#FIRST 300 FILES ARE HEADER INFO
	ice_conc = reshape(data, [rows, cols])
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
	

	print 'Num conc files:', size(files), 'in month:'+month_str
	ice_conc = ma.masked_all((size(files), rows, cols))
	
	for x in xrange(size(files)):
		fd = open(files[x], 'r')
		data = fromfile(file=fd, dtype=datatype)
		data = data[header:]
		#FIRST 300 FILES ARE HEADER INFO
		ice_conc[x] = reshape(data, [rows, cols])
		
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

	leapYrs=[1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]
	daysInYr=[365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366]
	#f1 = Dataset(dataPath+'REANALYSES/ERAI/EraInterim_Snow_2002_to_2016.nc', 'r')

	#f1 = Dataset(dataPath+'REANALYSES/ERAI/EraInterim_Snow_2002_to_2016.nc', 'r')
	
	f1 = Dataset(dataPath+'REANALYSES/ERAI/ERAI_'+varStr+'_'+str(yearT)+'.nc', 'r')

	# Units given in m of freshwater in a 12 hour period. 
	# So to convert to kg/m2/s multiply by den
	#var=var*1000./(60.*60.*12.)

	lon = f1.variables['longitude'][:]

	#0:60 as just want to extrat the higher latitiudes. Repeat this in the variable extraction too.
	lat = f1.variables['latitude'][:]
	#print lat[1]-lat[0]
	lowerLatidx=int(45./(lat[0]-lat[1]))
	#print lowerLatidx
	lat=lat[0:lowerLatidx]
	xpts, ypts=m(*np.meshgrid(lon, lat))

	
	numday=dayT

	#in units of m of water so times by 1000, the density of water, to express this as kg/m2
	# data is every 12-hours, so need to multiply numdays by 2, then also sum over the first two time intervals
	varT=f1.variables[varStr][numday*2:(numday*2)+2, 0:lowerLatidx, :].astype(float16)*1000.
	var=sum(varT, axis=0)

	return xpts, ypts, lon, lat, var


def get_ERA_wind_days(m, dataPath, yearT, dayT):

	leapYrs=[1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]
	daysInYr=[365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366]
	#f1 = Dataset(dataPath+'REANALYSES/ERAI/EraInterim_Snow_2002_to_2016.nc', 'r')

	print dataPath+'REANALYSES/ERAI/ERAI_'+'winds'+'_'+str(yearT)+'.nc'
	f1 = Dataset(dataPath+'REANALYSES/ERAI/ERAI_'+'winds'+'_'+str(yearT)+'.nc', 'r')


	lon = f1.variables['longitude'][:]

	#0:60 as just want to extrat the higher latitiudes. Repeat this in the variable extraction too.
	lat = f1.variables['latitude'][0:60]
	xpts, ypts=m(*np.meshgrid(lon, lat))


	#numYears = yearT-2009
	numday= dayT
	print numday

	# data is every 6-hours, so need to multiply numdays by 4, then also sum over the first four time intervals
	# only pick out the first 60 rows as want Arctic data only
	u10=f1.variables['u10'][numday*4:(numday*4)+4, 0:60, :].astype(float16)
	v10=f1.variables['v10'][numday*4:(numday*4)+4, 0:60, :].astype(float16)
	mag=mean(sqrt((u10**2)+(v10**2)), axis=0)

	return xpts, ypts, lon, lat, mag


def plotSnow(m, xpts , ypts, var_mag, shading='flat', out='./figure', units_lab='units', units_vec=r'm s$^{-1}$',
 minval=1., maxval=1., base_mask=0,res=1, scale_vec=1, vector_val=1, date_string='year', month_string='months', extra='',cbar_type='both', cmap_1=plt.cm.RdBu_r, norm=0):

    #PLOT SCALAR FIELD WITH OVERLYING VECTORS. 
    #VAR MAG MAY NOT NECCESARRILY BE THE MAGNITUDE OF THE VECTORS (E.G. IN THE CASE OF WIND CURL)
	rcParams['ytick.major.size'] = 2
	rcParams['axes.linewidth'] = .25
	rcParams['lines.linewidth'] = .25
	rcParams['patch.linewidth'] = .25
	rcParams['ytick.labelsize']=8
	rcParams['legend.fontsize']=8
	rcParams['font.size']=8 
	rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

	#print 'Plotting'
	fig = figure(figsize=(4.,4.))
	ax1=gca()

	cmap=cmap_1
	#cmap.set_under('w')
	#var_mag=ma.masked_where(var_mag<1e8, var_mag)
	if (norm==1):
		#print 'Norm:', norm
		norm=MidPointNorm_Good(midpoint=0)
		im1 = m.pcolormesh(xpts , ypts, var_mag, norm=norm, cmap=cmap,vmin=minval, vmax=maxval,shading=shading, edgecolors='None', zorder=4, rasterized=True)
	else:
		im1 = m.pcolormesh(xpts , ypts, var_mag, cmap=cmap,vmin=minval, vmax=maxval,shading=shading, edgecolors='None', zorder=4, rasterized=True)
	
	# LOWER THE SCALE THE LARGER THE ARROW
	#im1c = m.contour(xptsC , yptsC, conc_mean, levels=[0.15], zorder=5)


	m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=10)
	#m.drawmapboundary(fill_color='0.3')
	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	if base_mask==1:
	#m.drawmapboundary(fill_color='0.4' , zorder=1)
	    m.fillcontinents(color='0.7',lake_color='grey', zorder=5)
	
	#m.drawcoastlines(linewidth=0.25, zorder=10)

	ax1.annotate(date_string, xy=(0.4, 0.93),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=10, zorder=10)


	#ADD COLORBAR TO MAP
	cax = fig.add_axes([0.72, 0.95, 0.22, 0.03])

	#cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend=cbar_type, use_gridspec=True)
	cbar.set_label(units_lab)
	cbar.set_ticks([minval, maxval])

   

	subplots_adjust(bottom=0.0, left=0.0, top = 0.99, right=1.0)
	#print 'Saving..'
	savefig(out+'.png', dpi=150)
	close(fig)
#plot vector map (with vectors in x/y directions)

def get_region_maskPSsnow(datapath, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'
	
	region_lonlat = [80, 180, 70, 82]
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskCA=np.zeros((lons_mask.shape))
	mask = where((lons_mask>region_lonlat[0]) & (lons_mask<region_lonlat[1]) & (lats_mask>region_lonlat[2]) & (lats_mask<region_lonlat[3])& (region_mask<11))
	region_maskCA[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskCA, xpts, ypts
	else:
		return region_maskCA


def get_region_maskNAsnow(datapath, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'
	
	region_lonlat = [-20, 55, 72, 85]
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskCA=np.zeros((lons_mask.shape))
	mask = where((lons_mask>region_lonlat[0]) & (lons_mask<region_lonlat[1]) & (lats_mask>region_lonlat[2]) & (lats_mask<region_lonlat[3])&(region_mask<11))
	region_maskCA[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskCA, xpts, ypts
	else:
		return region_maskCA

class MidPointNorm_Good(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if mpl.cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint
def get_psnlatslons(data_path):
	mask_latf = open(data_path+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(data_path+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	return lats_mask, lons_mask

def get_region_mask_sect(datapath, mplot, xypts_return=0):
	datatype='uint8'
	file_mask = datapath+'/OTHER/sect_fixed_n.msk'
	# 1   non-region oceans
	# ;           = 2   Sea of Okhotsk and Japan
	# ;           = 3   Bering Sea
	# ;           = 4   Hudson Bay
	# ;           = 5   Gulf of St. Lawrence
	# ;           = 6   Baffin Bay/Davis Strait/Labrador Sea
	# ;           = 7   Greenland Sea
	# ;           = 8   Barents Seas
	# ;           = 9   Kara
	# ;           =10   Laptev
	# ;           =11   E. Siberian
	# ;           =12   Chukchi
	# ;           =13   Beaufort
	# ;           =14   Canadian Archipelago
	# ;           =15   Arctic Ocean
	# ;           =20   Land
	# ;           =21   Coast
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask, [448, 304])

	#mask_latf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lats_v3.dat', 'rb')
	#mask_lonf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lons_v3.dat', 'rb')
	#lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	#lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	#xpts, ypts = mplot(lons_mask, lats_mask)
	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
		lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask


def get_day_concSN_daily(datapath, year, month, day, alg=0, pole='A', vStr='v03', mask=1, maxConc=0, lowerConc=0, monthMean=0):
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
	files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/daily/'+vStr+'/'+str(year)+'/'+team_s+'_'+str(year)+month_str+day_str+'*'+vStr+'*')
	if (size(files)>0):
		print ('Same day conc file exists:')

	if (size(files)==0):
		# first try day before
		day_str = '%02d' % (day)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/daily/'+vStr+'/'+str(year)+'/'+team_s+'_'+str(year)+month_str+day_str+'*'+vStr+'*')
		if (size(files)>0):
			print ('Using day before file:')

	# If still nothing try day after
	if (size(files)==0):
		# try day after
		day_str = '%02d' % (day+2)
		files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/daily/'+vStr+'/'+str(year)+'/'+team_s+'_'+str(year)+month_str+day_str+'*'+vStr+'*')
		if (size(files)>0):
			print ('Using day after file:')

	
	fd = open(files[0], 'r')
	data = fromfile(file=fd, dtype=datatype)
	data = data[header:]
	#FIRST 300 FILES ARE HEADER INFO
	ice_conc = reshape(data, [rows, cols])
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
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskCA=np.zeros((lons_mask.shape))
	mask = where((region_mask<11) & (lats_mask>latN))
	region_maskCA[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskCA, xpts, ypts
	else:
		return region_maskCA

def get_region_maskCAsnow(datapath, mplot, xypts_return=0):
	header = 300
	datatype='uint8'
	file_mask = datapath+'/OTHER/region_n.msk'
	
	region_lonlat = [-150, -20, 78, 88]
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	region_maskCA=np.zeros((lons_mask.shape))
	mask = where((lons_mask>region_lonlat[0]) & (lons_mask<region_lonlat[1]) & (lats_mask>region_lonlat[2]) & (lats_mask<region_lonlat[3])& (region_mask==8))
	region_maskCA[mask]=1

	if (xypts_return==1):

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_maskCA, xpts, ypts
	else:
		return region_maskCA

def getGrid(outPath, dx):
  """Get model grid data"""

  dxStr=str(int(dx/1000))+'km'
  gridData=load(outPath+'gridData'+dxStr+'.txt')
  lonG=gridData[0]
  latG=gridData[1]
  xptsG=gridData[2]
  yptsG=gridData[3]

  return lonG, latG, xptsG, yptsG

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
    print fileT
    data = genfromtxt(fileT, delimiter=delim, skip_header=0, dtype=float)
    # data is a table-like structure (a numpy recarray) in which you can access columns and rows easily
    lats = data[:, 1].astype(float)
    lons = data[:, 2].astype(float)
    snowRange = np.round(data[:, 3].astype(float), decimals=3)
    snowDepth= snowRange*convFactor

    if (mask_hs==1):
      goodhs=where((snowDepth>=0.)&(snowDepth<=1.5))
      lats = array(lats[goodhs])
      lons = array(lons[goodhs])
      
      snowDepth = array(snowDepth[goodhs])

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

  print snowTypeT, yearT
  folders = glob(dataPath+'/ICEBRIDGE/STOSIWIG/'+snowTypeT+'/'+str(yearT)+'*')
  print 'folders', folders
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
  
  iceConcDay=array(data['iceConc'][dayT])
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
      snowDataT = sum(precipT/200., axis=0)

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
  
  if (size(outStrings)>1):
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
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask[header:], [448, 304])

	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
		lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask

