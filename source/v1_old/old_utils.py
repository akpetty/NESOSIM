
def correlateVars(var1, var2):
#correlate two variables
	trend, intercept, r_a, prob, stderr = stats.linregress(var1, var2)
	sig = 100.*(1.-prob)
	return trend, sig, r_a, intercept 

def get_nesosim_day(outPath, end_year, day, converttocm=1):

	files=glob(outPath+'*'+str(end_year)+'.nc')
	data=xr.open_dataset(files[0]) 
	snow_depth=data['snowDepth'][day].values
	lat=data['latitude'][:].values
	lon=data['longitude'][:].values

	if (converttocm==1):
		snow_depth=snow_depth*100.

	return snow_depth, lat, lon

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


def defGrid(m, dxRes=50000):
    nx = int((m.xmax-m.xmin)/dxRes)+1; ny = int((m.ymax-m.ymin)/dxRes)+1
    gridStr=str(int(dxRes/1000))+'km'
    lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)

    return lonsG, latsG, xptsG, yptsG, nx, ny

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
	region_mask = np.fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask, [448, 304])

	#mask_latf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lats_v3.dat', 'rb')
	#mask_lonf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lons_v3.dat', 'rb')
	#lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	#lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	#xpts, ypts = mplot(lons_mask, lats_mask)
	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
		lats_mask = reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask



class MidPointNorm_Good(plt.Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        plt.Normalize.__init__(self,vmin, vmax, clip)
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
