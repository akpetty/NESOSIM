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

	print 'Plotting'
	fig = figure(figsize=(4.,4.))
	ax1=gca()

	cmap=cmap_1
	#cmap.set_under('w')
	#var_mag=ma.masked_where(var_mag<1e8, var_mag)
	if (norm==1):
		print 'Norm:', norm
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
	
	m.drawcoastlines(linewidth=0.25, zorder=10)

	ax1.annotate(date_string, xy=(0.04, 0.93),xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', fontsize=10, zorder=10)


	#ADD COLORBAR TO MAP
	cax = fig.add_axes([0.7, 0.95, 0.26, 0.03])

	#cax = fig.add_axes([0.1, 0.1, 0.8, 0.04])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend=cbar_type, use_gridspec=True)
	cbar.set_label(units_lab)
	cbar.set_ticks(np.linspace(minval, maxval, 4))

   

	subplots_adjust(bottom=0.0, left=0.0, top = 0.99, right=1.0)
	print 'Saving..'
	savefig(out+'.png', dpi=150)
	close(fig)
#plot vector map (with vectors in x/y directions)

