""" gridInitialConds.py
	
	Script to grid the temperature-scaled initiail (august) snow depths
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: Temp-scaled snow depths, provided by Melinda Webster
	Output: binned initiail snow depths

	Python dependencies:
		See below for the relevant module imports
		Also some function in commongFuncs.py

	Update history:
		10/01/2018: Version 1
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


outPath='/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/'
dataPath='../../Data/'
figpath='/Volumes/PETTY_PASSPORT3/NESOSIM/Figures/InitialConds/'

m = Basemap(projection='npstere',boundinglat=60,lon_0=-45, resolution='l', round=False)
dx=100000.
dxStr=str(int(dx/1000))+'km'
print dxStr
region_maskG=load(outPath+'regionMaskG'+dxStr)
lonG, latG, xptsG, yptsG, nx, ny= cF.defGrid(m, dxRes=dx)


def main(year):
	lons, lats, snow=loadtxt(dataPath+'InitialConditions/'+'Initial_'+str(year)+'.txt', delimiter=',', unpack=True, skiprows=1)
	xpts, ypts=m(lons, lats)
	snow=snow/100.


	snowG = griddata((xpts, ypts), snow, (xptsG, yptsG), method='linear')

	snowG[where(latG<70)]=0
	snowG[where(region_maskG>8.2)]=0
	snowG[where(region_maskG<=7.8)]=0
	snowG[where(snowG<0)]=0

	day=226
	dayStr=str(day) #226 is middle of August

	iceConcDayG=load(outPath+'/IceConc/'+'/'+str(year)+'/iceConcG'+dxStr+'-'+str(year)+'_d'+dayStr)
	snowG[where(iceConcDayG<0.15)]=0
	
	snowG.dump(outPath+'InitialConds/ICsnow'+str(year)+'-'+dxStr)

	rcParams['ytick.major.size'] = 2
	rcParams['axes.linewidth'] = .5
	rcParams['lines.linewidth'] = .5
	rcParams['patch.linewidth'] = .5
	rcParams['ytick.labelsize']=8
	rcParams['legend.fontsize']=8
	rcParams['font.size']=8
	rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


	fig = figure(figsize=(3, 3))
	ax1=gca()
	minval=0
	maxval=0.1
	cints=0.025
	im1 = m.contourf(xptsG , yptsG, ma.masked_where(snowG==0, snowG), levels=np.arange(minval, maxval+cints*0.1, cints), cmap=cm.Reds ,extend='both', shading='gouraud', edgecolors='None', zorder=4, rasterized=True)
		
	#im1 = m.pcolormesh(xptsG , yptsG, AugSnow, vmin=minval, vmax=maxval, cmap=cm.Reds, zorder=1, rasterized=True)
	#m.fillcontinents(color='0.9',lake_color='grey', zorder=4)
	m.drawcoastlines(linewidth=.25,zorder=5)
	m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)

	bbox_args = dict(fc="white")
	ax1.annotate('.                      \n                \n        ', xy=(0.98, 0.98), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', zorder=10)
	
	cax = fig.add_axes([0.76, 0.95, 0.19, 0.03])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
	cbar.set_label('Snow depth (m)', labelpad=4)
	cbar.set_ticks(np.arange(minval, maxval+0.1, 0.1))
	cbar.solids.set_rasterized(True)
	ax1.annotate('Temp scaling '+str(year), xy=(0.03, 0.97),xycoords='axes fraction', backgroundcolor = 'w', horizontalalignment='left', verticalalignment='top', zorder=10)


	subplots_adjust(bottom=0.0, left=0.0, top = 0.99, right=1.0)
	savefig(figpath+'/ICSnowDepth_'+str(year)+'.png', dpi=300)


#-- run main program
if __name__ == '__main__':
	for y in xrange(2016, 2016+1, 1):
		print y
		main(y)

