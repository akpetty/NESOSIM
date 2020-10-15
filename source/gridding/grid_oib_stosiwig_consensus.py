
""" grid_oib_stosiwig_consensus.py
	
	Script to grid the STOSIWIG raw OIB snow picks to the NESOSIM grid and generate a consensus (median) gridded value
	Model code written by Alek Petty (05/01/2020)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: OIB quicklook data
	Output: Gridded OIB snow depths

	Python dependencies:
		See below for the relevant module imports
		Also reads in some functions in utils.py

"""


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import griddata
import sys
sys.path.append('../')
import utils as cF
import os
import pyproj
import cartopy.crs as ccrs

from config import forcing_save_path, figure_path, oib_data_path, anc_data_path


def main(year, dx=100000, extraStr='v11'):
	
	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	#print(xptsG)
	#print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	#print(dxStr)

	files = glob(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+'*'+dxStr+extraStr+'')
	print(files)
	oibdates=[file[-16:-8] for file in files]

	print(oibdates)
	for x in range(len(oibdates)):
		# Loop through dates of each flight. I want to keep separate to compre to the daily NESOSIM data.
		
		try:
			oib_dayG1=np.load(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+'GSFC'+extraStr, allow_pickle=True)
			oib_dayG2=np.load(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+'SRLD'+extraStr, allow_pickle=True)
			oib_dayG3=np.load(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+'JPL'+extraStr, allow_pickle=True)

			oib_dayGC=np.median((oib_dayG1, oib_dayG2, oib_dayG3), axis=0)
			cF.plot_gridded_cartopy(lonG, latG, oib_dayGC, proj=ccrs.NorthPolarStereo(central_longitude=-45), out=figure_path+'/OIB/'+oibdates[x]+dxStr+'MEDIAN'+extraStr, date_string=str(year), month_string='', varStr='OIB snow depth MEDIAN', units_lab=r'm', minval=0, maxval=0.6, cmap_1=plt.cm.viridis)
		
			oib_dayGC.dump(forcing_save_path+dxStr+'/OIB/'+str(year)+'/'+oibdates[x]+dxStr+'MEDIAN'+extraStr)
		except:
			print('All thre algs dont exist for this date:', oibdates[x])

#-- run main program
if __name__ == '__main__':
	for y in range(2012, 2019+1, 1):
		print (y)
		main(y)
