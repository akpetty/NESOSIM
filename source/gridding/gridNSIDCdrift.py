""" gridNSIDCdrift.py
	
	Script to grid the NSIDC v4 derived Arctic ice drifts
	Model written by Alek Petty (10/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov)

	Input: NSIDC v4 daily ice drifts
	Output: Gridded NSIDC ice drifts

	Python dependencies:
		See below for the relevant module imports
		Also some function in utils.py

	Update history:
		10/01/2020: Version 1
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
from scipy.spatial import Delaunay

# replace with config_[user name] and create your own config file
from config_ap import nsidc_raw_path
from config_ap import forcing_save_path
from config_ap import figure_path


def main(year, startMonth=0, endMonth=11, extraStr='v11', dx=100000, data_path=nsidc_raw_path, out_path=forcing_save_path, fig_path=figure_path+'IceDrift/NSIDCv4d/', anc_data_path='../../AncData/'):

	print(data_path)
	xptsG, yptsG, latG, lonG, proj = cF.create_grid(dxRes=dx)
	# create source and target cartopy projections
	srcProj = cF.EASE_North() 
	newProj = cF.P3413() 

	print(xptsG)
	print(yptsG)

	dxStr=str(int(dx/1000))+'km'
	print(dxStr)

	numDays=cF.getLeapYr(year)
	if (numDays>365):
		monIndex = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	else:
		monIndex = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	if not os.path.exists(out_path+'/'+dxStr+'/IceDrift/NSIDCv4d/'+str(year)):
		os.makedirs(out_path+'/'+dxStr+'/IceDrift/NSIDCv4d/'+str(year))

	if not os.path.exists(fig_path):
		os.makedirs(fig_path)

	calc_weights=1 # start as one to calculate weightings then gets set as zero for future files

	for month in range(startMonth, endMonth+1):
		print(month)
		numDays=monIndex[month]

		for x in range(numDays):
			dayT=sum(monIndex[0:month])+x
			dayStr='%03d' %dayT

			print('Drift, year:', year, 'day:', dayStr)

			xeasedrift, yeasedrift, lont, latt = cF.get_nsidc_driftv4(data_path, year, dayT, 'daily')
			xpts, ypts=proj(lont, latt)
			u_rot, v_rot = newProj.transform_vectors(srcProj, lont, latt, xeasedrift, yeasedrift)
			drift_day_xy=np.stack((u_rot, v_rot))
			xpts, ypts=proj(lont, latt)

			if calc_weights==1:
				# calculate Delaunay triangulation interpolation weightings for first file of the year
				print('calculating interpolation weightings')
				ptM_arr = np.array([xpts.flatten(),ypts.flatten()]).T
				tri = Delaunay(ptM_arr) # delaunay triangulation
				calc_weights=0
			
			
			drift_xyG = cF.int_smooth_drifts_v3(tri, xptsG, yptsG, xpts, ypts, latt, drift_day_xy, sigma_factor=0.5)
			print(drift_xyG)	

			cF.plot_drift_cartopy(lonG , latG , xptsG, yptsG, drift_xyG[0], drift_xyG[1], np.sqrt(drift_xyG[0]**2+drift_xyG[1]**2) , out=fig_path+str(year)+'_d'+dayStr+dxStr+extraStr, units_lab='m/s', units_vec=r'm s$^{-1}$',
				minval=0, maxval=0.5, vector_val=0.1, date_string=str(year)+'_d'+dayStr, month_string='', varStr='NSIDCv4 ice drift ',cbar_type='max', cmap_1=plt.cm.viridis)
				
			drift_xyG.dump(out_path+'/'+dxStr+'/IceDrift/NSIDCv4d/'+str(year)+'/NSIDCv4d_driftG'+dxStr+'-'+str(year)+'_d'+dayStr+extraStr)

#-- run main program
if __name__ == '__main__':
	
	for year in range(2009, 2020+1, 1):
		print(year)
		main(year)
	


	#years=np.arange(2019, 2019+1, 1)
	#from itertools import repeat
	#import concurrent.futures
	#with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:

		#args=((campaign, beam) for beam in beams)
		#print(args)
		#esult=executor.map(main, years)






