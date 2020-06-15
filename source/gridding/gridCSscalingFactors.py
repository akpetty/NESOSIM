''' gridCsscalingFactors.py

	script to grid regional scaling factors for scaling of reanalysis snowfall
	input to CloudSat observations (cf. Cabaj et al 2020, GRL)

	Alex Cabaj, adapted from code by Alek Petty

	input: CloudSat-scaling coefficients by quadrant
	output:

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp2d
import sys
sys.path.append('../')
import utils as cF
import os
import pyproj
import cartopy.crs as ccrs

import pandas as pd


xptsG, yptsG, latG, lonG, proj = cF.create_grid()

print(xptsG, yptsG)

print('grid generated')
print(latG.shape)

# print(latG[0,89])
# print(latG[155,90])
# print(latG[89,155])

# grab subset up to 60 degrees on all sides for interpolation 
# (to be consistent with publication)
latG_subset = latG[23:156,23:156]
lonG_subset = lonG[23:156,23:156]
print(latG_subset.shape)

nx,ny = latG.shape
nxs,nys = latG_subset.shape

# interpolate over this subset

REANs = ['ERAI','ERA5','MERRA_2']
R_FN = ['EI','E5','M2']
R_IDX = 1

# cloudsat scaling factors (as given in Cabaj et al 2020)
cs = pd.read_csv('weights_{}.csv'.format(R_FN[R_IDX]),index_col='time')


# interpolate over centre of grid (to match previous version domain);
# then extend edges south of 60 N and fill corners with constant values

x0=[0,1]
y0=[0,1]

# coordinates of quadrants

x = np.linspace(0,1,nxs)
y = np.linspace(0,1,nys)


scale_factors = np.zeros((12,nx,ny))

for i in range(12):
# for i in range(1):

	# calibration factors
	z0 = np.array([[cs.loc[i+1, 'q4'], cs.loc[i+1, 'q1']], [cs.loc[i+1, 'q3'], cs.loc[i+1, 'q2']]])
	# interpolate over centre (square bounded by 60 N)
	f = interp2d(x0, y0, z0, kind='linear')
	# fill centre with interpolation
	scale_factors[i,23:156,23:156] = f(x,y)

	# complete the rectangle by extending to edges
	# top edge
	print(scale_factors[i,:23,23:156].shape)
	scale_factors[i,:23,23:156] = scale_factors[i,23,23:156]
	# bottom edge
	scale_factors[i,156:,23:156] = scale_factors[i,155,23:156]
	# left edge
	scale_factors[i,23:156,:23] = np.transpose(np.tile(scale_factors[i,23:156,23],(23,1)))
	scale_factors[i,23:156,156:] = np.transpose(np.tile(scale_factors[i,23:156,155],(23,1)))

	#top left corner (fill with constant)
	scale_factors[i,:23,:23]=scale_factors[i,23,23]
	# bottom right corner
	scale_factors[i,156:,156:]=scale_factors[i,155,155]
	# top right
	scale_factors[i,:23,156:]=scale_factors[i,23,155]
	# bottom left
	scale_factors[i,156:,:23]=scale_factors[i,155,23]





	plt.imshow(scale_factors[i],origin='lower',interpolation=None)
	plt.title('Scaling factors for month {}'.format(i+1))
	plt.colorbar()
	# plt.savefig('scaling_factors_full_{}_{}.png'.format(i+1,R_FN[R_IDX]))
	plt.show()
