''' gridCsscalingFactors.py

	script to grid regional scaling factors for scaling of reanalysis snowfall
	input to CloudSat observations (cf. Cabaj et al 2020, GRL)

	Alex Cabaj, adapted from code by Alek Petty

	input: CloudSat-scaling coefficients by quadrant
	output:

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import sys
sys.path.append('../')
import utils as cF
import cartopy.crs as ccrs
import pandas as pd
import xarray as xr


ancDataPath = '../../anc_data/'
xptsG, yptsG, latG, lonG, proj = cF.create_grid()

# print(xptsG, yptsG)

print('grid generated')
print(latG.shape)


# use centre of grid up to 60 degrees latitude on all sides for interpolation 
latG_centre = latG[23:156,23:156] # hardcoding these indices for now
print('centre grid shape')
print(latG_centre.shape)

nx,ny = latG.shape
nxs,nys = latG_centre.shape

REANs = ['ERAI','ERA5','MERRA_2']
R_FN = ['EI','E5','M2']
R_IDX = 1 # which reanalysis to select

# cloudsat scaling factors (as given in Cabaj et al 2020)
cs = pd.read_csv('{}weights_{}.csv'.format(ancDataPath,R_FN[R_IDX]),index_col='time',comment='#')


# interpolate over centre of grid (to match previous NESOSIM version domain);
# then extend edges south of 60 N and fill corners with constant values

# values for interpolation
x0=[0,1]
y0=[0,1]

# coordinates of quadrants for interpolation
x = np.linspace(0,1,nxs)
y = np.linspace(0,1,nys)

# grid of monthly scaling factors for each month
scale_factors = np.zeros((12,nx,ny))

for i in range(12):

	# monthly scaling factors; load into an array with values in corners; 
	# values to be interpolated over the centre of the model domain
	z0 = np.array([[cs.loc[i+1, 'q4'], cs.loc[i+1, 'q3']], [cs.loc[i+1, 'q1'], cs.loc[i+1, 'q2']]])
	# interpolate over centre (square bounded by 60 N)
	f = interp2d(x0, y0, z0, kind='linear')
	# fill centre with interpolation
	scale_factors[i,23:156,23:156] = f(x,y)

	# pad the rest of the domain by extending to the edges
	# todo: would be cleaner with numpy.pad, but this works for now

	# top edge
	scale_factors[i,:23,23:156] = scale_factors[i,23,23:156]
	# bottom edge
	scale_factors[i,156:,23:156] = scale_factors[i,155,23:156]
	# left edge
	scale_factors[i,23:156,:23] = np.transpose(np.tile(scale_factors[i,23:156,23],(23,1)))
	# right edge
	scale_factors[i,23:156,156:] = np.transpose(np.tile(scale_factors[i,23:156,155],(23,1)))

	# extend to corners; fill corners with constants
	#top left corner
	scale_factors[i,:23,:23]=scale_factors[i,23,23]
	# bottom right corner
	scale_factors[i,156:,156:]=scale_factors[i,155,155]
	# top right
	scale_factors[i,:23,156:]=scale_factors[i,23,155]
	# bottom left
	scale_factors[i,156:,:23]=scale_factors[i,155,23]

	# plot to double-check (optional)
	plt.imshow(scale_factors[i],origin='upper',interpolation=None)
	plt.title('Scaling factors for month {}'.format(i+1))
	plt.colorbar()
	# plt.savefig('scaling_factors_full_{}_{}.png'.format(i+1,R_FN[R_IDX]))
	plt.show()

# output gridded scaling factors

# convert to dataarray for saving as netcdf
scale_da = xr.DataArray(scale_factors, dims=['time','x','y'],coords={'time':cs.index},name='scale_factors')
print(scale_da)

# save to netcdf file
scale_da.to_netcdf('{}scale_coeffs_{}.nc'.format(ancDataPath,REANs[R_IDX]))