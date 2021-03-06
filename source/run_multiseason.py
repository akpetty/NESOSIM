""" run_multiseason.py
	
	Run script for the NESOSIM model included in NESOSIM.py over multiple years
	Model written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov) or add a query to the GitHub repo (www.github.com/akpetty/NESOSIM)

	Update history:
		03/01/2018: Version 1

"""

import matplotlib
matplotlib.use("AGG")
from pylab import *


import subprocess
import shlex
import sys
from calendar import monthrange

from config import forcing_save_path, model_save_path, figure_path
print('Forcing file path:' forcing_save_path)
print('Output path:', model_save_path)
print('Figure save path:', figure_path)

import NESOSIM


yearS=2009
yearE=2019

monthS = 8 # August = 7, September = 8
monthE = 3 # 3 = April

dayS=0 # 1st day of month
dayE=monthrange(yearE, monthE+1)[1] # last day of the month

for y in range(yearS, yearE+1):
	print(y, monthS, dayS, y+1, monthE, dayE)

	NESOSIM.main(year1=y, month1=monthS, day1=dayS, year2=y+1, month2=monthE, day2=dayE-1,
		outPathT=model_save_path, 
		forcingPathT=forcing_save_path, 
		figPathT=figure_path+'Model/',
		precipVar='ERA5', windVar='ERA5', driftVar='OSISAF', concVar='CDR', 
		icVar='ERA5', densityTypeT='variable', extraStr='v11', outStr='oct28', IC=2, 
		windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=1.45e-7,
		dynamicsInc=1, leadlossInc=1, windpackInc=1, atmlossInc=1, plotBudgets=1, plotdaily=0,
		scaleCS=True, dx=100000)




