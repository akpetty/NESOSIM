""" run_oneseason.py
	
	Run script for the NESOSIM model included in NESOSIM.py 
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

from config import forcing_save_path
from config import model_save_path
from config import figure_path

print('Forcing file path:', forcing_save_path)
print('Output path:', model_save_path)
print('Figure save path:', figure_path)

yearS=2018
monthS=8 # August = 7
dayS=0

yearE=2019
monthE=3 # April = 7
dayE=29

print(yearS, monthS, dayS, yearE, monthE, dayE)

import NESOSIM	
NESOSIM.main(year1=yearS, month1=monthS, day1=dayS, year2=yearE, month2=monthE, day2=dayE,
	outPathT=model_save_path, 
	forcingPathT=forcing_save_path, 
	figPathT=figure_path,
	precipVar='ERA5', windVar='ERA5', driftVar='NSIDCv4', concVar='CDR', 
	icVar='ERA5', densityTypeT='variable', extraStr='v11', outStr='test', IC=2, 
	windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
	dynamicsInc=1, leadlossInc=1, windpackInc=1,scaleCS=True, dx=100000,
	plotdaily=1)




