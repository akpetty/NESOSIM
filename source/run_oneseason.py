""" run.py
	
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

yearS=2010
monthS=8 # August = 7
dayS=1

yearE=2016
monthE=4 # 4 = april
from calendar import monthrange
# Find last day in that month
dayE=monthrange(yearE, monthE+1)[1]
print(yearS, monthS, dayS, yearE, monthE, dayE)

import NESOSIM	
NESOSIM.main(year1=yearS, month1=monthS, day1=dayS, year2=yearE, month2=monthE, day2=dayE-1,
	outPathT=model_save_path, 
	forcingPathT=forcing_save_path, 
	figPathT=figure_path+'Model/',
	precipVar='ERAI', windVar='ERAI', driftVar='OSISAF', concVar='CDR', 
	icVar='ERAI', densityTypeT='variable', extraStr='v11', outStr='4x_v2_s03', IC=2, 
	windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
	dynamicsInc=1, leadlossInc=1, windpackInc=1,scaleCS=False, dx=100000)




