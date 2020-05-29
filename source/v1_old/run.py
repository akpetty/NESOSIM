""" run.py
	
	Run script for the NESOSIM model included in NESOSIM.py 
	Model written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov) or add a query to the GitHub repo (www.github.com/akpetty/NESOSIM)
	
	To do:
		Add a check to see if data exists

	Update history:
		03/01/2018: Version 1
		03/01/2018: Version 2

"""

import matplotlib
matplotlib.use("AGG")
from pylab import *


import subprocess
import shlex
import sys

rawDataPath = '/Users/aapetty/Data/'
griddedForcingPath = '/Users/aapetty/Data/Forcings/'
savePath = '/Users/aapetty/Data/NESOSIM/'
figPath = '/Users/aapetty/Data/NESOSIM/Figures/'

for y in range(2010, 2018+1):
	#===== GET DATA =========

#----- Snowfall/winds ---------
#getERAIsf.main(year, monthEnd, rawDataPath) # last number if the end month, so 12=December
#getERAIwinds.main(year, monthEnd, rawDataPath) # last number if the end month, so 12=December
#getERAI2mt.main(year, monthEnd, rawDataPath) # last number if the end month, so 12=December

#----- Sea ice Concentration ---------
#subprocess.call(["./GetData/getCDRnrt.sh", str(year)])

#----- Sea ice Drift ---------
#subprocess.call(["./GetData/getOSISAF.sh", str(year)])

#===== GRID DATA =========

#----- Reanalyses ---------
#gridERAIsf.main(year, grid_res, startMonth=3, endMonth=monthEnd-1, extraStr='temp', grid_res=50)
#gridERAIwinds.main(year, grid_res, startMonth=3, endMonth=monthEnd-1, extraStr='temp', grid_res=50)
#gridERAIt2m.main(year, grid_res, startMonth=0, endMonth=monthEnd-1, extraStr='temp', grid_res=50)

#----- Ice Conc ---------
#gridCDR.main(year, startMonth=3, endMonth=monthEnd-1, grid_res=50)

#----- Ice Drift ---------
#gridOSISAFdays.main(year, startMonth=3, endMonth=monthEnd-1, grid_res=50)



import NESOSIM

#outPath='/Volumes/PETTY_PASSPORT3/NESOSIM/Output/'
#forcingPath='../Forcings/'

# To loop over a number of years
for y in range(2010, 2018+1):
	if (y==1987):
		# Sea ice concentration data gap
		continue
	print (y)
	NESOSIM.main(year1=y, month1=7, day1=14, year2=y+1, month2=4, day2=0,
		outPathT=outPath, 
		forcingPathT=forcingPath, 
		figPathT=figPath,
		reanalysisP='ERAI', varStr='sf', driftP='OSISAFsig150', team_s='bt', densityTypeT='variable', 
		outStr='t2mWindOut', IC=1, windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
		dynamicsInc=1, leadlossInc=1, windpackInc=1)



