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

sys.path.append('../Scripts/Gridding/')
sys.path.append('../Scripts/GetData/')
import getERAIsf
import getERAIwinds
import getERAI2mt
import gridERAIwinds
import gridERAIsf
import gridERAIt2m
import gridCDR
import gridOSISAFdays


# TO DO: ADD A CHECK TO SEE IF DATA EXISTS FIRST!
yearS=2018
monthS=7
dayS=14

year=2019
monthEnd=4 # 12=December
from calendar import monthrange
dayEnd=monthrange(year, monthEnd)[1]

rawDataPath = '/Users/aapetty/Data/'
griddedForcingPath = '/Users/aapetty/Data/Forcings/'
savePath = '/Users/aapetty/Data/NESOSIM/'
figPath = '/Users/aapetty/Data/NESOSIM/Figures/'

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
NESOSIM.main(year1=yearS, month1=monthS, day1=dayS, year2=year, month2=monthEnd-1, day2=dayEnd-1,
	outPathT=savePath, 
	forcingPathT=griddedForcingPath, 
	figPathT=figPath,
	reanalysisP='ERAI', varStr='sf', driftP='OSISAFsig150', team_s='CDR', densityTypeT='variable', 
	outStr='v56', IC=4, windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
	dynamicsInc=1, leadlossInc=1, windpackInc=1)


# To loop over a number of years
# for y in range(2010, 2018+1):
# 	if (y==1987):
# 		continue
# 	print (y)
# 	NESOSIM.main(year1=y, month1=7, day1=14, year2=y, month2=9, day2=30,
# 		outPathT='/Volumes/PETTY_PASSPORT3/NESOSIMdev/Output/', 
# 		forcingPathT='/Volumes/PETTY_PASSPORT3/NESOSIMdev/Forcings/', 
# 		figPathT='/Volumes/PETTY_PASSPORT3/NESOSIMdev/Figures/',
# 		reanalysisP='ERAI', varStr='sf', driftP='OSISAF', team_s='OSISAF', densityTypeT='variable', 
# 		outStr='t2mWindOut', IC=1, windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
# 		dynamicsInc=1, leadlossInc=1, windpackInc=1)









# OLD
# for y in xrange(2016, 2017+1):
# 	if (y==1987):
# 		continue
# 	print y
# 	NESOSIM.main(y, 7, 14, month2=2, day2=30, reanalysisP='ERAI', varStr='sf', driftP='OSISAF', team_s='bt', densityTypeT='variable', 
# 		outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
# 		dynamicsInc=1, leadlossInc=1, windpackInc=1)

# NESOSIM.main(2013, 7, 14, reanalysisP='ERAI', varStr='sf', driftP='NSIDCv3', team_s='bt', densityTypeT='variable', 
# 	outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
# 	dynamicsInc=1, leadlossInc=1, windpackInc=1)
