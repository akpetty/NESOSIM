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


# TO DO: ADD A CHECK TO SEE IF DATA EXISTS FIRST!


# DOWNLOAD ERAI SF DATA
#sys.path.append('../Scripts/GetData/')

#import getERAIsf
#getERAIsf.main(year, 12) # last number if the end month, so 12=December

#sys.path.append('../Scripts/GetData/')
#import getERAIsf
#getERAIsf.main(year, 12) # last number if the end month, so 12=December


import NESOSIM

#outPath='/Volumes/PETTY_PASSPORT3/NESOSIM/Output/'
#forcingPath='../Forcings/'

# To loop over a number of years
for y in range(2010, 2016+1):
	if (y==1987):
		continue
	print (y)
	NESOSIM.main(y, 7, 14, 
		outPathT='/Volumes/PETTY_PASSPORT3/NESOSIM/Output/', 
		forcingPathT='/Volumes/PETTY_PASSPORT3/NESOSIM/Forcings/', 
		reanalysisP='ERAI', varStr='sf', driftP='NSIDCv3', team_s='bt', densityTypeT='variable', 
		outStr='BSThresh', IC=1, windPackFactorT=5.8e-7, windPackThreshT=5, leadLossFactorT=2.9e-7,
		dynamicsInc=1, leadlossInc=1, windpackInc=1)











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
