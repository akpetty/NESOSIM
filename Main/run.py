""" run.py
	
	Run script for the NESOSIM model included in NESOSIM.py 
	Model written by Alek Petty (03/01/2018)
	Contact me for questions (alek.a.petty@nasa.gov) or add a query to the GitHub repo (ADD THIS)

	Update history:
		03/01/2018: Version 1

"""

import matplotlib
matplotlib.use("AGG")
from pylab import *
import NESOSIM


NESOSIM.main(2010, 7, 14, reanalysisP='ERAI', varStr='sf', driftP='NSIDCv3', team_s='bt', densityTypeT='variable', 
	outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
	dynamicsInc=1, leadlossInc=1, windpackInc=1)


# To loop over a number of years
#for y in xrange(2000, 2014):
# 	if (y==1987):
# 		continue
# 	print y
# 	NESOSIM.main(y, 7, 14, reanalysisP='MEDIAN', varStr='sf', driftP='NSIDCv3', team_s='nt', densityTypeT='variable', 
# 		outStr='', IC=1, windPackFactorT=0.05, windPackThreshT=5, leadLossFactorT=0.025,
# 		dynamicsInc=1, leadlossInc=1, windpackInc=1)


