#!/usr/bin/env python
import cdsapi
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
import sys
sys.path.append('../')
from config import reanalysis_raw_path

# Need to set this up using https://cds.climate.copernicus.eu/api-how-to

# In ERAI we needed to forecast how much precip comes down over the next 12 hours (hence the use of step)
# but I think in ERA-5 we get the hourly accumulation from the previous step (so for time 01:00 this would be the precip from 00:00 to 01:00)

# Note this on expver for recent data when there is a cross over between some experimental and final data: https://confluence.ecmwf.int/pages/viewpage.action?pageId=173385064

def main(year, month2=12):
    
    day2=calendar.monthrange(year, month2)[1]
    date1 = datetime(year, 1, 1)
    date2 = datetime(year, month2, day2)
    dateStr1=date1.strftime('%Y-%m-%d')
    dateStr2=date2.strftime('%Y-%m-%d')
    print dateStr1, dateStr2
    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "expver": "1",
        "levtype": "sfc",
        #"number": "0",
        # for codes see http://apps.ecmwf.int/codes/grib/param-db
        # 228.128=tp, 144.128=sf
        "param": "144.128",
        # as we are getting total precip, we need to forecast how much precip comes down over 
        # the next 12 hours (hence the step)
        "step": "12",
        # in 0.5 degrees lat/lon
        "grid": "0.5/0.5",
        "stream": "oper",
        "date": ""+dateStr1+"/to/"+dateStr2,
        # see above comment. 
        "time": "00:00:00/12:00:00",
        "type": "fc",
        "format" : "netcdf",
        # set an output file name
        "target": ""+reanalysisDataPath+"ERAI_sf_"+str(year)+"temp.nc",

    })

if __name__ == '__main__':
    
    for x in range(2017, 2019+1):
        main(x, month2=12)





