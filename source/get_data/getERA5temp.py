#!/usr/bin/env python
import cdsapi
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
import sys
sys.path.append('../')
from config import reanalysis_save_path

# in ERAI we needed to forecast how much precip comes down over the next 12 hours (hence the use of step)
# but I think in ERA-5 we get the hourly accumulation from the previous step (so for time 01:00 this would be the precip from 00:00 to 01:00)

# Need to switch over to https://cds.climate.copernicus.eu/api-how-to eventually as other method being discontinued
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+How+to+calculate+daily+total+precipitation
c = cdsapi.Client()

def main(year, month):
   
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':'2m_temperature', 
            'year':str(year),
            'month':month,
            'day':[
                '01','02','03',
                '04','05','06',
                '07','08','09',
                '10','11','12',
                '13','14','15',
                '16','17','18',
                '19','20','21',
                '22','23','24',
                '25','26','27',
                '28','29','30',
                '31'
            ],
            'time':[
                '00:00','06:00',
                '12:00','18:00'
            ]
        },
        reanalysis_save_path+"ERA5_temp6hour_"+str(year)+month+'cds.nc')

if __name__ == '__main__':
    
    for x in range(1985, 1985+1):
        for m in range(0, 1+1):
            mStr='%02d' %(m+1)
            main(x, month=mStr)


    #years=np.arange(2018, 2019+1, 1)
    #from itertools import repeat
    #import concurrent.futures
    #from itertools import repeat
    #with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    
        #args=((campaign, beam) for beam in beams)
        #print(args)
        #result=executor.map(main, years, repeat('JRA55'))
        #result2=executor.map(main, years, repeat('MERRA'))
        #result3=executor.map(main, years, repeat('MERRA_2'))
        #result4=executor.map(main, years, repeat('NCEP_R1'))
        #result5=executor.map(main, years, repeat('NCEP_R2'))
        #result6=executor.map(main, years, months=['01'])
        #result0=executor.map(main, years,repeat('01'))
        #result1=executor.map(main, years,repeat('02'))
        #result2=executor.map(main, years,repeat('03'))
        #result3=executor.map(main, years,repeat('04'))
        #result4=executor.map(main, years,repeat('05'))
        #result5=executor.map(main, years,repeat('06'))
        #result6=executor.map(main, years,repeat('07'))
        #result7=executor.map(main, years, repeat('08'))
        #result8=executor.map(main, years, repeat('09'))
        #result9=executor.map(main, years, repeat('10'))
        #result10=executor.map(main, years, repeat('11'))
        #result11=executor.map(main, years, repeat('12'))
        #result11=executor.map(main, years, '07')
        #result12=executor.map(main, years, '08')
        #result13=executor.map(main, years, '09')
        #result14=executor.map(main, years, '10')
        #result15=executor.map(main, years, '11')
        #result16=executor.map(main, years, '12')
        


#for x in range(1990, 2000):
#    for m in range(0, 11+1):
#        mStr='%02d' %(m+1)
#        getERA5cds(x, months=[mStr])






