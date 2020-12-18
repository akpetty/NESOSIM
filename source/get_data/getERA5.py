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

c = cdsapi.Client()

def main(year, month):
   
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':'snowfall', 
            'year':year,
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
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ]
        },
        reanalysis_raw_path+'ERA5_sf_'+year+month+'cds.nc')

if __name__ == '__main__':
    

    x='2020'
    for m in range(2, 3+1):
        mStr='%02d' %(m+1)
        print(x, mStr)
        main(x, mStr)



    #months=['01', '02', '03', '04']
    #import concurrent.futures
    #from itertools import repeat
    #with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
    #    result=executor.map(main, repeat(2020), months)
    #years=np.arange(2000, 2016+1, 1)
    #from itertools import repeat
    #import concurrent.futures
    #from itertools import repeat
    #with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
    
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
    






