#!/usr/bin/env python
import cdsapi
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
import sys
sys.path.append('../')
from config import reanalysis_raw_path

# Need to set this up using https://cds.climate.copernicus.eu/api-how-to


c = cdsapi.Client()

def main(year, month):
   
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':['10m_u_component_of_wind','10m_v_component_of_wind'],
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
        reanalysis_raw_path+"ERA5_winds_"+year+month+'cds.nc')

for x in range(2010, 2018+1):
    for m in range(0, 11+1):
        mStr='%02d' %(m+1)
        main(str(x), mStr)






