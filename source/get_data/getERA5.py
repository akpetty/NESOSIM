""" getERA5.py
    
    Script to grab relevant ERA5 forcings
    Model written by Alek Petty (03/01/2022)
    Contact me for questions (alek.a.petty@nasa.gov) or add a query to the GitHub repo (www.github.com/akpetty/NESOSIM)
    
    - Need to set this up using https://cds.climate.copernicus.eu/api-how-to
    - In ERAI we needed to forecast how much precip comes down over the next 12 hours (hence the use of step). But in ERA-5 we get the hourly accumulation from the previous step (so for time 01:00 this would be the precip from 00:00 to 01:00)
    - Note this on expver for recent data when there is a cross over between some experimental and final data: https://confluence.ecmwf.int/pages/viewpage.action?pageId=173385064

    Update history:
        03/01/2022: Version 1

"""


import cdsapi
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
import sys
sys.path.append('../')
from config_cooler import reanalysis_raw_path


c = cdsapi.Client()

def main(year, month, var='snowfall', save_path='./'):
    """ Main download function

    Use a variable leadlossfactor parameter. This is relatively unconstrained!

    Args:
        year (str): year date
        month (str): month date
        var (str): variable of interest
        save_path (str): download path

    returns:
        raw data saved to save_path

    """
    
    if var=='snowfall':
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
        save_path+'ERA5_sf_'+year+month+'cds.nc')

    if var=='winds':
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
        save_path+"ERA5_winds_"+year+month+'cds.nc')

    if var=='temp':
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
        save_path+"ERA5_temp6hour_"+str(year)+month+'cds.nc')


if __name__ == '__main__':    

    x='2021'
    for m in range(7, 11+1):
        mStr='%02d' %(m+1)
        print(x, mStr)
        main(x, mStr, var='snowfall', save_path=reanalysis_raw_path)







