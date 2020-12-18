#!/bin/bash
year='2020'
SAVEPATH='/Volumes/PETTY_E/DATA/ICE_CONC/CDR/'
# nrt ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G10016/north/daily/
# final ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V3/north/daily/
FILE=ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G10016/north/daily/"$year"/
echo $FILE
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=6 --directory-prefix=$SAVEPATH  --reject "index.html*" \
-np -e robots=off -r -A "*$year*.nc" $FILE
