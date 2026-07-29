#!/bin/bash
# Download NOAA/NSIDC CDR daily sea ice concentration (G02202 Version 6).
# Usage: ./getCDR.sh /path/to/ICECONC/CDR/daily/ YYYY MM
# Saves to: $save_path/final/v6/$year/

save_path=$1
year=$2
month=$3

echo $save_path
echo $year
echo $month

FILE=https://noaadata.apps.nsidc.org/NOAA/G02202_V6/north/daily/"$year"/
echo $FILE

mkdir -p "$save_path"final/v6/"$year"/

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=5 --directory-prefix=$save_path"final/v6/$year/" --reject "index.html*" \
-np -e robots=off -r -A "*$year$month*.nc" $FILE
