#!/bin/bash
# Download Near-Real-Time NOAA/NSIDC CDR daily sea ice concentration (G10016 Version 4).
# Usage: ./getCDRnrt.sh /path/to/ICECONC/CDR/daily/ YYYY MM
# Saves to: $save_path/nrt/v4/$year/

save_path=$1
year=$2
month=$3

echo $save_path
echo $year
echo $month

FILE=https://noaadata.apps.nsidc.org/NOAA/G10016_V4/north/daily/"$year"/
echo $FILE

mkdir -p "$save_path"nrt/v4/"$year"/

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=5 --directory-prefix=$save_path"nrt/v4/$year/" --reject "index.html*" \
-np -e robots=off -r -A "*$year$month*.nc" $FILE
