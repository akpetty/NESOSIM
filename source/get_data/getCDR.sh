#!/bin/bash

save_path=$1
year=$2
month=$3

echo $save_path
echo $save_path"nrt/$year/$month/"
echo $year
echo $month

FILE=ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V4/north/daily/"$year"/
echo $FILE

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=7 --directory-prefix=$save_path"final/$year/$month/"  --reject "index.html*" \
-np -e robots=off -r -A "*$year$month*.nc" $FILE
