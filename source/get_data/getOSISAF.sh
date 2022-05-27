#!/bin/bash
# FINAL data from this ftp site ftp://osisaf.met.no/archive/ice/ice/drift_lr/merged/
# NRT data from this ftp site ftp://osisaf.met.no/prod/ice/drift_lr/merged/
# Only grabbing NH data (remove -A option if you want both hemispheres)

# run as ./getOSISAF.sh /cooler/sea_ice_pso/aapetty/raw_data/ 2020 01

save_path=$1
year=$2
month=$3

echo $save_path
echo $year
echo $month


wget -r -np -nH --cut-dirs=5 --reject "index.html*" --directory-prefix=$save_path"$year"/ -A "*nh*.nc" ftp://osisaf.met.no/archive/ice/drift_lr/merged/"$year"/"$month"/
#cp -r $save_path$year/osisaf.met.no/ $save_path$year/
#rm -r $save_path$year/osisaf.met.no
