#!/bin/bash
#ftp ftp://osisaf.met.no/prod/ice/drift_lr/merged/
save_path=''

year='2009'
echo $year

wget -r -np -nH --cut-dirs=5 --reject "index.html*" --directory-prefix=$save_path"$year"/ ftp://osisaf.met.no/archive/ice/drift_lr/merged/"$year"/
cp -r $save_path$year/osisaf.met.no/ $save_path$year/
rm -r $save_path$year/osisaf.met.no
