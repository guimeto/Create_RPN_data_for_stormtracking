#! /bin/sh
set -xa

#dir=/BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_100KM/
dir=/BIG1/dueymes/STORM_TRACK_DATA/REANALYSES/NARR/3_hrs/DATA_850_500/NAM/
dir1=/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_OUT_GUI/NARR/
yini=2018
yfin=2018
#liste_month='01 02 03 11 12'
liste_month='07 08 09 10 11 12'


y=$yini

while [ $y -le $yfin ]; do
  for m in $liste_month ; do  
  rm -rf ${dir1}out*.rpn
  rm -rf  *.fst
  make_tracks_narr  -i  ${dir}NARR_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn -o ${dir1}out_narr_uv850_$y${m} -level 850 -vc 1.5 -npt 8. -silent

  rm -rf *.fst
  mv fort.18 ${dir1}out_narr_uv850_$m${y}_gis.txt

  rm -rf fort.*
  rm -rf *.rpn
  rm -rf *gis
  done
  y=`echo $y | awk '{printf "%04d", $1+1}'`
done

echo " $jobname " ' se termine normalement '           

echo '===> ok <===='

