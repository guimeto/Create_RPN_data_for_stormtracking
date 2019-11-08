#guillaume Dueymes 16 Septembre 2015
#script pour telecharger les donnees de Reanalyse NARR
set -eax

ymin=2018
ymax=2019

#######VARIABLES AUX 3HEURES NIVEAUX PRESSION
list_var='hgt uwnd vwnd'
#######VARIABLES AUX 3HEURES MONOLEVEL
#list_var2='air.2m'
#######VARIABLES AUX 3HEURES NIVEAUX PRESSION
for var in $list_var ;do
y=$ymin
while [ $y -le $ymax ]; do
if [ "$y" == "2018" ]; then
	echo "ANNEE COURANTE NON TERMINEE!"
list_mon='06 07 08 09 10 11 12'
else
	echo "ANNEE COMPLETE!"
list_mon='01 02 03 04 05 06 07 08 09'
fi
for m in $list_mon ; do 
wget ftp://ftp.cdc.noaa.gov/Datasets/NARR/pressure/${var}.${y}${m}.nc
##cdo sellevel,500  ${var}.${y}${m}.nc ${var}500.${y}${m}.nc
##cdo sellevel,850  ${var}.${y}${m}.nc ${var}850.${y}${m}.nc
mv ${var}.${y}${m}.nc  ./temporaire_3hrs/NARR_${var}_lc_${y}_${m}_3hrs.nc
##mv  ${var}850.${y}${m}.nc ./temporaire_3hrs/NARR_${var}850_lc_${y}_${m}_3hrs.nc
##mv  ${var}500.${y}${m}.nc ./temporaire_3hrs/NARR_${var}500_lc_${y}_${m}_3hrs.nc
#cdo sellevel,1000  ${var}.${y}${m}.nc ${var}1000.${y}${m}.nc

#rm -rf ${var}.${y}${m}.nc
 done
y=`echo $y | awk '{printf "%04d", $1+1}'`
done
done


#######VARIABLES AUX 3HEURES MONOLEVEL
for var2 in $list_var2 ;do
y=$ymin
while [ $y -le $ymax ]; do

wget ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/${var2}.${y}.nc
cdo splitmon ${var2}.${y}.nc  ./${var2}.${y}
rm ${var2}.${y}.nc 
if [ "$y" == "2019" ]; then
	echo "ANNEE COURANTE NON TERMINEE!"
list_mon='01 02 03 04 05 06'
else
	echo "ANNEE COMPLETE!"
list_mon='01 02 03 04 05 06 07 08 09 10 11 12'
fi

for m in $list_mon ; do 
#mv ./${var2}.${y}${m}.nc /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_ORIGINAL/temporaire_3hrs/NARR_${var2}_lc_${y}_${m}_3hrs.nc
mv ./${var2}.${y}${m}.nc /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_ORIGINAL/temporaire_3hrs/NARR_tas_lc_${y}_${m}_3hrs.nc

done

y=`echo $y | awk '{printf "%04d", $1+1}'`
done
done

