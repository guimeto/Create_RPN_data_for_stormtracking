#guillaume Dueymes 30 Septembre 2015
#script pour selection un niveau de NARR et preparer le fichier a la conversion netcdf to RPN (IDL)
#
#
set -eax

ymin=2019
ymax=2019

list_mon='01 02 03 04 05 06 07 08 09'


y=$ymin
while [ $y -le $ymax ]; do
 for m in $list_mon ; do 
   
r.diag select /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_ORIGINAL/temporaire_3hrs/hgt/VAR100km/NARR_GZ_ps_100km_${y}_${m}_3hrs.rpn ./tmp.rpn -lv1 500
echo "desire(-1,-1,-1,-1,-1,-1,-1)" >> directive_name
echo "zap(-1,-1,-1,-1,500,-1,-1)" >> directive_name
editfst -s ./tmp.rpn -d /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_100KM/NARR_GZ500_ps100km_${y}_${m}_3hrs.rpn -i directive_name 
rm directive_name  ./tmp.rpn


r.diag select /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_ORIGINAL/temporaire_3hrs/hgt/VAR100km/NARR_GZ_ps_100km_${y}_${m}_3hrs.rpn ./tmp.rpn -lv1 850
echo "desire(-1,-1,-1,-1,-1,-1,-1)" >> directive_name
echo "zap(-1,-1,-1,-1,850,-1,-1)" >> directive_name
editfst -s ./tmp.rpn -d /BIG1/dueymes/WORK/REANALYSES/Amerique_du_nord/NARR/NARR_100KM/NARR_GZ850_ps100km_${y}_${m}_3hrs.rpn -i directive_name 
rm directive_name  ./tmp.rpn
 
 done
y=`echo $y | awk '{printf "%04d", $1+1}'`
done


