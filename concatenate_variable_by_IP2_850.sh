diag #! /bin/sh

## Guillaume Dueymes 24 Juillet 2015 
#
# script qui effectue les operations suivantes sur un fichier RPN
# - 1: interpolation des niveaux sigma sur les niveaux pression 
# - 2: ecriture des IP2 en fonction de la date 
#  INPUT: le present code travaille avec des fichiers mensuels 
#
set -eax

y_ini=2018
y_end=2018

List_month='07 08 09 10 11 12'
model='NARR'

################################## PARTIE ############################
# Partie  : cette section  permet d ouvrir les fichiers sur les IP2 
 
y=$y_ini
 while [ $y -le $y_end ]; do          #debut de la boucle sur les annees
ip2=0
echo 'year' $y

for m in $List_month;  do        #debut de la boucle sur les mois


 DATE=${y}${m}0100000000
  CMCSTAMP=`r.date -Sn ${DATE}`  #definition de la date CMC: date qui est lue dans le fichier RPN 
ip2=0
if [[ $m -eq 01 || $m -eq 03 || $m -eq 05 || $m -eq 07 || $m -eq 08 || $m -eq 10 || $m -eq 12 ]]; then
ipmax=741
fi
if [[ $m -eq 04 || $m -eq 06 || $m -eq 11 ]]; then
ipmax=717
fi
if [[ $m -eq 09  ]]; then
ipmax=717
fi
if [ $m -eq 02 ]; then
if [[ $y -eq 1980 || $y -eq 1984 || $y -eq 1988 || $y -eq 1992 || $y -eq 1996 || $y -eq 2000 || $y -eq 2004 || $y -eq 2008 || $y -eq 2012 || $y -eq 2016  ]]; then
echo $(( ipmax=693))
else
echo $(( ipmax=669))
fi
fi

  while [ $ip2 -le $ipmax ]; do          #debut de la boucle sur les annees
echo 'ip2' $ip2
   
echo "desire(-1,GZ,-1,-1,850,$ip2,-1)" >> directive_name
editfst -s ./GZ850/${model}_GZ850_ps100km_${y}_${m}_3hrs.rpn  -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn   -i directive_name
rm directive_name
echo "desire(-1,UU,-1,-1,850,$ip2,-1)" >> directive_name
editfst -s ./UU850/${model}_UU850_ps100km_${y}_${m}_3hrs.rpn  -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn   -i directive_name
rm directive_name
echo "desire(-1,UU,-1,-1,500,$ip2,-1)" >> directive_name
editfst -s ./UU500/${model}_UU500_ps100km_${y}_${m}_3hrs.rpn  -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn   -i directive_name
rm directive_name
echo "desire(-1,VV,-1,-1,850,$ip2,-1)" >> directive_name
editfst -s ./VV850/${model}_VV850_ps100km_${y}_${m}_3hrs.rpn  -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn  -i directive_name
rm directive_name
echo "desire(-1,VV,-1,-1,500,$ip2,-1)" >> directive_name
editfst -s ./VV500/${model}_VV500_ps100km_${y}_${m}_3hrs.rpn  -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn  -i directive_name
rm directive_name 

ip2=`expr $ip2 + 3`
done ###boucle sur les IP2
#on ajoute la topographie
echo "desire(-1,MX,-1,-1,-1,-1,-1)" >> directive_name
editfst -s ./NARR_OROG_100km.rpn -d ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn  -i directive_name
rm directive_name 

#on ajoute la date 
 echo "desire(-1,-1,-1,-1,-1,-1,-1)" >> directive
 echo "zap(-1,-1,-1,${CMCSTAMP},-1,-1,-1)" >> directive
 editfst -s ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn  -d ./${model}_GZ_UU_VV_lc_${y}${m}_new_3hrs_100km.rpn -i directive
rm ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn 
mv ./${model}_GZ_UU_VV_lc_${y}${m}_new_3hrs_100km.rpn ./${model}_GZ_UU_VV_850_lc_${y}${m}_3hrs_100km.rpn 
 rm -rf directive

done ###boucle sur les mois
y=`expr $y + 1`
done ###boucle sur les annees
echo "STOP"
################################################################################################################
