#################
# Guillaume Dueymes 30 Septembre 2015
# interpolation des donnees NARR de 32km a 100km
# utilisation de pgsm et gritable
#
#################
set -xa

res=100km
ymin=2018
ymax=2018

liste_month='07 08 09 10 11 12'

y=$ymin
while [ $y -le $ymax ] ; do
for m in $liste_month ; do
 echo $y

############ Partie de code VV
echo " VOIRENT=OUI" >> directive
echo " VOIRSRT=OUI" >> directive
echo " SORTIE(STD,5000,A)" >> directive
echo " GRILLE(PS,85,72,27.5,81.,100000.0,25.0,NORD) " >> directive
echo " COMPAC=-16" >> directive
echo " IP2(TOUT)" >> directive
#echo " EXTRAP(VOISIN)" >> directive   ############################################################
echo " SETINTX(LINEAIR)" >> directive
echo " HEURE(-1)" >> directive
echo " CHAMP(GZ,ALL)" >> directive
echo " END" >> directive

pgsm  -iment ./NARR_hgt_lc_${y}_${m}_3hrs.fst   -ozsrt ./VAR100km/NARR_GZ_ps_100km_${y}_${m}_3hrs.rpn  -i directive
rm -rf  directive


done

y=`echo $y | awk '{printf "%04d", $1+1}'`
done
