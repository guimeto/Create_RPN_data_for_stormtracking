#!/bin/sh


month_list="01 02 03 04 05 06 07 08 09 08 09"

year_min=2019
year_max=2019




for month in ${month_list}

  do
       year=${year_min}  
       while [ ${year} -le ${year_max} ]
         do 
            MOIS=${month}
            ANNEE=${year}
            export MOIS
            export ANNEE 
            echo "$ANNEE - $MOIS"
            /unique/logiciels/beluga/gdl-0.9.2/bin/gdl -e statistique
            year=$[year+1]
         done


  done


