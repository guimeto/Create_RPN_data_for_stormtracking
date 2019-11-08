#!/usr/bin/env python
#! 
#!
#!
#! Guillaume Dueymes 23/10/2017
#! programme qui lit les donnees de type netcdf
#! et ecrit les donnees en format type RPN
#!
#! Convetion des variables:
#!     netcdf: air       =   RPN: TT
#!     netcdf: hgt       =   RPN: GZ
#!     netcdf: uwnd      =   RPN: UU;     
#!     netcdf: vwnd      =   RPN: VV
#!     netcdf: apcp      =   RPN: PR
#!     netcdf: acpcp     =   RPN: PC
#!     netcdf: shum      =   RPN: HU 
#!     netcdf: omega     =   RPN: WW 
#!     netcdf: prate     =   RPN: RT      
#!     netcdf: prmsl     =   RPN: PN 
#!     netcdf: pres.sfc  =   RPN: P0 
#!     netcdf: shtfl     =   RPN: FC 
#! 
#! first release of the python_rpn package (interfaces to rmnlib and vgrid)
#! . s.ssmuse.dot rpnpy
#! import rpnpy.vgd.all as vgd
#! import rpnpy.librmn.all as rmn  
# 
#from __future__ import division
import rpnpy.vgd.all as vgd
import rpnpy.librmn.all as rmn
import matplotlib.patches as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
import netCDF4
from netCDF4 import Dataset
import rpnstd
import sys
import re
import Fstdc
import os
##################################################
# Housekeeping
plt.close('all')
# Fst constants 
DEET = 0
NPAS = 0
GridType = "!"
DATEO = 0
IP3 = 0
IG1 = 221
IG2 = 0
IG3 = 0
IG4 = 0
VNAME ='GZ'
#---------------open netcdf and fst file with same file name--------------------------------
# Open netcdf file
mois =  ['07','08','09','10','11','12']
yi=2018
yf=2018

for year in range(yi,yf+1):
    for moi in mois: 

        filecdf_path = 'NARR_hgt_lc_'+str(year)+'_'+moi+'_3hrs.nc'

        filecdf_read = Dataset(filecdf_path,'r')
#levels =  [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 600, 550, 500, 450, 400, 350, 300, 275, 250, 225, 200, 175, 150, 125, 100]
# Create fst file name
        filefst_path = re.sub('nc','fst',filecdf_path)
        iunit = Fstdc.fstouv(0,filefst_path,'RND+R/W')
        
        latitude=filecdf_read.variables['lat'][:]
        longitude=filecdf_read.variables['lon'][:]
        var=filecdf_read.variables['hgt'][:]
        lev=filecdf_read.variables['level'][:]
        temps=filecdf_read.variables['time'][:] 

        dt=len(temps)
        dz=len(lev)

        for t in range(0,dt): 
            for lv in range(0,dz): 
                temp = var[t,lv,:,:]   
                var2 = np.transpose(temp)
                tmp = np.array(var2).squeeze()
                tmp=np.where(abs(tmp)>1000000,float('nan'),tmp)
                IP1= lev[lv]
                IP2=(t+1)*3
# ier = Fstdc.fstecr(tmp,iunit,'GZ','P','FCST',IP1,IP2,IP3,DATEO,GridType,int(IG1),IG2,IG3,IG4,DEET,NPAS,16,1)     
                print(tmp.shape)
                ier = Fstdc.fstecr(tmp,iunit,VNAME,'P','FCST',IP1,IP2,IP3,DATEO,GridType,int(IG1),IG2,IG3,IG4,DEET,NPAS,32,1)

#print var_data3.shape, var_data2.shape                  
# Close file
        filecdf_read.close()
        ier = Fstdc.fstfrm(iunit)
