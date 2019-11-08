 
PRO statistique

;--------------------------------------------------------------------------
;                               
;  
;   
;       
;   Revise:par   Rabah Aider  ESCER/UQAM  07/05/2010
;   Adapte par E. D. POAN le 25 Fev 2015                                       
;
;                             
;--------------------------------------------------------------------------

; Ce programme calcule certaines caracteristiques des tempetes a
; partir des sorties de l<algorithme de Sinclair (2007)
;

; data_i:         Premier jour_heure d'analyse: jjhh = 100
; data_o:         Dernier jour_heure d'analyse: jjhh = 3118
; di:             Direction du déplacement de cyclones
; direc:          Option de calculer les provenances
; file_cont:      Fichier de sortie contenant les statistiques
; file_data:      Fichier d'entree (sortie du code de Sinclair)
; file_dir:       Fichier de sortie contenant la densité de cyclones selon la provenance
; file_mdir:      Fichier de sortie contenant l'intensité moyenne selon la provenance
; grille:         28 points x 144 points couvrent l'hemisphere de 20 degrees latitudes
;                (distance entre les points de grille est par 2.5 degrees)
; mi:             Intensité moyenne de cyclones en fonction da sa provenance
; month_i:        Premier mois de la période étudie; mm
; month_o:        Dernier mois de la période étudie; mm
; moy_cir:        Circulation totale 
; moy_vel:        Intensité totale 
; moy_wat:        Eau précipitable totale 
; st(0, *, *):    Vie moyenne
; st(1, *, *):    Nombre de cyclones
; st(2, *, *):    Nombre de genèses
; st(3, *, *):    Nombre de dissolutions
; st(4, *, *):    Nombre de trajectoires
; st(5, *, *):    Nombre de cyclones intenses
; st(6, *, *):    Vitesse moyenne
; st(7, *, *):    Circulation moyenne 
; st(8, *, *):    Intensité moyenne
; st(9, *, *):    Eau précipitable moyenne
; st(10, *, *):   Vitesse moyenne du vent 
; v(*,*,i):       Trajectoire courante
; v(*,j,*):       Position actuelle du cyclone sur la trajectoire
; v(0,*,*):       Latitude
; v(1,*,*):       Longitude
; v(2,*,*):       Eau précipitable
; v(3,*,*):       Tourbillon du vent de gradient
; v(4,*,*):       Circulation

;month=['01', '02', '03','11','12']
;month=['04', '05', '06', '07', '08', '09', '10']
;for imo=0,7  do begin
;for iyear=1981,2010 do begin

year=FIX(GETENV('ANNEE'))
month=GETENV('MOIS')
print, '---------------------------------------------------'
print, '---------------------------------------------------'
print, '-----------Traitement pour ',year, '--', month,'--------------'
print, '---------------------------------------------------'
for imo=0,0  do begin
for iyear=year,year do begin

NMAX=100000
n=0L
N1=0
N2=0
N3=0L
N4=0
nn=0L

ctrl=intarr(4,nmax)
dat_i=intarr(3118)
dat_o=intarr(3118)
month_i=intarr(12)
month_o=intarr(12)
ntj=intarr(4032)
di=intarr(14,91,361)
mq=fltarr(14,91,361)
mi=fltarr(14,91,361)
div=intarr(16,91,361)
mqv=fltarr(16,91,361)
mv=fltarr(16,91,361)
no=intarr(14)
nq=intarr(16)
nlat=29
nlon=145
lon=fltarr(nlon,nlat)
lat=fltarr(nlon,nlat)
g_lon=fltarr(nlon+164,nlat+6)
g_lat=fltarr(nlon+164,nlat+6)
st=fltarr(11,nlat+62,nlon+216)

str=lonarr(nlat+62,nlon+216)
moy_cir=lonarr(nlat+62,nlon+216)
moy_vel=lonarr(nlat+62,nlon+216)
moy_wat=lonarr(nlat+62,nlon+216)
moy_vit=lonarr(nlat+62,nlon+216)
dat_i=100
dat_o=3118
month_i=01
month_o=01
print, 'direc, 0=no and 1=yes'
direc=1
;fic_in='/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_OUT/CRCM5_MPI/BC_RCP85/out_crcm5_mpi_uv850_'
;fic_out='/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_STAT/CRCM5_MPI/BC_RCP85/stat_crcm5_mpi_uv850_'


fic_in='/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_OUT_GUI/NARR/out_narr_uv850_'
fic_out='/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_OUT_GUI/NARR/stat_narr_uv850_'




print,fic_in +strcompress(iyear,/remove_all)+ month(imo)  + '.txt'

openr,1,fic_in +strcompress(iyear,/remove_all)+ month(imo)  + '.txt'
openw,4,fic_out + month(imo) +strcompress(iyear,/remove_all) 


while not eof(1) DO begin
  n=n+1
  readf,1,format='(I4,i5,I6,I3)', N1,N2,N3,N4
  ctrl(0,n)=N1
  ctrl(1,n)=N2
  ctrl(2,n)=n
  ctrl(3,n)=N4

ENDWHILE


close,1

nn=n
d=fltarr(nn+1)
dmin=fltarr(nn+1)
ddd=fltarr(nn+100)
vd=fltarr(nn+1)
print,'nn',nn
nt=0L
ntrj=0L
ncyc=0L
qtrj=0L
qcyc=0L
ikod=0L
if n gt nmax then begin
 print, 'il faut changer nmax'
 stop
endif

ctrl=ctrl(*,0:n)

;print, 'ctrl == ', ctrl
nb_pts_max=max(ctrl(3,*))
print, 'nb_pts_max = ',nb_pts_max
v=fltarr(7,nb_pts_max,n+1)
dat=lonarr(nb_pts_max,n+1)
date=intarr(2,n+1)
date_years=fltarr(2,1)
day=lonarr(n+1)
openr,1,fic_in +strcompress(iyear,/remove_all)+ month(imo)  + '.txt'
lun1=1
for i=1L,nn do begin
 np=ctrl(3,i)
;print, 'ctrl == ', ctrl(*,i)


 nd=ctrl(2,i)
print, 'np == ',np
 vector=fltarr(7,np)
if (iyear lt 2000) then begin
print, iyear 
;readf,1,format='(i4,(i2,i2),I6,I3,1000(2I4,i6,i6,i6,i6,i6))',date_day, date_years, NO_TRAJ,NB_PTS,vector
readf,1,format='(i5,(i2,i2),I6,I3,1000(2I4,i6,i6,i6,i6,i6))',date_day, date_years, NO_TRAJ,NB_PTS,vector
endif else readf,1,format='(i4,(i2,i3),I6,I3,1000(2I4,i6,i6,i6,i6,i6))',date_day, date_years, NO_TRAJ,NB_PTS,vector
 
;print, 'vecteur == ',date_day, date_years, NO_TRAJ,NB_PTS,vector
 dat(0,i)=(date_years(0)*100+date_years(1))
 day(i)=date_day
 v(0,0:np-1,i)=reform(vector(0,*)/10.0)
 v(1,0:np-1,i)=reform(vector(1,*)/10.0)
;print, 'lat == ',v(0,0:np-1,i) , 'lon == ', v(1,0:np-1,i)
 v(2,0:np-1,i)=reform(vector(2,*))
 v(3,0:np-1,i)=reform(vector(3,*)/10.0)
 v(4,0:np-1,i)=reform(vector(4,*)/10)
 v(5,0:np-1,i)=reform(vector(5,*)/10.0)
 v(6,0:np-1,i)=reform(vector(6,*)/10.0)
 date(0,i)=reform(date_years(0,*))
 date(1,i)=reform(date_years(1,*))
 nt=nt+NB_PTS
endfor
close,1
print,iyear,'-',month(imo)
for i=1L,nn do begin
 d(i)=0
 dmin(i)=map_2points(v(1,0,i),v(0,0,i),v(1,ctrl(3,i)-1,i),v(0,ctrl(3,i)-1,i), /meters)/1000
 for j=0,ctrl(3,i)-2 do begin
   if v(1,j,i) ne v(1,j+1,i) or v(0,j,i) ne v(0,j+1,i) then begin
    diss=map_2points(v(1,j,i),v(0,j,i),v(1,j+1,i),v(0,j+1,i), /meters)/1000
    d(i)=d(i)+diss
   endif
endfor

endfor

for j=0,28 do begin
 for i=0,144 do begin
  lon(i,j)=2.5*i
  lat(i,j)=20+2.5*j
endfor
endfor

for j=3,nlat+2 do begin

 for i=82,nlon+81 do begin
  g_lon(i,j)=lon(i-82,j-3)
  g_lat(i,j)=lat(i-82,j-3)
 endfor

 for k=1,82 do begin
  g_lat(226+k,j)=lat(k,j-3)
  g_lon(226+k,j)=lon(k,j-3)
  g_lat(k-1,j)=lat(62+k,j-3)
  g_lon(k-1,j)=lon(62+k,j-3)
 endfor

endfor

for l=0,2 do begin
 for i=0,nlon-1 do begin
   g_lon(i+82,32+l)=lon(i,27-l)
   g_lat(i+82,32+l)=lat(i,27-l)
   g_lon(i+82,l)=lon(i,0)
   g_lat(i+82,l)=lat(i,0)
 endfor
 for k=0,81 do begin
   g_lat(k,l)=g_lat(k,3)
   g_lon(k,l)=g_lon(k,3)
   g_lat(k,32+l)=g_lat(144+k,32+l)
   g_lon(k,32+l)=g_lon(144+k,32+l)
   g_lat(226+k,32+l)=g_lat(82+k,32+l)
   g_lon(226+k,32+l)=g_lon(82+k,32+l)
   g_lat(226+k,l)=g_lat(226+k,3)
   g_lon(226+k,l)=g_lon(226+k,3)
 endfor
endfor

for l=82,nlon+81 do begin
 for s=3,31 do begin
  for k=0,9 do begin
   st(k,g_lat(l,s),g_lon(l,s))=0L
  endfor
  for ss=0,7 do begin
   di(ss,g_lat(l,s),g_lon(l,s))=0L
   mq(ss,g_lat(l,s),g_lon(l,s))=0L
 endfor
for is=0,15 do begin
   div(is,g_lat(l,s),g_lon(l,s))=0L
   mqv(is,g_lat(l,s),g_lon(l,s))=0L
;   ntj(is,g_lat(l,s),g_lon(l,s))=0L
 endfor

for ii=0,nn do begin
ntj(ii)=0
endfor

ikod=0L

 moy_cir(g_lat(l,s),g_lon(l,s))=0L
 moy_vel(g_lat(l,s),g_lon(l,s))=0L
 moy_wat(g_lat(l,s),g_lon(l,s))=0L
 moy_vit(g_lat(l,s),g_lon(l,s))=0L
 endfor
endfor

M=[0,0]
N=[0,3]
Q=[0,7.5]
dcrit_3=map_2points(M[0],M[1],N[0],N[1], /meters)

dcrit_7=map_2points(M[0],M[1],Q[0],Q[1], /meters)
jcod=0L
  for s=3,31 do begin
  for l=82,226 do begin
  ;    if g_lat(l,s) ge 50 and g_lat(l,s) le 70 then begin
  ;    if g_lon(l,s) ge 265 and g_lon(l,s) le 295 then begin	


    d0=map_2points(g_lon(l,s),g_lat(l,s),g_lon(l+1,s),g_lat(l+1,s),/meters)
    if d0 gt 0 then distx=dcrit_7/d0
    distx=uint(distx)
    i1=g_lat(l,s-3)
    i2=g_lat(l,s+3)
    j1=g_lon(l-distx,s)
    j2=g_lon(l+distx+1,s)
    error0=0
    error1=0
    if j1 gt j2 then error1=1
    if s ge 29 then error0=1
    if s eq 31 and l ge 83 then error0=2
    cont=0
    cont_q=0
    for ll=0,13 do begin
     no(ll)=0
    endfor
    for is=1,15 do begin
    nq(is)=0
    endfor
   

  for i=1L,nn do begin
      
     if d(i) ge 1200 and dmin(i) ge 600 then begin
      bingo=0
      bingo5=0
      bingof=0
      np=ctrl(3,i)
      if error0 eq 0 then begin
       if i1 le v(0,0,i) and v(0,0,i) le i2 then begin
        if error1 eq 0 then begin
         if j1 le v(1,0,i) and v(1,0,i) le j2 then goto, cal_g
        endif
        if error1 eq 1 then begin
         if j1 le v(1,0,i) or v(1,0,i) le j2 then goto, cal_g
        endif
       endif
      endif 
      if error0 eq 1 then begin
       if i1 le v(0,0,i) then goto, cal_g
      endif
cal_g: if g_lon(l,s) ne v(1,0,i) or g_lat(l,s) ne v(0,0,i) then begin
        d_gen=map_2points(g_lon(l,s),g_lat(l,s),v(1,0,i),v(0,0,i),/meters)
        if d_gen le dcrit_3 then begin
         st(2,g_lat(l,s),g_lon(l,s))=st(2,g_lat(l,s),g_lon(l,s))+1
         if ctrl(0,i) eq dat_i and date(1,i) eq month_i then st(2,g_lat(l,s),g_lon(l,s))=st(2,g_lat(l,s),g_lon(l,s))-1
        endif
       endif
      dd=0
      bau=0
      for j=0,np-1 do begin
      if error0 eq 0 then begin
       if i1 le v(0,j,i) and v(0,j,i) le i2 then begin
        if error1 eq 0 then begin
         if j1 le v(1,j,i) and v(1,j,i) le j2 then goto, cal_c
        endif
        if error1 eq 1 then begin
         if j1 le v(1,j,i) or v(1,j,i) le j2 then goto, cal_c
        endif
       endif
      endif 
      if error0 eq 1 then begin
       if i1 le v(0,j,i) then goto, cal_c
      endif
 cal_c: dista=map_2points(g_lon(l,s),g_lat(l,s),v(1,j,i),v(0,j,i),/meters)

        if j lt np-1 then begin
         dista_n=map_2points(g_lon(l,s),g_lat(l,s),v(1,j+1,i),v(0,j+1,i),/meters)
     	 if v(1,j,i) ne v(1,j+1,i) or v(0,j,i) ne v(0,j+1,i) then begin
	  if dista le dcrit_3 then begin
	   st(1,g_lat(l,s),g_lon(l,s))=st(1,g_lat(l,s),g_lon(l,s))+1
 	   moy_cir(g_lat(l,s),g_lon(l,s))=moy_cir(g_lat(l,s),g_lon(l,s))+v(4,j,i)
 	   moy_vel(g_lat(l,s),g_lon(l,s))=moy_vel(g_lat(l,s),g_lon(l,s))+v(3,j,i)
 	   moy_wat(g_lat(l,s),g_lon(l,s))=moy_wat(g_lat(l,s),g_lon(l,s))+v(2,j,i)
           moy_vit(g_lat(l,s),g_lon(l,s))=moy_vit(g_lat(l,s),g_lon(l,s))+v(5,j,i)
	   if bingo eq 0 then bingo=1
           bingo5=bingo
 	    ddd(j)=map_2points(v(1,j,i),v(0,j,i),v(1,j+1,i),v(0,j+1,i),/meters)/1000
 	    dd=dd+ddd(j)
	    bau=bau+6
;  	    if v(3,j,i) ge 6 then st(5,g_lat(l,s),g_lon(l,s))=st(5,g_lat(l,s),g_lon(l,s))+1
       	   endif
           if bingo eq 0 then begin
       	    if dista le dcrit_7 and dista_n le dcrit_7 then begin
 	     if v(1,j,i) ge 334 and v(1,j+1,i) le 20 and g_lon(l,s) le 334 then v(1,j,i)=v(1,j,i)-360
 	     if v(1,j,i) ge 334 and v(1,j+1,i) le 20 and g_lon(l,s) ge 334 then v(1,j+1,i)=v(1,j+1,i)+360
 	     try=[PNT_LINE([g_lon(l,s),g_lat(l,s)],[v(1,j,i),v(0,j,i)],[v(1,j+1,i),v(0,j+1,i)],Pl),pl]
 	     m0=min([v(1,j,i),v(1,j+1,i)])
 	     mx0=max([v(1,j,i),v(1,j+1,i)])
 	     m1=min([v(0,j,i),v(0,j+1,i)])
  	     mx1=max([v(0,j,i),v(0,j+1,i)])
  	     if m0 gt 360 then try(1)=try(1)+360
 	     if m0 le try(1) and mx0 ge try(1) then begin
 	      if m1 le try(2) and mx1 ge try(2) then begin
               d_try=map_2points(g_lon(l,s),g_lat(l,s),try(1),try(2),/meters)
               if d_try le dcrit_3 then bingo=1
 	      endif
 	    endif
 	   endif
 	  endif
 	 endif
        endif
	if direc eq 1 then begin
            
        if j ge 1 then begin
         if v(1,j-1,i) ne v(1,j,i) or v(0,j-1,i) ne v(0,j,i) then begin
 	  if j eq 1 then dir=map_2points(v(1,j-1,i),v(0,j-1,i),v(1,j,i),v(0,j,i))
          if j eq 2 then dir=map_2points(v(1,j-2,i),v(0,j-2,i),v(1,j,i),v(0,j,i))
          if j eq 3 then dir=map_2points(v(1,j-3,i),v(0,j-3,i),v(1,j,i),v(0,j,i))
	  if j eq 4 then dir=map_2points(v(1,j-4,i),v(0,j-4,i),v(1,j,i),v(0,j,i))
;	  if j eq 4 then dir1=dir
;	  if j ge 5 then dir=dir1
          if j ge 4 then dir=map_2points(v(1,j-4,i),v(0,j-4,i),v(1,j,i),v(0,j,i))
; 	 	print,dir  
 	  if dista le dcrit_3 then begin
		
	 for ss=-5,7,2 do begin
	
             if dir(1) lt ss*22.5 and dir(1) ge (ss-2)*22.5 then begin
	     
 	     ll=ss+5		   
	   di(ll,g_lat(l,s),g_lon(l,s))=di(ll,g_lat(l,s),g_lon(l,s))+1		
           mq(ll,g_lat(l,s),g_lon(l,s))=mq(ll,g_lat(l,s),g_lon(l,s))+v(3,j,i)            
	     no(ll)=no(ll)+1	     
	     endif 	   
  	  endfor
		
 	      if dir(1) lt -7.*22.5 and dir(1) ge 7.*22.5 then begin	     
	    
	    	di(13,g_lat(l,s),g_lon(l,s))=di(13,g_lat(l,s),g_lon(l,s))+1	

             mq(13,g_lat(l,s),g_lon(l,s))=mq(13,g_lat(l,s),g_lon(l,s))+v(3,j,i)	    
 	     no(13)=no(13)+1
 	     endif
	    
             dirv=v(6,j,i)
	     
             for is=1,13,2 do begin
             if dirv gt is*22.5 and dirv le (is+2)*22.5 then begin
             div(is,g_lat(l,s),g_lon(l,s))=div(is,g_lat(l,s),g_lon(l,s))+1
 	     mqv(is,g_lat(l,s),g_lon(l,s))=mqv(is,g_lat(l,s),g_lon(l,s))+v(5,j,i)    
             nq(is)=nq(is)+1
             endif
             endfor

             if dirv gt 15*22.5 or dirv le 22.5 then begin
             div(15,g_lat(l,s),g_lon(l,s))=div(15,g_lat(l,s),g_lon(l,s))+1
 	     mqv(15,g_lat(l,s),g_lon(l,s))=mqv(15,g_lat(l,s),g_lon(l,s))+v(5,j,i)   
             nq(15)=nq(15)+1
             endif

           endif
         endif
        endif
	endif
      endfor

	if dista le dcrit_3 then ikod=ikod+1

      if bau ne 0 then vd(i)=float(dd)/float(bau)
      if error0 eq 0 then begin
       if i1 le v(0,np-1,i) and v(0,np-1,i) le i2 then begin
        if error1 eq 0 then begin
         if j1 le v(1,np-1,i) and v(1,np-1,i) le j2 then goto, cal_d
        endif
        if error1 eq 1 then begin
         if j1 le v(1,np-1,i) or v(1,np-1,i) le j2 then goto, cal_d
        endif
       endif
       endif 
       if error0 eq 1 then begin
        if i1 le v(0,np-1,i) then goto, cal_d
       endif
cal_d:if g_lon(l,s) ne v(1,(np-1),i) or g_lat(l,s) ne v(0,(np-1),i) then begin
       d_dis=map_2points(g_lon(l,s),g_lat(l,s),v(1,(np-1),i),v(0,(np-1),i),/meters)
       if d_dis le dcrit_3 then begin 
           st(3,g_lat(l,s),g_lon(l,s))=st(3,g_lat(l,s),g_lon(l,s))+1
           if ctrl(0,i) eq dat_o and date(1,i) eq month_o then st(3,g_lat(l,s),g_lon(l,s))=st(3,g_lat(l,s),g_lon(l,s))-1
       endif
       if v(1,np-2,i) ne v(1,np-1,i) or v(0,np-2,i) ne v(0,np-1,i) then begin
        if d_dis le dcrit_3 then begin
         st(1,g_lat(l,s),g_lon(l,s))=st(1,g_lat(l,s),g_lon(l,s))+1
         bingof=1
         moy_cir(g_lat(l,s),g_lon(l,s))=moy_cir(g_lat(l,s),g_lon(l,s))+v(4,np-1,i)
         moy_vel(g_lat(l,s),g_lon(l,s))=moy_vel(g_lat(l,s),g_lon(l,s))+v(3,np-1,i)
         moy_wat(g_lat(l,s),g_lon(l,s))=moy_wat(g_lat(l,s),g_lon(l,s))+v(2,np-1,i)
         moy_vit(g_lat(l,s),g_lon(l,s))=moy_vit(g_lat(l,s),g_lon(l,s))+v(5,np-1,i)
;         if v(3,np-1,i) ge 6 then st(5,g_lat(l,s),g_lon(l,s))=st(5,g_lat(l,s),g_lon(l,s))+1
        endif
       endif
   endif

  if bingo5 eq 1 then begin
   st(6,g_lat(l,s),g_lon(l,s))=st(6,g_lat(l,s),g_lon(l,s))+vd(i)
   str(g_lat(l,s),g_lon(l,s))=str(g_lat(l,s),g_lon(l,s))+1
  endif
  if bingo eq 1 or bingof eq 1 then st(4,g_lat(l,s),g_lon(l,s))=st(4,g_lat(l,s),g_lon(l,s))+1
   if bingo eq 1 or bingof eq 1 then if v(3,np-1,i) ge 6 then st(5,g_lat(l,s),g_lon(l,s))=st(5,g_lat(l,s),g_lon(l,s))+1
  if bingo eq 1 then st(0,g_lat(l,s),g_lon(l,s))=st(0,g_lat(l,s),g_lon(l,s))+np


  if bingo ne 1 and bingof eq 1 then st(0,g_lat(l,s),g_lon(l,s))=st(0,g_lat(l,s),g_lon(l,s))+np
  endif
  endfor

  if str(g_lat(l,s),g_lon(l,s)) gt 0 then st(6,g_lat(l,s),g_lon(l,s))=st(6,g_lat(l,s),g_lon(l,s))/str(g_lat(l,s),g_lon(l,s))
  
if st(4,g_lat(l,s),g_lon(l,s)) gt 0 then begin
   num=float(4)*float(st(4,g_lat(l,s),g_lon(l,s)))
   st(0,g_lat(l,s),g_lon(l,s))=10*st(0,g_lat(l,s),g_lon(l,s))/float(num)
   st(0,g_lat(l,s),g_lon(l,s))=uint(st(0,g_lat(l,s),g_lon(l,s)))
  endif
   if direc eq 1 then begin
    for ss=0,13 do begin
     if no(ss) gt 0 then mi(ss,g_lat(l,s),g_lon(l,s))=float(mq(ss,g_lat(l,s),g_lon(l,s)))/float(no(ss))
     if no(ss) eq 0 then mi(ss,g_lat(l,s),g_lon(l,s))=0
	
    endfor

    for is=1,15 do begin
     if nq(is) gt 0 then mv(is,g_lat(l,s),g_lon(l,s))=float(mqv(is,g_lat(l,s),g_lon(l,s)))/float(nq(is))
     if nq(is) eq 0 then mv(is,g_lat(l,s),g_lon(l,s))=0
    endfor
   endif 
  if st(1,g_lat(l,s),g_lon(l,s)) gt 0 then begin
    st(7,g_lat(l,s),g_lon(l,s))=float(moy_cir(g_lat(l,s),g_lon(l,s)))/float(st(1,g_lat(l,s),g_lon(l,s)))
    st(8,g_lat(l,s),g_lon(l,s))=float(moy_vel(g_lat(l,s),g_lon(l,s)))/float(st(1,g_lat(l,s),g_lon(l,s)))
    st(9,g_lat(l,s),g_lon(l,s))=float(moy_wat(g_lat(l,s),g_lon(l,s)))/float(st(1,g_lat(l,s),g_lon(l,s)))
    st(10,g_lat(l,s),g_lon(l,s))=float(moy_vit(g_lat(l,s),g_lon(l,s)))/float(st(1,g_lat(l,s),g_lon(l,s)))
  endif
  if error0 eq 2 then goto, calcul

endfor

endfor

calcul: for s=83,226 do begin
	 for i=0,9 do begin
	  st(i,90,g_lon(s,31))=st(i,90,0)
	 endfor
        endfor

finalul=0L
fing=0L
find=0L
for s=3,31 do begin
  for l=82,226 do begin
    if l lt 226 and s lt 31 then finalul=finalul+st(1,g_lat(l,s),g_lon(l,s))
    if l lt 226 and s lt 31 then fing=fing+st(2,g_lat(l,s),g_lon(l,s))
    if l lt 226 and s lt 31 then find=find+st(3,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then finalul=finalul+st(1,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then  fing=fing+st(2,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then find=find+st(3,g_lat(l,s),g_lon(l,s))
    printf,4,format='(f5.1,f8.1,i6,i6,i6,i6,i6,i6,i6,f8.1,f8.1,f8.1)',$
 	g_lat(l,s),g_lon(l,s),st(2,g_lat(l,s),g_lon(l,s)),st(3,g_lat(l,s),g_lon(l,s)),$
         st(4,g_lat(l,s),g_lon(l,s)),st(5,g_lat(l,s),g_lon(l,s)),st(0,g_lat(l,s),g_lon(l,s)),st(6,g_lat(l,s),g_lon(l,s)),$
        st(7,g_lat(l,s),g_lon(l,s)),st(8,g_lat(l,s),g_lon(l,s)),st(9,g_lat(l,s),g_lon(l,s)),st(10,g_lat(l,s),g_lon(l,s))
endfor
endfor
for s=3,30 do begin
  for l=82,226 do begin
    if l lt 226 and s lt 31 then finalul=finalul+st(1,g_lat(l,s),g_lon(l,s))
    if l lt 226 and s lt 31 then fing=fing+st(2,g_lat(l,s),g_lon(l,s))
    if l lt 226 and s lt 31 then find=find+st(3,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then finalul=finalul+st(1,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then  fing=fing+st(2,g_lat(l,s),g_lon(l,s))
    if s eq 31 and l eq 82 then find=find+st(3,g_lat(l,s),g_lon(l,s))
     	
     if direc eq 1 then begin 
    ;  printf,5,format='(f5.1,f8.1,i6,i6,i6,i6,i6,i6,i6,i6)',$  ; E,SE,S,SO,O,NO,N,NE
    ;    g_lat(l,s),g_lon(l,s),di(2,g_lat(l,s),g_lon(l,s)),di(4,g_lat(l,s),g_lon(l,s)),$
    ;    di(6,g_lat(l,s),g_lon(l,s)),di(8,g_lat(l,s),g_lon(l,s)),$
    ;   di(10,g_lat(l,s),g_lon(l,s)),di(12,g_lat(l,s),g_lon(l,s)),$
    ;   di(13,g_lat(l,s),g_lon(l,s)),di(0,g_lat(l,s),g_lon(l,s))
     
    ;    printf,2,format='(f5.1,f8.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1)',$
    ;     g_lat(l,s),g_lon(l,s),mi(2,g_lat(l,s),g_lon(l,s)),mi(4,g_lat(l,s),g_lon(l,s)),mi(6,g_lat(l,s),g_lon(l,s)),$
     ;    mi(8,g_lat(l,s),g_lon(l,s)),mi(10,g_lat(l,s),g_lon(l,s)),mi(12,g_lat(l,s),g_lon(l,s)),mi(13,g_lat(l,s),g_lon(l,s)),$
      ;   mi(0,g_lat(l,s),g_lon(l,s))
      
      ;   printf,3,format='(f5.1,f8.1,i6,i6,i6,i6,i6,i6,i6,i6)',$  ; E,SE,S,SO,O,NO,NE
      ;  g_lat(l,s),g_lon(l,s),div(7,g_lat(l,s),g_lon(l,s)),div(5,g_lat(l,s),g_lon(l,s)),div(3,g_lat(l,s),g_lon(l,s)),$
       ; div(1,g_lat(l,s),g_lon(l,s)),div(15,g_lat(l,s),g_lon(l,s)),div(13,g_lat(l,s),g_lon(l,s)),div(11,g_lat(l,s),g_lon(l,s)),$
      ;  div(9,g_lat(l,s),g_lon(l,s))
      ;    printf,7,format='(f5.1,f8.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1,f5.1)',$
      ;   g_lat(l,s),g_lon(l,s),mv(7,g_lat(l,s),g_lon(l,s)),mv(5,g_lat(l,s),g_lon(l,s)),mv(3,g_lat(l,s),g_lon(l,s)),$
       ;  mv(1,g_lat(l,s),g_lon(l,s)),mv(15,g_lat(l,s),g_lon(l,s)),mv(13,g_lat(l,s),g_lon(l,s)),mv(11,g_lat(l,s),g_lon(l,s)),$
       ; mv(9,g_lat(l,s),g_lon(l,s))


    endif
  endfor
endfor

print,ntrj,nn,ncyc,nt
print,'Resultat', finalul,fing,find

;for i=1,7 do begin
; close,i
;endfor

close, 1
close, 4
print, 'Nombre total des trajectoires est:',nn
print, 'Nombre total des cyclones est:',nt

print, 'bingo'
close,lun1
endfor
endfor
end
