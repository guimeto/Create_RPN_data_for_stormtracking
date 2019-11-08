	PROGRAM MAKE_TRACK
*
c**********************************************************************
c
c       ---------------------------------------------------------------
c       Programme pour determiner de facon objective
c       les trajectoires des depressions et anticyclones
c       ---------------------------------------------------------------
c       Version originale (1.0): Mark Sinclair (ERAU)
c       Reference:
c       - Mark R. Sinclair, 1997: "Objective Identification of Cyclones 
c       and Their Circulation Intensity, and Climatology", Weather and
c       Forecasting, 12, 3, 595-612. 
c       ---------------------------------------------------------------
c       Adaptation au format RPN: Corina Rosu Costea (SCA/UQAM)
c       ---------------------------------------------------------------
c
c**********************************************************************
c
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c       > Version 1.3 - JF Caron, ARMA/MRD/EC, Printemps-Ete 2007	
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c       -----------------------------------------
c       *     Nouveautees de cette version      *
c       -----------------------------------------
c       - Ajout de plusieurs sorties en fichier STD
c       - Filtrage deplace APRES les calculs;
c       - Cles d'appel pour fichiers input/output;
c       - Cles d'appel pour les options principales;
c       - Nouvelle S-R (ecrit_fst) pour ecriture en fichier STD
c       - Regroupement des S-R dans la librairie LIBTRACKS
c       -----------------------------------------
*
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::       
c       > Version 1.2 - Milka Radojevic, SCA/UQAM, Octobre 2004
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c       Reference: 
c       "Activité des cyclones extra-tropicaux simulés par le modèle 
c       couplé canadian de circulation générale (MCCG3)" - Mémoire de 
c       maitrise, UQAM, 2006.
c       contact: milka.radojevic@ec.gc.ca ou radojevic.milka@ouranos.ca
c       -----------------------------------------
c       *     Nouveautees de cette version      *
c       -----------------------------------------
c       - Adaptation du code pour lire plusieurs fichiers à la fois et
c         créer autant de fichiers de sorties;		   
c       - Correction du calcul de la circulation;
c       -----------------------------------------
*
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c       > Version 1.1 - Corina Rosu Costea, SCA/UQAM, Nov 2003-Avr 2004
c       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c       Reference: 
c       "Les caractéristiques des cyclones et l`apport d`eau dans les 
c       bassins versants du Québec" - Mémoire de maitrise, UQAM, 2005.
c       contact: rosu@sca.uqam.ca ou rosu.corina@ouranos.ca
c       -----------------------------------------
c       *     Nouveautees de cette version      *
c       -----------------------------------------
c      - Augmentation de dist_crit de 4.5 à 7 pour 4 analyses/jour 
c        et de 2.25 à 3.75 pour 8 analyses/jour; 
c      - Augmentation de NG a 9 et changement nécessaire dans le vecteur iorder
c      - Introduction du champ ZS dans la S-R read_orog; 
c      - Correction du calcul du point estimé;
c      - A la suggestion de Peter Zwack, on prend la moitié du vent à 500 mb 
c        au lieu du vent climatologique pour estimer la vitesse de deplacement;
c      - Introduction de tf pour inclure toutes les trajectoires qui ne sont 
c        pas fini au dernier pas de temps du fichier d`entre;
c      - Création ecrit*.f  pour visualiser le champ du tourbillon calculé
c        par l`algorithme et introduction des paramètres pastime, npasc et 
c        dateo;
c      - Correction du calcul du tourbillon par l`introduction de la force
c        gravitationnelle, g
c       -----------------------------------------
*
c	> Version 1.0 - Août 2003
c       -----------------------------------------
c	- Reception du code de Sinclair
c       -----------------------------------------
*
c**********************************************************************
c
c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c----------------------------------------------------------------------
c       Definitions des variables:
c----------------------------------------------------------------------
c       SRAD      = Cressman smoothing radius (km)
c       naj       = nombre de donnees par jour (24 hrs)
c	npt       = temps de vie minimum (en pas de temps) d'une trajectoire
c	dist_crit = diametre du cercle de recherche (en deg)
c       nggg      = ???
c       N_INCR    = time increment of your data (h)
c----------------------------------------------------------------------
c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c      
	PARAMETER(IJD=40000,NTW=2000,NT=2000,NC=2000)
	INTEGER*4 TIME(NTW)
	CHARACTER*72 BASE,FCHAR
	CHARACTER*10 NUMB
	integer npas, naj, npt, dateo, cressman
	integer level,filtre,isign2,prcnt,pnmi,anal,nofst
	integer n_moins, n_fev, tf, ndz, nrl, nggg
	real work1(ijd),work2(ijd),work3(ijd)
	INTEGER*4 YEAR(NTW),MONTH(NTW),DAY(NTW),HOUR(NTW),NNT(NC),NTIME(
     %NT,NC)
	integer*4 directie1(ntw),lungime1(ntw)
	integer*4 NG(NT),ICM(NC),ITM(NC),ICG(NC,NT),ITG(NC,NT),INTS(4)
	integer*4 IORDER(9,362880),ICL(20)
	integer*4 pression(ntw),pression1(ntw),tourbillon(ntw),presiune(
     %ntw), presiune1(ntw)
	real dist_crit,vc
	REAL*4 CLAT(NT,NC),CLON(NT,NC),CVORT(NT,NC),CPRESS(NT,NC),CU(NT,
     %NC)
        REAL*4,CWIND(NT,NC),CDIRC(NT,NC),RLON(NC),RLAT(NC)
	real*4 CV(NT,NC),CCIRC(NT,NC),TVORT(NC),TPRESS(NC),TLAT(NC)
	real*4 TLON(NC),TCIRC(NC),TKE(NC),TU(NC),TV(NC),PLAT(NC),PLON(NC
     %)
        real*4 TU1(NC),TV1(NC),DIRC(NC),DIRC1(NC),WIND(NC),PHI(NC),BETA(
     %NC),TAN(NC)
	real*4 TUg(NC),TVg(NC)
        real*4 RPRES(NC),RDIR(NC),RDIR1(NC),RWIND(NC)
 	REAL*4 RLONN(NT,NC),RLATT(NT,NC),RWINDD(NT,NC),RDIRR(NT,NC),RPRESS(NT,
     %NC)
	real*4 PMN(NC,NC),BLAT(NT,NC),BLON(NT,NC),XLAT(NC),XLON(NC),RVOR
     %T(NC),RU1(NC),RV1(NC),COR(NT,NC)
	real*4 CLAT1(300,300,1200),CLON1(150,150,1200),CPRESS1(150,150,1
     %200)
	real*4 FLON(50,50,70),FLAT(50,50,70),FLON1(50,50,70),FLAT1(50,50
     %,70)
	real*4 CWIND1(150,150,1200),CDIRC1(150,150,1200)
	real*4 CLAT2(150,1200),CLON2(150,1200),CPRESS2(150,1200),CWIND2(
     %150,1200),CDIRC2(150,1200)	
	REAL*4 SLAT(IJD),SLON(IJD),H(IJD),VORT(IJD),PRESS(IJD)
	real*4 U(IJD),V(IJD),U1(IJD),V1(IJD),UG(IJD),VG(IJD),FB(IJD),FC(
     %IJD),OROG(IJD)
	real*4 gg(ijd),vct(ijd),pn(ijd),ctg(ijd),term(ijd),high_terrain(
     %ijd)
	real pressnf(ijd),hnf(ijd)
	LOGICAL FOUND,TUSED(NC),CUSED(NC),VORT_C,IJFOUND,INS,DO_CIRC
	logical TCLOSED(NC),CLOSED(NT,NC),ALL_THERE,BAD
	INTEGER*2 YYMM,NTR,DDHH,ILAT(NC),ILON(NC),IPRESS(NC),IVORT(NC)
	integer*2 ICIRC(NC),IDIRC(NC),IWIND(NC),IDIAM(NC),NUM1(NC),NUM2(
     %NC),NUM3(NC,NC),NM1(NC),NUM4(NC,NC),NM2(NC)
c
        PARAMETER(NCC=80000,LX=200,LY=200)	
	COMMON/GRID/NI,NJ,MAP,PMAP(10)
	COMMON/NEWGRID/NX,NY
c
        REAL*4 PARAMS(50),PMAP_WF(10),C_VORT(NCC),C_PRESS(NCC)
	real*4 PDS_VORT(6),PDS_PRESS(6),X(LX),Y(LY),XDRAW(100),YDRAW(100
     %)
c
	CHARACTER*3 MNTH(12)
	CHARACTER*10 SEASON
	CHARACTER*1 ANS
	CHARACTER*256 WF_NAME,TXT_FILE, name
	CHARACTER*256 infile, outfile
	character*2 getvar
	character*1 typvar, grtyp
c
	integer ier, ip1, ip3, silent, hi, hf, ip2
	integer fnom, fstouv, fstopc, fstfrm, fclos
	integer real_len
	integer iun01, iun51, glbkey, get_dimgrille, fstprm_dio
	integer nstamp(ntw)
	integer nrecs(3)
	integer deet,ig1,ig2,ig3,ig4
	DATA RAD/6370000./,PI/3.14159265/,R/287./,OME/7.292E-5/
c       NI, NJ, MAP, PMAP(1)...PMAP(10) define the map projection parameters 
c       for the computational grid. This is different from the data grid
c       The computational grid, defined in the next line, is a polar 
c       stereographic 61 x 61 grid covering the hemisphere with 180 long 
c       pointing in positive y-direction. 
c	DATA NI/61/,NJ/61/,MAP/0/,PMAP/31.,31.,180.,43.2,-60.,0.,0.,0.,0
c       .,0./
c	initialized later
c
	DATA MNTH/'JAN','FEB','MAR','APR','MAY','JUN',
     -            'JUL','AUG','SEP','OCT','NOV','DEC'/
C
c	LOGICAL NCEP
	RD=PI/180.
C
 100    FORMAT(A)
c**********************
c       Traitement des cles dappel
c**********************
	call get_options(silent,ip3,getvar,infile,outfile,naj,npt,
	1    dist_crit,level,isign2,filtre,prcnt,pnmi,vc,anal,nofst)
       print *, "okkkkkkkkkkkkkkkkkkkkkk1"	
c**********************
c       Initialisations
c**********************
	iun01 = 1
	iun51 = 51
	WF_NAME = infile(1:real_len(infile)) 
        print*
	name = outfile(1:real_len(outfile))//'.rpn' 
	TXT_FILE = outfile(1:real_len(outfile))//'.txt' 
c        TXT_FILE2 = outfile(2:real_len(outfile))//'.txt' ! fichier sortie TXT
c***  Ouverture fichier input RPN
       print *, "okkkkkkkkkkkkkkkkkkkkkk2"	
      ier = fnom(iun01, WF_NAME, 'RND', 0)
       print *, "okkkkkkkkkkkkkkkkkkkkkk3"	
      nrecs(1) = fstouv(iun01, 'RND')
       print *, "okkkkkkkkkkkkkkkkkkkkkk3"	
	if (nofst.ne.0) then
c***    Ouverture fichier output RPN
	   ier = fnom(iun51, name, 'RND', 0)
	   nrecs(2) = fstouv(iun51, 'RND')
	endif
c***  On ne tolere pas les erreurs a partir de warning et plus
      ier = fstopc('TOLRNC','WARNIN',.false.)
c***  Obtenir dimension grille, pas de temps, niveaux verticaux
	glbkey = get_dimgrille(ni_wf, nj_wf, iun01, getvar)
	if ( (ni_wf*nj_wf).gt.ijd) then
	   print*
	   print*,'Error : size of input grid greater than IJD'
	   print*,'Stopping'
	   print*
	   stop
	endif
	call get_proj(iun01, glbkey, map_wf, params)
	DO K=1,10
	   PMAP_WF(K)=PARAMS(K)
	END DO
c	Dans cette version, la grille de travail (ni,nj) est la meme 
c       que la grille des donnees en input (ni_wf,nj_wf)
	ni = ni_wf
	nj = nj_wf
        print*, 'ni_wf = ', ni
        print*, 'nj_wf = ', nj
	NIJ=NI*NJ
	map = map_wf
	do K = 1,10
	   pmap(K) = pmap_wf(K)
	enddo
c Find the mid-point of the data file domain and see if it is NH or SH
	CALL LL(MAP_WF,PMAP_WF,0.5*NI_WF,0.5*NJ_WF,SLAT_M,SLON_M)
	IF(SLAT_M.GT.0.0)THEN
	   IHEM=1
	   print*
	   PRINT *,'Processing for Northern Hemisphere ...'
	   print*
	ELSE
	   IHEM=-1
	   print*
	   PRINT *,'Processing for Southern Hemisphere ...'
	   print*
	ENDIF
c***
      ier=fstprm_dio(dateo,deet,typvar,grtyp,ig1,ig2,ig3,ig4,glbkey)
c***************************
c       Costal Erosion: pour l`étude des  érossions côtières il 
c                       faut mettre erosion=1
c***************************
        print*
c 	 write (*,*)'Costal Erosion Study, erosion=',erosion 
c        read (*,*) erosion 
	 erosion=1
cc***************************
c       Wind units
c***************************
c        print*
c 	 write (*,*)'Speed wind units : 1 for kt', unit
c        read (*,*) unit ! 1 pour noueds
c	 unit=1
***************************
c       Cyclones or Anticyclones ?
c***************************
	if (isign2.eq.1) then
  	   ISIGN=1
	   print*
 	   PRINT*,'Tracking is for CYCLONES'
	   print*
 	ELSE IF(isign2.eq.-1)THEN
 	   ISIGN=-1
	   print*
 	   PRINT*,'Tracking is for ANTICYCLONES'
	   print*
  	ELSE
 	   PRINT *,'ISIGN2: Invalid option'
 	   stop
 	ENDIF
c***************************
c       Track vorticity or pressure centres [V/P] ?
c***************************
	print*,'prcnt',prcnt
	IF(prcnt.eq.-1)THEN
	   VORT_C=.TRUE.
	   PRINT *,'Tracking VORT extrema'
	   DO_CIRC=.TRUE.
	ELSE IF(prcnt.eq.0)THEN
	   VORT_C=.FALSE.
	   PRINT *,'Tracking PRESSURE extrema'
	   DO_CIRC=.FALSE.
	ELSE
	   PRINT *,'PRCNT: Invalid answer'
	   stop
	ENDIF
c***************************
c       Smoothing parameters
c***************************
c       select smoothing radius
c	Ex: SRAD=800 !Cressman smoothing within a radius of 800. km
	SRAD=float(FILTRE)*1.E3
c***************************
c       Work out latitude and Coriolis parameters
c***************************
c       save X and Y index arrays
	DO J=1,NJ
  	   Y(J)=FLOAT(J)
	END DO
	DO I=1,NI
  	   X(I)=FLOAT(I)
	END DO
c       work out latitude and Coriolis parameters
c       here LL gets a lat,lon from index X,Y (X:1,NI, Y:1,NJ)
	DO J=1,NJ
	DO I=1,NI
	   IJ=I+(J-1)*NI
	   CALL LL(MAP,PMAP,X(I),Y(J),SLAT(IJ),SLON(IJ))
	   SLATT=SLAT(IJ)
	   IF(ABS(SLATT).LT.20.) THEN
	     IF(SLAT(IJ).LE.0.)SLATT=-30.
	     IF(SLAT(IJ).GT.0.)SLATT=30.
	   ENDIF
	   FC(IJ)=2.*OME*SIND(SLATT)	
	END DO
	END DO
       print *, "okkkkkkkkkkkkkkkkkkkkkk3"	
c**********************************************************************
c       Read orography
c**********************************************************************
	print*
	print*,'Read orography ...'
	print*
	call read_orog(orog,iun01,ijd)
        print*,'--------------------------------------------------------
     %----------------------------------'
        print*,'OROG MIN == ',minval(orog)
        print*,'OROG MAX == ',maxval(orog)
        print*,'--------------------------------------------------------
     %----------------------------------'
	nofst=0
	if (nofst.ne.0) then
	   call ecrit_fst(iun51,orog,'MX','TOPO_NF ',0,0,
	1	0,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	endif
c       smooth orography to 400. km
	CALL SMOOTH_C (OROG,NI_WF,NJ_WF,MAP_WF,PMAP_WF,400.E3,OROG)
c	if (nofst.ne.0) then
c	   call ecrit_fst(iun51,orog,'MX','TOPO    ',0,0,
c	1	0,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
c	endif
	HIGH_TERRAIN_FACT=1.0
	IF(LEVEL.EQ.500)THEN
	   HIGH_TERRAIN_FACT=0.0
	   PRINT *,'No attenuation near high terrain ...'
	ENDIF
c**********************************************************************
c       Time parameters
c**********************************************************************
c       Fill integer arrays YEAR,MONTH,DAY,HOUR, dimensioned NTIMES 
c       with all the times of the data. You'll use these later
	call get_times(year,month,day,hour,nstamp,ntimes,nincr,hi,hf,
	1    iun01,'GZ',level,ntw)
	print*
	print*,'YEAR 1=',YEAR(1)
	print*
	N_INCR=24/naj
	IF(N_INCR.LE.0.OR.N_INCR.GT.24) then
	   print*
	   print*,'MAKE_TRACKS: Invalid time increment'
	   STOP
	else
	   print*
	   print*,'Les donnees sont a toutes les :',N_INCR,'heures'
	   print*
	endif
c	print*, NTIMES
	IF(NTIMES.EQ.0.OR.NTIMES.GT.NTW) then
	   print*
	   print*,'MAKE_TRACKS: Invalid end time found'
	   STOP
	endif
c       last time is YEAR(NTIMES), MONTH(NTIMES) etc
c	Definition du premier IP2
	if (anal.ne.0) then
c	   hi = 0
c	   hf = 0
c	   call get_hihf(hi, hf, tstep, iun01, getvar)
	   ip2 = hi - N_INCR
	else
	   ip2 = 0
	endif
c	ip2=0
c**********************************************************************
c       OPEN OUTPUT FILE - this contains the track data you will create
c**********************************************************************
c
	print*
	print*,'Open unformatted output file and loop through times ...'
	print*
        OPEN(UNIT=2,FILE=TXT_FILE,STATUS='unknown',FORM='FORMATTED',
	1    RECL=20000)
	IYYMMDDHH_START=(YEAR(1)-1900)*1000000+MONTH(1)*10000+
	1    DAY(1)*100+HOUR(1)
c       header record contains hemisphere (IHEM), cyclones or a/c (ISIGN),
c       time increment (N_INCR), start time (IYYMMDD etc), and name of data 
c       file (WF_NAME). You'll need these later after you've made the tracks,
c       for track processing
c       Unformatted write
	print*
	print*,IHEM,ISIGN,N_INCR,IYYMMDDHH_START,WF_NAME
	print*
c----------------------------------------------------------------------
c       initialse counts
c----------------------------------------------------------------------
	NTRACKS=0
	N_TRACK=0
	N_CENTERS=0
	nggg=0
	DO IC=1,NC	
	   ICM(IC)=0
	   TUSED(IC)=.FALSE.
	END DO
	IMIN=0
	ISEC=0
	IYEAR_PREV=9999
	tf=0
	do n=1,ntimes
	   tf=tf+1
	enddo
c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c----------------------------------------------------------------------
c                     Cyclone finding algorithm
c----------------------------------------------------------------------
c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	print*,'NTIMES', NTIMES
	DO N=1,NTIMES
c**********************************************************************
c----------------------------------------------------------------------
c       1) Lecture et traitements des champs
c----------------------------------------------------------------------
c**********************************************************************
c----------------------------------------------------------------------
c       Definition de IP2 et NPAS pour ECRIT_FST
c----------------------------------------------------------------------	
	  if (anal.ne.0) then
	     ip2 = ip2 + N_INCR
	     if (deet.ne.0) then
		npas = nint(float(ip2)*(3600.0/float(deet)))
	     else
		npas = 0
	     endif
	     if (anal.eq.2) then
		
		dateo = nstamp(1)
	     endif
	  else
	     
	     npas = 0
	     dateo = nstamp(n)
	  endif
	  print*
	  print*,'==============>>>  Pas de temps :', ip2
	  print*
c----------------------------------------------------------------------
c       Definition de lheure et dela date: HH DD-MMM-YYYY
c----------------------------------------------------------------------
c       tell data file what time you want (substitute your call)
c       CALL FSTIME(YEAR(N),MONTH(N),DAY(N),HOUR(N),IMIN,ISEC)
	  IF(YEAR(N).GT.IYEAR_PREV)THEN
	     PRINT *,'Starting year ',YEAR(N) 
	  ENDIF
c       IYEAR_PREV=YEAR(N)
c       just to make sure correct times are read
	  WRITE(*,110) HOUR(N),DAY(N),MNTH(MONTH(N)),YEAR(N)
	  TIME(N)=(YEAR(N)-1900)*1000000+MONTH(N)*10000+DAY(N)*100+HOUR(
     %N)
 110      FORMAT(' Time: ',I2.2,'Z ',I2.2,'-',A3,'-',I4)
 111      FORMAT(4I2.2)
c----------------------------------------------------------------------
c       Lecture du Geopotentiel et du Vent
c----------------------------------------------------------------------
c       Read Geopotential and upper level wind from data workfile
	facteur=1.00/2.00
c	facteur2=0.5147105/2.0
c        facteur2=1.
        facteur2=1/0.514444
c******************************************************************
c        Lecture GZ, UU, VV a 1000 hPa et a 500 hPA..... 
c        Modifier coeff si necessaire      
c******************************************************************
          print *, 'test1 -------------------------------'
          print *, 'iun = ',iun01
          print *, 'lev = ',level
c          call read_field(H,nstamp(n),iun01,'GZ',9.80665,level,ni,nj,ijd)
c          print *, 'test2 -------------------------------'
c
c          call read_field(Ug,nstamp(n),iun01,'UU',1.*facteur2,level, ni,nj,ijd)
c          call read_field(Vg,nstamp(n),iun01,'VV',1.*facteur2,level, ni,nj,ijd)
c          call read_field(U,nstamp(n),iun01,'UU',(1.*facteur2/2.0),500, ni,nj,ijd)
c          call read_field(V,nstamp(n),iun01,'VV',(1.*facteur2/2.0),500, ni,nj,ijd)
          call read_field(H,iun01,'GZ',9.80665,level,ip2,ni,nj,ijd)
          print *, 'test2 -------------------------------'
          call read_field(Ug,iun01,'UU',1.*facteur2,level,ip2,ni,nj,ijd)
          call read_field(Vg,iun01,'VV',1.*facteur2,level,ip2,ni,nj,ijd)
          call read_field(U,iun01,'UU',0.50000*facteur2,500,ip2,ni,nj,ij
     %d)
          call read_field(V,iun01,'VV',0.50000*facteur2,500,ip2,ni,nj,ij
     %d)
c**********************************************************************
c       Calcul du vent GRADIENT
c**********************************************************************
c          CALL FGRAD(H,V1,U1)
c          CALL GRADIENT_WIND_FACTOR(H,WORK1)
c             DO IJ=1,NIJ
c                V1(IJ)= V1(IJ)*WORK1(IJ)/FC(IJ)
c                U1(IJ)=-U1(IJ)*WORK1(IJ)/FC(IJ)
cc               WIND_SPD(IJ)=SQRT(UG(IJ)**2+VG(IJ)**2)
c             END DO
c       Ici, le vecteur WIND_SPD sert seulement de vecteur de travail...
c**********************************************************************
c       Calcul du TOURBILLON et de la PNM
c**********************************************************************
c       Compute VORT(IJ), PRESS, WIND
         U1=Ug
         V1=Vg
         if (prcnt.ne.-1)goto 1100
          CALL FGRAD(U1,WORK3,WORK1)
          CALL FGRAD(V1,WORK2,WORK3)
          DO J=1,NJ
          DO I=1,NI
             IJ=I+(J-1)*NI	
c            VORT(IJ)=gg(ij)	!  tourbillon geostrophique
             VORT(IJ)=WORK2(IJ)-WORK1(IJ)
c            print*,VORT(IJ),WORK2(IJ)
             IF(SLAT(IJ).LT.0.0)VORT(IJ)=-VORT(IJ)
          END DO
          END DO
c       Ici, les vecteurs PRESS et WIND_SPD servent seulement de vecteurs 
c       de travail...
1100	  if (pnmi.eq.0) then
	     call read_field(press,iun01,'PN',1.0,0,ip2,ni,nj,ijd)
	  else
c       This code is to get a temperature based on month/latitude
c       for the conversion of H 1000 to MSL Press (PRESS)
	     COS_MONTH=COS(((MONTH(N)-1)/12.)*2.0*3.14159) 
	     DO J=1,NJ
		   DO I=1,NI
		    IJ=I+(J-1)*NI
		    SLAT15=ABS(SLAT(IJ))/15.
		    C1=0.5*(5.+IHEM*COS_MONTH) 
		    TEMP=299.0-C1*(SLAT15*(SLAT15-1.))
		    if (pnmi.ne.0) PRESS(IJ)=level*(1.+H(IJ)/(R*TEMP)) 
		   END DO
		   END DO
	  endif
          if (SRAD.GT.0) then
	     CALL SMOOTH_C (PRESS,NI_WF,NJ_WF,MAP_WF,PMAP_WF,SRAD,PRESS)
	  endif
c**********************************************************************
c       Filtrage
c**********************************************************************
	  pressnf = press 
	  hnf = h 
          if (SRAD.GT.0) then
	     CALL SMOOTH_C (VORT,NI_WF,NJ_WF,MAP_WF,PMAP_WF,SRAD,VORT)
	     CALL SMOOTH_C (PRESS,NI_WF,NJ_WF,MAP_WF,PMAP_WF,SRAD,PRESS)
	     CALL SMOOTH_C (h,NI_WF,NJ_WF,MAP_WF,PMAP_WF,SRAD,h)
	  endif
c**********************************************************************
c       Variables pour spline cubique et initialisation de compteurs
c**********************************************************************
c       continuous bicubic spline fit to pressure and vorticity field
      	  CALL BICUBIC(VORT,X,Y,NI,NJ,C_VORT)	
      	  CALL BICUBIC(PRESS,X,Y,NI,NJ,C_PRESS)
	  IC=0	
	  IP=0	
c**********************************************************************
c----------------------------------------------------------------------
c       2) Recherche des CENTRES
c----------------------------------------------------------------------
c**********************************************************************
c       Uses cubic splines to more accurately locate centers between 
c       grid points of computational domain
	  DO J=2,NJ-1
	  DO I=2,NI-1
	    IJ=I+(J-1)*NI
	    ALAT=ABS(SLAT(IJ))	
	    ALON=SLON(IJ)
C       Progressively inhibit finding centers over high terrain 
c       (only for H 1000)
C       This is because MSL P is unreliable over high terrain
C       HIGH_TERRAIN increases from 0 at and below 1000. m
 	    HIGH_TERRAIN(ij)=HIGH_TERRAIN_FACT*MAX(0.0,(OROG(IJ)-1000.)/300.)*
	1	 500.E3/MAX(SRAD,500.E3)
c	    IF (OROG(IJ).GT.1000.) then 
c               HIGH_TERRAIN=9.9E10 !or a sudden kill at 1000. m
c	    endif
c       Ajout pour garder tout le pouvoir de nettoyage en terrain 
c       montagneux lorsque l'on abaisse vc (pour plus de details en 
c       terrain plat)
	    if ( (HIGH_TERRAIN(ij).ne.0.0) .and. (vc.lt.2.0) ) then
	       HIGH_TERRAIN(ij) = 2.0 - vc + HIGH_TERRAIN(ij)
	    endif
c       Centers with (ABS(vort) < VCRIT are not counted
c       VCRIT gets bigger over high terrain, meaning only strong centers
c       get counted
	    VCRIT=(vc+HIGH_TERRAIN(ij))/1.E5
c       Creation d'un champ (vct) pour visualiser les cyclones suivient,
c       soit ceux ayant un max de vort > vcrit
	    if (vort(ij).gt.vcrit) then
	       vct(ij) = vort(ij)
	    else
	       vct(ij) = 0.0
	    endif
c       Note that VORT is positive for cyclones, negative for anticyclones,
c       both hemispheres
c       ---> Make VP always positive (ISIGN=1 cyclones, -1 anticyclones)
	    VP=VORT(IJ)*ISIGN
	    PP=PRESS(IJ)*ISIGN
c       Definition des 8 points de grille environnants
	    IP1J=IJ+1
	    IM1J=IJ-1
	    IJP1=I+J*NI
	    IJM1=I+(J-2)*NI
	    IP1JP1=IJP1+1
	    IP1JM1=IJM1+1
	    IM1JP1=IJP1-1	
	    IM1JM1=IJM1-1
	    FOUND=.FALSE.
c----------------------------------------------------------------------
c       2A) Recherche des CENTRES de PNM
c----------------------------------------------------------------------
c       Find extrema on IJ grid then use bicubic spline to refine
c       location of pressure extrema. Only search poleward of lat 20.
	    IF (ALAT.GT.20.0.AND.
	1	 PP.LT.PRESS(IP1J)*ISIGN.AND.
	2	 PP.LT.PRESS(IM1J)*ISIGN.AND.
	3	 PP.LT.PRESS(IJP1)*ISIGN.AND.
	4	 PP.LT.PRESS(IJM1)*ISIGN.AND.
	5	 PP.LT.PRESS(IP1JP1)*ISIGN.AND.
	6	 PP.LT.PRESS(IP1JM1)*ISIGN.AND.
	7	 PP.LT.PRESS(IM1JP1)*ISIGN.AND.
	8	 PP.LT.PRESS(IM1JM1)*ISIGN) THEN
	      IF(.NOT.VORT_C)FOUND=.TRUE.
	      IP=IP+1
	      XPP=FLOAT(I)
	      YPP=FLOAT(J)
	      PMIN=1.E30
	      DENOM=1.
	      DO IPASS=1,2
	         DENOM=2.*DENOM
	         DO II=2,4
	         DO JJ=2,4
	            XP=XPP+FLOAT(II-3)/DENOM
	            YP=YPP+FLOAT(JJ-3)/DENOM
                    CALL BICUBEV(PRESS,X,Y,NI,NJ,C_PRESS,XP,YP,PDS_PRESS
     %)
	            PTEST=PDS_PRESS(1)*ISIGN
	            IF(PTEST.LT.PMIN) THEN
	               PMIN=PTEST
	               XPF=XP
	               YPF=YP
	            ENDIF
	         END DO
	         END DO
	         XPP=XPF
	         YPP=YPF
	      END DO
	      CALL LL(MAP,PMAP,XPF,YPF,PLAT(IP),PLON(IP))
	   ENDIF
c----------------------------------------------------------------------
c       2B) Recherche des CENTRES de TOURBILLON
c----------------------------------------------------------------------
c       Find vorticity extrema
	    IF(VORT_C)THEN
	     IF(ALAT.GT.20.0.AND.VP.GT.vcrit.and.
	1	    vp.gt.VORT(IP1J)*ISIGN.AND.
	2	    VP.GT.VORT(IM1J)*ISIGN.AND.
	3	    VP.GT.VORT(IJP1)*ISIGN.AND.
	4	    VP.GT.VORT(IJM1)*ISIGN.AND.
	5	    VP.GT.VORT(IP1JP1)*ISIGN.AND.
	6	    VP.GT.VORT(IP1JM1)*ISIGN.AND.
	7	    VP.GT.VORT(IM1JP1)*ISIGN.AND.
	8	    VP.GT.VORT(IM1JM1)*ISIGN) THEN
	          FOUND=.TRUE.
		  XPP=FLOAT(I)
		  YPP=FLOAT(J)
		  VMAX=-1.E30
		  DENOM=1.
		  DO IPASS=1,2
		     DENOM=2.*DENOM
		     DO II=2,4
		     DO JJ=2,4
		        XP=XPP+FLOAT(II-3)/DENOM
		        YP=YPP+FLOAT(JJ-3)/DENOM
		        CALL LL(MAP,PMAP,XP,YP,PPLAT,PPLON)
      		        CALL BICUBEV(VORT,X,Y,NI,NJ,C_VORT,XP,YP,PDS_VORT)
		        VTEST=ISIGN*PDS_VORT(1)
		        IF(VTEST.GT.VMAX) THEN
			   VMAX=VTEST
			   XPF=XP
			   YPF=YP
	                   PLAT_SAV=PPLAT
	                   PLON_SAV=PPLON
		        ENDIF
		     END DO
		     END DO
		     XPP=XPF
		     YPP=YPF
		  END DO
	        ENDIF
	     ENDIF
C       write lat lon press, vort to arrays for all points found
	     IF(FOUND)THEN
		IC=IC+1
	        CALL LATLON(MAP,PMAP,XPF,YPF,TLAT(IC),TLON(IC),DDX,DDY,
	1	     DN,DE,DDA)
      		CALL BICUBEV(VORT,X,Y,NI,NJ,C_VORT,XPF,YPF,PDS_VORT)
		TVORT(IC)=IHEM*PDS_VORT(1)
      		CALL BICUBEV(PRESS,X,Y,NI,NJ,C_PRESS,XPF,YPF,PDS_PRESS)
		TPRESS(IC)=PDS_PRESS(1)
		IF(TPRESS(IC).LT.900.)PRINT *,'ERROR PRESS.LT.900'
		CALL BL1(TU(IC),XPF,YPF,U,NI,NJ)
		CALL BL1(TUg(IC),XPF,YPF,Ug,NI,NJ)
		CALL BL1(TU1(IC),XPF,YPF,U1,NI,NJ)
		CALL BL1(TV(IC),XPF,YPF,V,NI,NJ)
		CALL BL1(TVg(IC),XPF,YPF,Vg,NI,NJ)
		CALL BL1(TV1(IC),XPF,YPF,V1,NI,NJ)
	     ENDIF
	  END DO
	END DO
c       HERE, we have all the centers.
	N_NEW=IC
	PRINT *,'  # new CENTER positions:',N_NEW
c----------------------------------------------------------------------
c       Ecriture de differents champs en Fichier Standard
c----------------------------------------------------------------------
	goto 1236
	if (deet.eq.0)then
	npas=0
	else
	npas = nint(float(ip2)*(3600.0/float(deet)))
	endif
	if (nofst.ne.0) then
cc       1) Tourbillon du vent gradient
	call ecrit_fst(iun51,vct,'QR','GR-VCRIT',level,ip2,
	1    npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,vort,'QR','GRADIENT',level,ip2,
	1    npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	endif
c       2) Pression au niv de la mer
	if (pnmi.eq.0) then
	   call ecrit_fst(iun51,press,'PN','PNM_MD  ',0,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	else
	   call ecrit_fst(iun51,press,'PN','Z10002PN',0,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	endif
	call ecrit_fst(iun51,pressnf,'PN','PNM_NF  ',0,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
c       5) Vent geostropique à 1000 hPa et 500hpa
	call ecrit_fst(iun51,ug,'UU','GEOPT   ',level,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,vg,'VV','GEOPT   ',level,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,u,'UU','MRC   ',500,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,v,'VV','MRC   ',500,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,hnf,'GZ','MRC   ',level,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
	call ecrit_fst(iun51,gg,'GE','GEOP   ',level,ip2,
	1	npas,dateo,ni,nj,ig1,ig2,ig3,ig4,ip3,grtyp,typvar,deet)
1236    continue
c----------------------------------------------------------------------
c       3.1) Calcul de la direction du vent
c----------------------------------------------------------------------
	 DO IC=1,N_NEW
           TAN(IC)=0.0
           DIRC1(IC)=0.0
           PHI(IC)=0.0
           BETA(IC)=0.0
	   IG=IG3/100.
	     TAN(IC)=TVg(IC)/TUg(IC)	
	     DIRC1(IC)=atand(TAN(IC))
             if(TVg(IC).eq.0.and.TUg(IC).lt.0)DIRC1(IC)=0.0000
	     if(TUg(IC).eq.0.and.TVg(IC).gt.0)DIRC1(IC)=90.000
	     if(TVg(IC).lt.0.and.TUg(IC).lt.0)DIRC1(IC)=DIRC1(IC)+180.
	     if(TVg(IC).lt.0.and.TUg(IC).gt.0)DIRC1(IC)=DIRC1(IC)+360.
	     if(TUg(IC).lt.0.and.TVg(IC).gt.0)DIRC1(IC)=DIRC1(IC)+180.
        PHI(IC)=360.-TLON(IC)-IG
        BETA(IC)=90.-PHI(IC)
        DIRC(IC)=270.-DIRC1(IC)+BETA(IC)
        IF(DIRC(IC).LT.0.)DIRC(IC)=DIRC(IC)+360.0
        WIND(IC)=SQRT(TUg(IC)*TUg(IC)+TVg(IC)*TVg(IC))*3.6	
c	print*,IP2,DIRC(IC),WIND(IC),TLON(IC),TLAT(IC)
c	pause
        ENDDO
c**********************************************************************
c----------------------------------------------------------------------
c       3.2) Calcul de la circulation
c----------------------------------------------------------------------
c**********************************************************************
c       For each center found, work out its circulation
c       First, we need to define the boundary of the area 
c       to do the circulation calcls -l 
	  NTHETA=0
	  DO IC=1,N_NEW
	     IF(DO_CIRC)THEN
C       see if close to another center. We need that to check later for
C       overlapping areas
	        TDIST_MIN=3500.E3
	        NCL=0
		IF(IC.GE.2)THEN
		    DO ICC=1,IC-1
		       CALL GCPAZD(TLAT(ICC),TLON(ICC),TLAT(IC),TLON(IC),
	1		    TDAZ,TDIST)
		       IF(TDIST.LT.TDIST_MIN) THEN
			  NCL=NCL+1
			  ICL(NCL)=ICC
		       ENDIF
		    END DO
		ENDIF
	        NCLOSE=NCL  
C       set some constants
		DRAD=1.
		IDTHETA=20.
C       initialise AREA and FCIRC for inner points
	        AREA=PI*DRAD**2/4.
	        TAREA_20=AREA
		TAREA=AREA
c               FCIRC=1.E-5*TVORT(IC)*AREA
                TCIRC(IC)=TVORT(IC)*AREA
c       move azimuthally around polar grid 360 deg
c       actually go around a bit further - this helps the smoothing of 
c       the boundary
c       ignore for circulation calculation after IT > NTHETA
		NTHETA=0
                DO ITHETA=IDTHETA,360,IDTHETA
		   NTHETA=NTHETA+1
		END DO
C       counters
		IT=0	
		ITC=0
                ICT=0
                DO ITHETA=-IDTHETA,360,IDTHETA
		   IT=IT+1
		   IF(ITHETA.GE.IDTHETA)ITC=ITC+1
		   THETA=FLOAT(ITHETA)
		   RAD_TEST=99.
C       RAD_LIM is the radius to start looking outward from 
		   RAD_LIM=5.0
C		   VTHR=0.05*1000.E3/SRAD
		   VTHR=0.0
		   DAZ_DIFF_MIN=999.
C       step radially outward 
		   IR=0
		   RADIUS=0.0
	           DO WHILE (RAD_TEST.EQ.99.)
		       RADIUS=RADIUS+DRAD
		       IR=IR+1
c		       DAREA=2.*PI*RADIUS*DRAD*FLOAT(IDTHETA)/360.
                       DAREA_1=2.*PI*(RADIUS-DRAD)**2*FLOAT(IDTHETA)/360
     %.
                       DAREA_2=2.*PI*RADIUS**2*FLOAT(IDTHETA)/360
                       DAREA=DAREA_2-DAREA_1
		       CALL GCPMOV(TLAT(IC),TLON(IC),SSLAT,SSLON,THETA,
	1		    RADIUS*1.111E5)
		       CALL PXY(MAP,PMAP,SSLAT,SSLON,XP,YP)
C       see if <20° lat
      		       IF(ABS(SSLAT).GT.21..and.yp.ge.y(1).and.xp.ge.
	1		    x(1).and.yp.le.y(nj).and.xp.le.x(ni)) THEN
			  CALL BICUBEV(VORT,X,Y,NI,NJ,C_VORT,XP,YP,PDS_VORT)
			  CALL BICUBEV(PRESS,X,Y,NI,NJ,C_PRESS,XP,YP,PDS_PRESS)
			  CALL BL1(RU1(IR),XP,YP,U1,NI,NJ)
			  CALL BL1(RV1(IR),XP,YP,V1,NI,NJ)
		          TAREA_20=TAREA_20+DAREA
		       ELSE
		          RAD_TEST=RADIUS		
		       ENDIF
		       FVORT=PDS_VORT(1)
		       RVORT(IR)=FVORT
		       FPRESS=PDS_PRESS(1)
		       RPRES(IR)=FPRESS
		       TAREA=TAREA+DAREA
CCCCCCCCCCCCCCCCCC CALCUL DE LA DIRECTION DU VENT 
	  CALL LATLON(MAP,PMAP,XP,YP,XLAT(IR),XLON(IR),DDX,DDY,
	1	     DN,DE,DDA)
	   TAN(IR)=0.0
           DIRC1(IR)=0.0
           PHI(IR)=0.0
           BETA(IR)=0.0
           RDIR1(IR)=0.
           RDIR(IR)=0.
	     TAN(IR)=RV1(IR)/RU1(IR)	
	     RDIR1(IR)=atand(TAN(IR))
             if(RV1(IR).eq.0.and.RU1(IR).lt.0)RDIR1(IR)=0.0000
	     if(RU1(IR).eq.0.and.RV1(IR).gt.0)RDIR1(IR)=90.000
	     if(RV1(IR).lt.0.and.RU1(IR).lt.0)RDIR1(IR)=RDIR1(IR)+180.
	     if(RV1(IR).lt.0.and.RU1(IR).gt.0)RDIR1(IR)=RDIR1(IR)+360.
	     if(RU1(IR).lt.0.and.RV1(IR).gt.0)RDIR1(IR)=RDIR1(IR)+180.
        PHI(IR)=360.-XLON(IR)-25.0
        BETA(IR)=90.-PHI(IR)
        RDIR(IR)=270.-RDIR1(IR)+BETA(IR)
        IF(RDIR(IR).LT.0.)RDIR(IR)=RDIR(IR)+360.0
        RWIND(IR)=SQRT(RU1(IR)*RU1(IR)+RV1(IR)*RV1(IR))*3.6
CCCCCCCCCCCCCCCC
c       see if FVORT < 0
		       IF (RADIUS.GE.RAD_LIM.AND.FVORT.LT.0.0) then
			  RAD_TEST=RADIUS-DRAD
		       endif
c       see vort increasing more than VTHR - 
c       dont start looking until RADIUS > RAD_LIM
		       VCH=FVORT_PREV-FVORT  
		       IF(RADIUS.GE.RAD_LIM.AND.VCH.LT.VTHR) then
			  RAD_TEST=RADIUS-DRAD
		       endif
C       see if inside another cyclone's domain	
		       IF (RADIUS.GE.RAD_LIM.AND.
	1		    IC.GT.1.AND.
	2		    NCLOSE.GT.0.AND.
	3		    RAD_TEST.EQ.99.) THEN 
			   DO NCL=1,NCLOSE
			      ICC_T=ICL(NCL)
			      INS=.FALSE.
			      DO ITT=1,NTHETA
				 CALL PXY(MAP,PMAP,BLAT(ITT,ICC_T),
	1			      BLON(ITT,ICC_T),XDRAW(ITT),YDRAW(ITT))
			      END DO
			      CALL INSIDE(XDRAW,YDRAW,NTHETA,XP,YP,INS,XC,YC)
			      IF(INS)RAD_TEST=RADIUS-DRAD
			   END DO
	               ENDIF
		       FVORT_PREV=FVORT
		   END DO   
C       boundary is RAD_TEST
		   RAD_BOUND=RAD_TEST
C       boundary radius can only change by BDD per azimuth step
		   BDD=NINT(RAD_PREV*FLOAT(IDTHETA)/60.)*DRAD
		   IF(DAZ_DIFF_MIN.GT.45.0.AND.IT.GT.1)THEN
		      IF((RAD_BOUND-RAD_PREV).GT.BDD)RAD_BOUND=RAD_PREV+BDD
		      IF((RAD_PREV-RAD_BOUND).GT.BDD)RAD_BOUND=RAD_PREV-BDD
		   ENDIF
	           RAD_PREV=RAD_BOUND
C       compute blat blon at boundary
		   CALL GCPMOV(TLAT(IC),TLON(IC),BLATT,BLONN,THETA,
	1		RAD_BOUND*1.111E5)
C       compute circulation out to this boundary
		   IF(ITC.GT.0)THEN
		      IR=0
                       FCIRC=0
	              DO RADIUS=DRAD,RAD_BOUND,DRAD
		         IR=IR+1
                       DAREA_1=2.*PI*(RADIUS-DRAD)**2*FLOAT(IDTHETA)/360
     %.
                       DAREA_2=2.*PI*RADIUS**2*FLOAT(IDTHETA)/360.
                       DAREA=DAREA_2-DAREA_1
                         FCIRC=FCIRC+DAREA*IHEM*RVORT(IR)
ccccccccccccccccccccccccccccccc
		   CALL GCPMOV(TLAT(IC),TLON(IC),BLATT,BLONN,THETA,
	1		RADIUS*1.111E5)
		      ICT=ICT+1
		      RLATT(IC,ICT)=BLATT
		      RLONN(IC,ICT)=BLONN
                      RWINDD(IC,ICT)=RWIND(IR)
		      RDIRR(IC,ICT)=RDIR(IR)
                      RPRESS(IC,ICT)=RPRES(IR)
cccccccccccccccccccccccccccccccc
		      END DO
		      BLAT(IC,ITC)=BLATT
		      BLON(IC,ITC)=BLONN
		   ENDIF
                TCIRC(IC)=TCIRC(IC)+FCIRC
		END DO	
		NUM1(IC)=ICT
		NUM2(IC)=ITC
                TCIRC(IC)=1.2345E10* TCIRC(IC)
c	        TCIRC(IC)=1.E-7*1.2345E10*FCIRC*TAREA/TAREA_20    !normalize
	     ELSE	
c	        CALL POLAR_AV(VORT,NI,NJ,MAP,PMAP,TLAT(IC),TLON(IC),
c	1	     5.,FCIRC,IERR)
C	        TCIRC(IC)=FCIRC*IHEM
	     ENDIF
	     CALL CHLALO(TLAT(IC),TLON(IC),FCHAR)
	  END DO
c**********************************************************************
c----------------------------------------------------------------------
c       4) Systeme OUVERT ou FERME 
c----------------------------------------------------------------------
c**********************************************************************
c       See which vorticity centers are close to pressure extrema. 
c       Call these closed
	  DO IC=1,N_NEW
		TCLOSED(IC)=.FALSE.
	  END DO
	  NPP=IP
	  DO IP=1,NPP
	     TDIST_MIN=1.E30
	     DO IC=1,N_NEW
		CALL GCPAZD(PLAT(IP),PLON(IP),TLAT(IC),TLON(IC),TDAZ,TDIST)
		IF(TDIST.LT.TDIST_MIN) THEN
		    TDIST_MIN=TDIST
		    IC_C=IC
		ENDIF
	     END DO
	     IF(TDIST_MIN.LT.3.3E5)TCLOSED(IC_C)=.TRUE.
	  END DO
	  IF(.NOT. VORT_C)THEN	
	     DO IC=1,N_NEW
		TCLOSED(IC)=.TRUE.
	     END DO
	  ENDIF
c**********************************************************************
c----------------------------------------------------------------------
c       5) Contruction des trajectoires 
c----------------------------------------------------------------------
c**********************************************************************
	  IF(N.EQ.1) GOTO 1000
c       Go through all CENTERs for which a track is being constructed ...
	  IPASS=0
2000      CONTINUE
	  IPASS=IPASS+1
	  DIST_CRIT_FACT=9.-1.0*(IPASS-1)  
	  DO IC=1,NC	
	     CUSED(IC)=.FALSE.
	     TUSED(IC)=.FALSE.
	  END DO
	  SPEED_FACT=FLOAT(N_INCR)/12. 
	  WP=0.55**SPEED_FACT 
	  WV=0.65**SPEED_FACT 
	  N24=NINT(24./FLOAT(N_INCR))	
c
	  IMATCH=0
	  DO IT=1,NTRACKS   
	      PMAX=0.0
	      NN=NNT(IT)
	      AALAT=ABS(CLAT(NN,IT))
	      dd=0.0
	      degtor=PI/180.0
	      IF(CU(NN,IT).NE.0.0.OR.CV(NN,IT).NE.0.0) then
		 DD=ATAN2(CU(NN,IT),CV(NN,IT))/degtor+180.0
	      endif
	      IF(Dd.LT.0.0) dd=dd+360.0
	      IF(DIR.GT.360.0) dd=dd-360.0
	      ff=SQRT(CU(NN,IT)*CU(NN,IT)+CV(NN,IT)*CV(NN,IT))
	      DAZ_AV=DD+180.
	      IF(DAZ_AV.GT.360.)DAZ_AV=DAZ_AV-360.
	      DIST_AV=SPEED_FACT*FF*4.32E4 
	      X_av =-1.0*DIST_av*SIN(DAZ_av*DEGTOR)
	      Y_av=-1.0*DIST_av*COS(DAZ_av*DEGTOR)
C       Previous motion over 24 h
	      IF(NN.GT.1) THEN
	         NPREV=MAX(1,NN-N24)
	         WM=0.70**SPEED_FACT 
	         NDIF=NN-NPREV
	         CALL GCPAZD(CLAT(NPREV,IT),CLON(NPREV,IT),CLAT(NN,IT),
     -                        CLON(NN,IT),DAZ_PREV,DIST_PREV)
		 X_PREV =-1.0*DIST_PREV*SIN(DAZ_PREV*DEGTOR)
		 Y_PREV=-1.0*DIST_PREV*COS(DAZ_PREV*DEGTOR)
	         X_EST=WM*X_PREV/FLOAT(NDIF)+(1.-WM)*X_AV
	         Y_EST=WM*Y_PREV/FLOAT(NDIF)+(1.-WM)*Y_AV
	         P_EST=CPRESS(NN,IT)+WP*(CPRESS(NN,IT)-CPRESS(NN-1,IT))
	         V_EST=CVORT(NN,IT)+WV*(CVORT(NN,IT)-CVORT(NN-1,IT))
	      ELSE 
	         X_PREV=0.0
		 Y_PREV=0.0
	         X_EST=X_AV
	         Y_EST=Y_AV
	         P_EST=CPRESS(NN,IT)	
	         V_EST=CVORT(NN,IT)     
	      ENDIF
	      DAZ_EST=0.0
	      IF(X_EST.NE.0.0.OR.Y_EST.NE.0.0) then
		 DAZ_EST=ATAN2(X_EST,Y_EST)/degtor+180.0
	      endif
	      IF(DAZ_EST.LT.0.0) then
		 DAZ_EST=DAZ_EST+360.0
	      endif
	      IF(DAZ_EST.GT.360.0) then
		 DAZ_EST=DAZ_EST-360.0
	      endif
	      DIST_EST=SQRT(X_EST*X_EST+Y_EST*Y_EST)
	      if (DIST_EST.eq.0.E0) then
		 PLAT_EST=CLAT(NN,IT)
		 PLON_EST=CLON(NN,IT)-360.
	      else
		 CALL GCPMOV(CLAT(NN,IT),CLON(NN,IT),PLAT_EST,PLON_EST,
	1	      DAZ_EST,DIST_EST)
	      endif
	      RPSQ_CRIT=DIST_CRIT**2
	      DO IC=1,N_NEW   
	         PMN(IC,IT)=-999.
	         CALL GCPAZD(PLAT_EST,PLON_EST,TLAT(IC),TLON(IC),TDAZ,TD
     %IST)
c       Examine all CENTERs within DIST_CRIT of predicted next track position
c       and assign probability of association
	         TDIST=TDIST/1.11E5
	         IF(TDIST.LT.DIST_CRIT)THEN
		     SINP=ABS(SIND(TLAT(IC)))
		     PSCALE=0.55/SINP
		     PDIST=PSCALE*(TPRESS(IC)-P_EST)
	             IF(VORT_C)THEN
			 SCALE_INTENS=1.0+ABS(TVORT(IC))/5.
		         VSCALE=1.6
			 VDIST=VSCALE*(TVORT(IC)-V_EST)
			 PV_SQ=0.5*(PDIST**2+VDIST**2)/SCALE_INTENS
		     ELSE
			 PV_SQ=PDIST**2
		     ENDIF
	             RPSQ=TDIST**2+PV_SQ
	             PROB=1.0-RPSQ/RPSQ_CRIT
C                    probable match found
	             IF(PROB.GT.0.0)THEN
	                 PMN(IC,IT)=PROB
		         IMATCH=IMATCH+1
	                 ICM(IMATCH)=IC
	                 ITM(IMATCH)=IT
	             ENDIF
	         ENDIF
	      END DO		
	   END DO		
	  N_MATCHES=IMATCH
c         print*, 'N_MATCHES=',N_MATCHES, ic, it
c         write (*, 123) nn, ic, it
 123	  format (2x, 3(i4))
c------------------------------------------------------------------------
c       Organize matches into groups ICG,ITG(I,NN): I=group #, NN=element
c------------------------------------------------------------------------
	  IG=0
	  DO IM=1,N_MATCHES	
	    IC=ICM(IM)
	    IT=ITM(IM)
	    ITM(IM)=0
	    ICM(IM)=0
	    FOUND=.FALSE.
	    DO I=1,IG	
	        NGG=NG(I)
	        DO NN=1,NGG
	          IF(IC.EQ.ICG(I,NN).AND.IT.EQ.ITG(I,NN))STOP 'two group
     %s same'
	          IF(.NOT.FOUND.AND.(IC.EQ.ICG(I,NN).OR.IT.EQ.ITG(I,NN))
     %)THEN
	              FOUND=.TRUE.
	              NG(I)=NG(I)+1
		      ICG(I,NG(I))=IC
		      ITG(I,NG(I))=IT
	          ENDIF
	        END DO	
	    END DO
	    IF(.NOT.FOUND) THEN 
		IG=IG+1
		NG(IG)=1
		ICG(IG,1)=IC
		ITG(IG,1)=IT
	    ENDIF
	  END DO    
	  N_GROUPS=IG
c       Merge groups with elements in common
 200	  NGP=N_GROUPS
	  DO IG=1,NGP-1
	    NGG=NG(IG)
	    DO NN=1,NGG
	      IC=ICG(IG,NN)
	      IT=ITG(IG,NN)
	      DO IG_T=IG+1,NGP
	        NGG_T=NG(IG_T)
	        DO NN_T=1,NGG_T
		  IC_T=ICG(IG_T,NN_T)
		  IT_T=ITG(IG_T,NN_T)
	          IF(IC_T.EQ.IC.AND.IT_T.EQ.IT)PRINT *,'ERROR IN DOUBLE 
     %MATCH'
	          IF(IC_T.EQ.IC.OR.IT_T.EQ.IT)THEN
	             IG_FOUND=IG
	             IG_MATCH=IG_T
	             DO N2=1,NG(IG_FOUND) 
			ICG(IG_MATCH,N2+NG(IG_MATCH))=ICG(IG_FOUND,N2)
			ITG(IG_MATCH,N2+NG(IG_MATCH))=ITG(IG_FOUND,N2)
		     END DO
	             NG(IG_MATCH)=NG(IG_MATCH)+NG(IG_FOUND)
	             DO N2=1,NG(N_GROUPS) 
			ICG(IG_FOUND,N2)=ICG(N_GROUPS,N2)
			ITG(IG_FOUND,N2)=ITG(N_GROUPS,N2)
		     END DO
	             NG(IG_FOUND)=NG(N_GROUPS)
	             N_GROUPS=N_GROUPS-1
	             GOTO 200
	          ENDIF
	        END DO
	      END DO
	    END DO
	 END DO
	  NG_S=0
	  DO IG=1,N_GROUPS
	     NG_S=NG_S+NG(IG)
	  END DO
	  IF(NG_S.NE.N_MATCHES) then
	     PRINT *,'ERROR: # in groups .NE. # MATCHES',NG_S
	  endif
c       now, look at all the combinations of matches and choose the one 
c       with the largest probability sum
	  DO IG=1,N_GROUPS
	     IF(NG(IG).GT.9)THEN
		print*,'***************', NG(IG)
		NG(IG)=9
		nggg=nggg+1
	     ENDIF
	  END DO
C       loop through all groups containing more than one match
	  DO IG=1,N_GROUPS
	    IF(NG(IG).GT.1) THEN
	       CALL perm(NG(IG),IORDER,N_PERM_GROUPS)
	       PMAX=-1.E30	
	       IG_PERM_SAVE=0
C       loop through these in the order of permutation IG_PERM 
	       DO IG_PERM=1,N_PERM_GROUPS
	          PSUM=0.0
	          DO NN=1,NG(IG)
	            IC=ICG(IG,NN)
	            IT=ITG(IG,NN)
		    CUSED(IC)=.FALSE.
		    TUSED(IT)=.FALSE.
	          END DO
	          DO N1=1,NG(IG)
	            NN=IORDER(N1,IG_PERM)
	            IC=ICG(IG,NN)
	            IT=ITG(IG,NN)
	            IF(.NOT.(CUSED(IC)).AND..NOT.(TUSED(IT)))THEN
	               CUSED(IC)=.TRUE.
	               TUSED(IT)=.TRUE.
		       PSUM=PSUM+PMN(IC,IT)	
	               IF(PMN(IC,IT).EQ.-999.)PRINT *,'ERROR IN PERM PRO
     %BS 1'
		    ENDIF
	          END DO    
	          IF(PSUM.GT.PMAX) THEN
	             PMAX=PSUM
	             IG_PERM_SAVE=IG_PERM
	             PSUM_SAVE=PSUM
	          ENDIF
	       END DO	    
C       reset
	       DO NN=1,NG(IG)
	          IC=ICG(IG,NN)
	          IT=ITG(IG,NN)
		  CUSED(IC)=.FALSE.
		  TUSED(IT)=.FALSE.
	       END DO
C       grab the ones with the largest total probability
	       PSUM=0.0
	       DO N1=1,NG(IG)
		  NN=IORDER(N1,IG_PERM_SAVE)
	          IC=ICG(IG,NN)
	          IT=ITG(IG,NN)
	          IF(.NOT.(CUSED(IC).OR.TUSED(IT))) THEN
		     ICM(IT)=IC
	             CUSED(IC)=.TRUE.
	             TUSED(IT)=.TRUE.
		     PSUM=PSUM+PMN(IC,IT)	
	             IF(PMN(IC,IT).EQ.-999.)PRINT *,'ERROR IN PERM PROBS
     % 2'
		  ENDIF
	       END DO
	       IF(PSUM.NE.PSUM_SAVE)PRINT *,'ERROR IN PSUM'
	    ELSE    
	       IT=ITG(IG,1)
	       IC=ICG(IG,1)
	       IF(TUSED(IT).OR.CUSED(IC))PRINT *,'ERROR - ALREADY USED'
	       TUSED(IT)=.TRUE.
	       CUSED(IC)=.TRUE.
	       ICM(IT)=IC
	       IF(PMN(IC,IT).LE.0.0)PRINT *,'ERROR IN PERM PROBS 3'
	    ENDIF   
	 END DO			
C
	  IDEF=0
	  DO IT=1,NTRACKS
	     NN=NNT(IT)
	     IC=ICM(IT)
	     CALL CHLALO(CLAT(NN,IT),CLON(NN,IT),FCHAR)
	     IF(IC.GT.0.AND.TUSED(IT)) THEN	
	        CALL GCPAZD(TLAT(IC),TLON(IC),CLAT(NN,IT),CLON(NN,IT),
	1	     CDAZ,CDIST)
	        ITDD=NINT(CDAZ/10.)	
	        ITFF=NINT(CDIST/(12*3600.*0.5144))
	        CALL CHLALO(TLAT(IC),TLON(IC),FCHAR(14:))
	        FCHAR(14:)='--->'//FCHAR(14:)	
		if (n.eq.tf) then
		   do ii=2,tf-1
		      if (NNT(IT).EQ.ii) then
			 PRINT*,it,'Trajectoire non terminee - fin forcee!'
			 TUSED(IT)=.FALSE.
		      endif
		   end do
		else
		   NNT(IT)=NNT(IT)+1
C       write information for next track point into track arrays
		   NN=NNT(IT)
		   CLAT(NN,IT)=TLAT(IC)
		   CLON(NN,IT)=TLON(IC)
		   CVORT(NN,IT)=TVORT(IC)
		   CPRESS(NN,IT)=TPRESS(IC)
		   CCIRC(NN,IT)=TCIRC(IC)
                   CWIND(NN,IT)=WIND(IC)
                   CDIRC(NN,IT)=DIRC(IC)
		   CU(NN,IT)=TU(IC)
		   CV(NN,IT)=TV(IC)
		   CLOSED(NN,IT)=TCLOSED(IC)
                   NUM3(NN,IT)=NUM1(IC)
		   NUM4(NN,IT)=NUM2(IC)
		   NTIME(NN,IT)=N
		   DO ICT=1,NUM1(IC)		
		   CLAT1(NN,IT,ICT)=0.0			
		   CLON1(NN,IT,ICT)=0.0
		   CPRESS1(NN,IT,ICT)=0.0
		   CWIND1(NN,IT,ICT)=0.0
                   CDIRC1(NN,IT,ICT)=0.0	
		   enddo
		   DO ICT=1,NUM3(NN,IT)		
		   CLAT1(NN,IT,ICT)=RLATT(IC,ICT)
		   CLON1(NN,IT,ICT)=RLONN(IC,ICT)
		   CPRESS1(NN,IT,ICT)=RPRESS(IC,ICT)
		   CWIND1(NN,IT,ICT)=RWINDD(IC,ICT)
                   CDIRC1(NN,IT,ICT)=RDIRR(IC,ICT)
		   ENDDO
		   DO ITC=1,NUM2(IC)		
		   FLAT(NN,IT,ITC)=0.0			
		   FLON(NN,IT,ITC)=0.0	
		   enddo
		   DO ITC=1,NUM4(NN,IT)		
		   FLAT(NN,IT,ITC)=BLAT(IC,ITC)
		   FLON(NN,IT,ITC)=BLON(IC,ITC)
		   ENDDO
	ENDIF
	     ELSE
	        WRITE(*,'(6X,A20,2F7.1,1x,'' track ends'')')
c     -             'Track: '//FCHAR,CPRESS(NN,IT),CVORT(NN,IT)
	     ENDIF
	  END DO
c**********************************************************************
c----------------------------------------------------------------------
c       6) Output track in TXT file
c----------------------------------------------------------------------
c**********************************************************************
	  DO IT=1,NTRACKS
	      IF(.NOT. TUSED(IT))THEN
	         NTR=NNT(IT)
				   DDHH=MOD(TIME(NTIME(1,IT)),10000)
				   YY=MOD(YEAR(NTIME(1,IT)),100)
                                   IYY=NINT(YY)
                                   IYEAR=YEAR(NTIME(1,IT))
	         DO NN=1,NTR
				 ILAT(NN)=NINT(CLAT(NN,IT)*10.)
		    		 ILON(NN)=NINT(CLON(NN,IT)*10.)
			         RLAT(NN)=CLAT(NN,IT)
                                 RLON(NN)=CLON(NN,IT)-360.
		    		 IPRESS(NN)=NINT((CPRESS(NN,IT)))
		    		 IF(CVORT(NN,IT).GT.100.)PRINT *,'VORT > 100'
		    		 IVORT(NN)=NINT(CVORT(NN,IT)*1000000.)
		    		 ICIRC(NN)=NINT(CCIRC(NN,IT)*1E-6)
				 IWIND(NN)=CWIND(NN,IT)
	                         IDIRC(NN)=CDIRC(NN,IT)
				 IDIAM(NN)=NINT((sqrt(ICIRC(NN)*1E12/(3.14*IVORT(NN))))/1000)
     %				
		NM1(NN)=NUM3(NN,IT)
		NM2(NN)=NUM4(NN,IT)
			DO ICT=1,NM1(NN)		
		                CLAT2(NN,ICT)=0.0		
		   		CLON2(NN,ICT)=0.0				
		   		CPRESS2(NN,ICT)=0.0
		   		CWIND2(NN,ICT)=0.0
                   		CDIRC2(NN,ICT)=0.0
		        END DO
		       DO ICT=1,NM1(NN)		
		                CLAT2(NN,ICT)=CLAT1(NN,IT,ICT)
		   		CLON2(NN,ICT)=CLON1(NN,IT,ICT)				
		   		CPRESS2(NN,ICT)=CPRESS1(NN,IT,ICT)
		   		CWIND2(NN,ICT)=CWIND1(NN,IT,ICT)
                   		CDIRC2(NN,ICT)=CDIRC1(NN,IT,ICT)
		        END DO
			DO ITC=1,NM2(NN)		
		                FLAT1(NN,ITC)=0.0		
		   		FLON1(NN,ITC)=0.0						   		
		        END DO
		       DO ITC=1,NM2(NN)		
		                FLAT1(NN,ITC)=FLAT(NN,IT,ITC)
		   		FLON1(NN,ITC)=FLON(NN,IT,ITC)				
		       END DO
		 END DO
				 N_CENTERS=N_CENTERS+NTR
c------------------------
c       Etape 1 
c------------------------
c       On ecrit Numero (N_TRACK), Duree (NTR) de la trajectoire
c       et Temps (DDHHYYMM) d'apparition de la trajectoire
		 if (ntr.ge.npt) then
		    N_TRACK=N_TRACK+1
		    YYMM=TIME(NTIME(1,IT))/10000 
		    write(2,'(I4,$)') ddhh
		    write(2,'(I5,$)') yymm
		    write (2,'(I6,i3,$)') n_track,ntr
                    write(16,124)n_track,MNTH(MONTH(NTIME(1,IT))),IYY,MO
     %NTH(NTIME(1,IT)),DAY(NTIME(1,IT)),YEAR(NTIME(1,IT))
		    write(17,124)n_track,MNTH(MONTH(NTIME(1,IT))),IYY,MONTH(NTI
     %ME(1,IT)),DAY(NTIME(1,IT)),YEAR(NTIME(1,IT))
c124           FORMAT(i4.4,'NAR_','             ',i2.2,'/',i2.2,'/',i4)
124           FORMAT(i4.4,' NAR_',A3,I2.2,'             ',i2.2,'/',i2.2,
     %'/',i4)
		 endif
c------------------------
c       Etape 2 
c------------------------
c       On ecrit Position, Tourbillon, Pression, etc 
c       pour chacun des pas de temps de la trajectoire
c       (sauf le dernier)
		do nn=1,ntr-1
		   nnn=NTIME(nn,IT)
		   if (ntr.ge.npt) then
	   if(DAY(nnn).lt.10)then
	   write(16,125)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(NN)),abs(
     %RLAT(NN)),IPRESS(NN),IWIND(NN)
	   write(17,125)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(NN)),abs(
     %RLAT(NN)),ICIRC(NN),IDIRC(NN)
          else
           write(16,126)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(NN)),abs(
     %RLAT(NN)),IPRESS(NN),IWIND(NN)
	   write(17,126)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(NN)),abs(
     %RLAT(NN)),ICIRC(NN),IDIRC(NN)
          endif
c Output GIS
	   write(18,127)IYY,MONTH(nnn),DAY(nnn),HOUR(nnn),n_track,RLAT(N
     %N),RLON(NN),IPRESS(NN),IVORT(NN),ICIRC(NN),IWIND(NN),IDIRC(NN),IDI
     %AM(NN)
	   write(2,'(2I4,i6,i6,i6,i6,i6$)') ILAT(NN),ILON(NN),IPRESS(NN)
     %,IVORT(NN),ICIRC(NN),IWIND(NN),IDIRC(NN)
125         FORMAT(I3,I2,I3.2,'00Z',F6.1,F5.1,i5,i5,'  *')
126         FORMAT(I3,I3,I3.2,'00Z',F6.1,F5.1,i5,i5,'  *')	
127         FORMAT(I2.2,I2.2,I2.2,I2.2,'0000',i4,F5.1,F7.1,i6,i6,i6,i6,i
     %6,i6)       	
             endif
		enddo
cccccccccc ecriture des positions des points appartenant a la tempete
		DO nn=1,ntr				
		   if (ntr.ge.npt) then	   		
c		write(25,'(I3)')
		ICT1=1
c           write(25,'(I5,I3,I3,I4,F7.1,F7.1,F7.1,F7.1,F7.1)')iyear,n_track,nn,ICT1,CLAT2(NN,1),CLON2(NN,1),CPRESS2(NN,1),CWIND2(NN,1),CDIRC2(NN,1)
		do ICT=1,NM1(NN)
           write(25,'(I5,I3,I3,I4,F7.1,F7.1,F7.1,F7.1,F7.1)')iyear,n_tra
     %ck,nn,ICT,CLAT2(NN,ICT),CLON2(NN,ICT),CPRESS2(NN,ICT),CWIND2(NN,IC
     %T),CDIRC2(NN,ICT)
		enddo
		do ITC=1,NM2(NN)
    		write(26,'(I5,I3,I3,I4,F7.1,F7.1)')iyear,n_track,nn,ITC,FLAT1(NN,I
     %TC),FLON1(NN,ITC)
		enddo
		   endif
             ENDDO
c------------------------
c       Etape 3 
c------------------------
c       On ecrit Position, Tourbillon, Pression, etc 
c       pour le dernier pas de temps de la trajectoire nn=ntr
		   if (ntr.ge.npt) then
 			 nnn=NTIME(ntr,IT)
		    if(DAY(nnn).lt.10)then
c          write(10,125)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(ntr)),abs(RLAT(ntr)),IPRESS(ntr),IVORT(ntr),ICIRC(ntr),IWIND(ntr),IDIRC(ntr)
	 write(16,125)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(ntr)),abs(R
     %LAT(ntr)),IPRESS(ntr),IWIND(ntr)
	 write(17,125)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(ntr)),abs(R
     %LAT(ntr)),IWIND(ntr),IDIRC(ntr)
          else
          write(16,126)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(ntr)),abs(
     %RLAT(ntr)),IPRESS(ntr),IWIND(ntr)
	 write(17,126)MONTH(nnn),DAY(nnn),HOUR(nnn),abs(RLON(ntr)),abs(R
     %LAT(ntr)),IWIND(ntr),IDIRC(ntr)
         endif
	 write(18,127)IYY,MONTH(nnn),DAY(nnn),HOUR(nnn),n_track,RLAT(ntr
     %),RLON(ntr),IPRESS(ntr),IVORT(ntr),ICIRC(ntr),IWIND(ntr),IDIRC(ntr
     %),IDIAM(ntr)
	      write(2,'(2I4,i6,i6,i6,i6,i6)') ILAT(Ntr),ILON(ntr),IPRESS
     %(ntr),IVORT(ntr),ICIRC(ntr),IWIND(ntr),IDIRC(ntr)
		   endif
	      ENDIF
	  END DO    
C just keep current tracks
	  ITT=0
	  DO IT=1,NTRACKS
	     IF(TUSED(IT)) THEN
	        ITT=ITT+1
	        NNT(ITT)=NNT(IT)
	        DO NN=1,NNT(IT)
	            CLAT(NN,ITT)=CLAT(NN,IT)
	            CLON(NN,ITT)=CLON(NN,IT)
	            CVORT(NN,ITT)=CVORT(NN,IT)
	  	    CPRESS(NN,ITT)=CPRESS(NN,IT)
	            CCIRC(NN,ITT)=CCIRC(NN,IT)
                    CWIND(NN,ITT)=CWIND(NN,IT)
                    CDIRC(NN,ITT)=CDIRC(NN,IT)
	            CU(NN,ITT)=CU(NN,IT)
	            CV(NN,ITT)=CV(NN,IT)
	            CLOSED(NN,ITT)=CLOSED(NN,IT)
	            NTIME(NN,ITT)=NTIME(NN,IT)
		    NUM3(NN,ITT)=NUM3(NN,IT)
		    NUM4(NN,ITT)=NUM4(NN,IT)
 		   DO ICT=1,NUM3(NN,ITT)		
		   CLAT1(NN,ITT,ICT)=CLAT1(NN,IT,ICT)		   	
		   CLON1(NN,ITT,ICT)=CLON1(NN,IT,ICT)
		   CPRESS1(NN,ITT,ICT)=CPRESS1(NN,IT,ICT)
		   CWIND1(NN,ITT,ICT)=CWIND1(NN,IT,ICT)
                   CDIRC1(NN,ITT,ICT)=CDIRC1(NN,IT,ICT)		
 		   ENDDO
		   DO ITC=1,NUM4(NN,ITT)		
		   FLAT1(NN,ITT,ITC)=FLAT1(NN,IT,ITC)		   	
		   FLON1(NN,ITT,ITC)=FLON1(NN,IT,ITC)		   	
 		   ENDDO
	        END DO
	     ENDIF
	  END DO
	  NTRACKS=ITT
 255      FORMAT(1X,A20,1X,F5.1,1X,F6.1,1X,F5.2,1X,A)
C
1000      CONTINUE
C
	  DO IC=1,N_NEW
	     IF(.NOT.CUSED(IC))THEN	
	         NTRACKS=NTRACKS+1
	         CLAT(1,NTRACKS)=TLAT(IC)
	         CLON(1,NTRACKS)=TLON(IC)
	         CVORT(1,NTRACKS)=TVORT(IC)
	         CPRESS(1,NTRACKS)=TPRESS(IC)
	         CCIRC(1,NTRACKS)=TCIRC(IC)
		 CWIND(1,NTRACKS)=WIND(IC)
                 CDIRC(1,NTRACKS)=DIRC(IC)		
	         CU(1,NTRACKS)=TU(IC)
	         CV(1,NTRACKS)=TV(IC)
	         CLOSED(1,NTRACKS)=TCLOSED(IC)
	         NTIME(1,NTRACKS)=N
	         NNT(NTRACKS)=1
		 NUM3(1,NTRACKS)=NUM1(IC)
		 NUM4(1,NTRACKS)=NUM2(IC)
		 do ICT=1,NUM3(1,NTRACKS)
		   CLAT1(1,NTRACKS,ICT)=RLATT(IC,ICT)
		   CLON1(1,NTRACKS,ICT)=RLONN(IC,ICT)
		   CPRESS1(1,NTRACKS,ICT)=RPRESS(IC,ICT)
		   CWIND1(1,NTRACKS,ICT)=RWINDD(IC,ICT)
                   CDIRC1(1,NTRACKS,ICT)=RDIRR(IC,ICT)			
		 enddo
		do ITC=1,NUM4(1,NTRACKS)
		   FLAT(1,NTRACKS,ITC)=BLAT(IC,ITC)
		   FLON(1,NTRACKS,ITC)=BLON(IC,ITC)		
		 enddo
	     ENDIF
	  END DO	
C =====================================================================
	END DO			
C	
 999	CONTINUE
	print*
	write(*,'(2(A13,I8))')'# CENTERs = ',N_CENTERS,'  # tracks = ',N
     %_TRACK
 5000	print*,'# nomber_great=',nggg,ddhh,n
c//////////////////////////////////////////////////////////////////////
c///////////////  Fin de la boucle temporelle  ////////////////////////
c//////////////////////////////////////////////////////////////////////
*
	ier=fstopc('MSGLVL','INFORM',.false.)
	ier = fstfrm(iun01)
	ier = fclos(iun01)
	if (nofst.ne.0) then
	   ier = fstfrm(iun51)
	   ier = fclos(iun51)
	endif
	CLOSE(UNIT=2)
*
	call message("Fin d'execution normale")
*
	end
c######################################################################
c#########################                  ###########################
c#########################   Sous-routines  ###########################
c#########################                  ###########################
c######################################################################
c***********************************************************************
c----  Sous-routine pour la lecture et la definition des options  ------
c
c       Attention: 
c       - l`algorithme fonction bien avec de fichiers PS rpn 
c       - il y a des problèmes avec des autres projections
c       - aux resolutions horizontales on peut avoir une dist_crit différente
c    
c***********************************************************************
      subroutine get_options(SILENT,IP3,GETVAR,INFILE,OUTFILE,NAJ,NPT,
	1    dist_crit,LEVEL,ISIGN2,FILTRE,PRCNT,PNMI,VC,ANAL,NOFST)
      implicit none
      integer silent, ip3, naj, npt, level, filtre
      integer anticycl,isign2,prcnt,pnmi,anal,nofst
      real vc,dist_crit
      character*2 getvar
      character*256 infile, outfile
      integer ier
      integer fstopc
      integer narg
      parameter (narg=16)
      character *8 cle(narg)
      character *256 def(narg),val(narg)
c- cles         
      data cle /'i.','o.','silent','ip3','getvar','naj','npt',
	1  'dist_crit','level','anticycl','filtre','prcnt','pnmi',
	2  'vc','anal','nofst'/
c- valeurs des cles si NON mentionnees a l'appel
      data val /'scrap','scrap','-1','-1',' ','-1','-1',
	1  '-1.','-1','-1','-1','-1','-1',
	2  '-1.','-1','-1'/
c- valeurs des cles si mentionnees a l'appel SANS valeurs
      data def /'scrap','scrap','0','0',' ','-1','8',
	1  '3.75','1000','0','800','0','0',
	2  '2.5','0','0'/
*-------------------------------------------------------------------
c On associe les parametres donnes a l'appel avec les 
c cles appropriees
      call ccard(cle,def,val,narg,-1)
*
      read (val(1),4) infile
      read (val(2),4) outfile
 4    format(a256)
      read (val(3),*) silent
      read (val(4),*) ip3
      read (val(5),3) getvar
 3    format(a2)
      read (val(6),*) naj
      read (val(7),*) npt
      read (val(8),*) dist_crit
      read (val(9),*) level
      read (val(10),*) anticycl
      read (val(11),*) filtre
      read (val(12),*) prcnt
      read (val(13),*) pnmi
      read (val(14),*) vc
      read (val(15),*) anal	
      read (val(16),*) nofst
      if (val(1).eq.'scrap') then
         print*
         print*,'******** Cles standards ********'
         print*
         print*,'-i : fichier standard d''entree'
         print*
         print*,'-o : fichier standard de sortie contenant les'
         print*,'     resultats du calcul'
      endif
*
c     On attribue les valeurs par defaut aux cles non-specifiees
*
      if (silent.eq.0) then
         ier=fstopc('MSGLVL','WARNIN',.false.)
      endif
*
      if (ip3.eq.-1) then
         read (def(4),*) ip3
      endif
*
      if (val(2).eq.'scrap') then
         print*, "Valeurs par defaut: output  = scrap"
      endif
*
      if (getvar.eq.' ') then
         getvar = 'GZ'
      endif
*
      if (naj.eq.-1) then
         naj = 8
      endif
*
      if (npt.eq.-1) then
	 if (naj.eq.8) then
	    npt = 8
	 elseif (naj.eq.4) then
	    npt = 4
	 else
	    print*
	    print*,'*** Erreur: Vous devez specifier une valeur pour NPT
     %'
	    print*
	    stop
	 endif
      endif
*
      if (dist_crit.eq.-1.) then
	 if (naj.eq.8) then
	    dist_crit = 3.
	 elseif (naj.eq.4) then
	    dist_crit = 7. 
	 else
	    print*
	    print*,'*** Erreur: Vous devez specifier une valeur pour DIS
     %T_CRIT'
	    print*
	    stop
	 endif
      endif
*
      if (vc.eq.-1.) then
	 vc = 2.5
      endif
*
      if (level.eq.-1) then
         level = 1000
      endif
*
      if (anticycl.eq.0) then
         isign2 = -1		
      else
	 isign2 = 1		
      endif
*
      if (filtre.eq.-1) then
         filtre = 800
      endif
*
	print*
	print*,'/////////////////////////////////////////////////////'
	print*,'>       Algorithme de calcul des trajectoires       <'
	print*,'/////////////////////////////////////////////////////'
	print*
	if (anal.lt.0) then
	   print*,' input: PREVISION'
	else
	   print*,' input: ANALYSES'
	   if (anal.eq.0) then
	      print*,' dateo = datev, ip2 = 0'
	   elseif (anal.eq.2) then
	      print*,' ip2 = heure a partir de la 1ere analyse'
	   else
	      print*,'ERREUR: anal = 0 (sans valeur) ou 2 ... Stop !'
	      stop
	   endif
	endif
	print*
	print*,'----------------------------------'
	print*,'-----       Parametres       -----'
	print*,'----------------------------------'
	print*,'  NAJ            : ',naj
	print*,'  NPT            : ',npt
	print*,'  Dist_crit (deg): ',dist_crit
	print*,'  SRAD (km)      : ',filtre
	print*,'  VC             : ',vc
	print*,'----------------------------------'
	print*,'  Niveau (hPa)   : ',level
	print*,'----------------------------------'
	print*
	if (isign2.eq.1) then
	   print*,'===>>> On suit les CYCLONES'
	   print*
	else
	   print*,'===>>> On suit les ANTICYCLONES'
	   print*
	endif
	if (prcnt.eq.-1) then
	   print*,'----------------------------------'
	   print*,'On suit les extremes de TOURBILLON'
	   print*,'----------------------------------'
	   print*
	else
	   print*,'----------------------------------'
	   print*,'On suit les extremes de PRESSION'
	   print*,'----------------------------------'
	   print*
	endif
	if (pnmi.eq.0) then
	   print*,'On lit la PNM dans le fichier d entree'
	   print*
	else
	   print*,'On calcul la PNM a partir de GZ(level)'
	   print*
	endif
	if (nofst.eq.0) then
	   print*,'**************************************'
	   print*,' PAS d ecriture en format FICHIER STD'
	   print*,'**************************************'
	   print*
	endif
      return
      end
