close all, clear all, clc

model={'narr_uv850'}
rep_in0={'/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_CLIM_GUI/NARR/ASCII/'}
rep_out0={'/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_CLIM_GUI/NARR/NETCDF/'}

mod=char(model(1))
rep_in=char(rep_in0(1))
rep_out=char(rep_out0(1))


month={'janvier','fevrier','mars','avril','mai','juin','juillet','aout','septembre'};
%month={'octobre','novembre','decembre'};

%month2={'10','11','12'};
month2={'01','02','03','04','05','06','07','08','09'};

vari={'var','anom','anom_std','clim'};

start_year=1981;
end_year=2018;
%% --------------------


for ivar=1:length(vari);

ivar
for im=1:length(month);

im

if ivar==4 
end_year=start_year;
end


for iy=start_year:end_year 

iy

dir=strcat(rep_out,char(vari(ivar)),'/')
%dir=strcat('/BIG1/emmanuel/TEMPETES/OTHER_TRACKS/CanESM2_RCP45/',rep_out,char(vari(ivar)),'/')
cd(dir)
if ivar ==4
filenc=[char(strcat(char(vari(ivar)), '_storm_',mod,'_',char(month(im)),'.nc'))];
%filenc=[char(strcat(char(vari(ivar)), '_storm_ncep2_nh_',char(month(im)),'.2nc'))];
else
filenc=[char(strcat(char(vari(ivar)), '_storm_',mod,'_',char(month(im)),'_',int2str(iy),'.nc'))];
%filenc=[char(strcat(char(vari(ivar)), '_storm_ncep2_nh_',char(month(im)),'_',int2str(iy),'.2nc'))];

end
%nt=end_year-start_year+1;

% AMNO-GRID

lon=220:2.5:310;
lat=22.5:2.5:80;

% NH-GRID
%lon=0.:2.5:360.;
%lat=20.:2.5:90.;


nx=size(lon,2);
ny=size(lat,2);

time=1:1;
%time=1:nt;
[mlat,mlon]=meshgrid(lat,lon);
disp('Creer le fichier Netcdf...')

ncid = netcdf.create(filenc,'NC_WRITE');
% Definir dimensiones
dimid_lon = netcdf.defDim(ncid,'lon',nx);
dimid_lat = netcdf.defDim(ncid,'lat',ny); 
dimid_time = netcdf.defDim(ncid,'time',1); 
%dimid_time = netcdf.defDim(ncid,'time',nt); 
% Definir les attributions des variables
varid_lon = netcdf.defVar(ncid,'lon','double',dimid_lon);
netcdf.putAtt(ncid,varid_lon,'units','degrees_east')

%AMNO grid

netcdf.putAtt(ncid,varid_lon,'actual_range','220.f, 310.f')
%NH  grid
%netcdf.putAtt(ncid,varid_lon,'actual_range','0.f, 360.f')
netcdf.putAtt(ncid,varid_lon,'long_name','Longitude')
netcdf.putAtt(ncid,varid_lon,'axis','X')

% 
varid_lat = netcdf.defVar(ncid,'lat','double',dimid_lat);
%netcdf.putAtt(ncid,varid_lat,'units','degrees_north')
% amno grid
netcdf.putAtt(ncid,varid_lat,'actual_range','22.5f, 80.f')
% NH GRID
netcdf.putAtt(ncid,varid_lat,'long_name','Latitude')
netcdf.putAtt(ncid,varid_lat,'axis=','Y')
%


varid_time = netcdf.defVar(ncid,'time','double',dimid_time);
%netcdf.putAtt(ncid,varid_time,'units','Years since 1979-01-01 00:00:00')
netcdf.putAtt(ncid,varid_time,'units','months')
netcdf.putAtt(ncid,varid_time,'long_name','Time')
%netcdf.putAtt(ncid,varid_time,'actual_range','1979.,2011.')


if ivar ~= 3 
% 
%varid_1 = netcdf.defVar(ncid,'gen_dens','double',[dimid_lon,dimid_lat]);
varid_1 = netcdf.defVar(ncid,'gen_dens','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_1,'long_name','cyclogenesis_density')
netcdf.putAtt(ncid,varid_1,'units','nbre/month')
netcdf.putAtt(ncid,varid_1,'missing_value',-9999)

varid_2 = netcdf.defVar(ncid,'lys_dens','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_2,'long_name','cyclolysis_density')
netcdf.putAtt(ncid,varid_2,'units','nbre/month')
netcdf.putAtt(ncid,varid_2,'missing_value',-9999)

varid_3 = netcdf.defVar(ncid,'strm_den','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_3,'long_name','storm track density')
netcdf.putAtt(ncid,varid_3,'units','nbre/month')
netcdf.putAtt(ncid,varid_3,'missing_value',-9999)

varid_4 = netcdf.defVar(ncid,'strm_den_6','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_4,'long_name','intense_storm_density')
netcdf.putAtt(ncid,varid_4,'units','nbre/month')
netcdf.putAtt(ncid,varid_4,'missing_value',-9999)

varid_5 = netcdf.defVar(ncid,'mean_centers','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_5,'long_name','mean_centers_per_storm')
netcdf.putAtt(ncid,varid_5,'units','nbre/storm track')
netcdf.putAtt(ncid,varid_5,'missing_value',-9999)

varid_6 = netcdf.defVar(ncid,'strm_speed','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_6,'long_name','storm_speed')
netcdf.putAtt(ncid,varid_6,'units','Km/h')
netcdf.putAtt(ncid,varid_6,'missing_value',-9999)

varid_7 = netcdf.defVar(ncid,'circ','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_7,'long_name','mean circulation')
netcdf.putAtt(ncid,varid_7,'units','s^-1xm^2x10^7')
netcdf.putAtt(ncid,varid_7,'missing_value',-9999)

varid_8 = netcdf.defVar(ncid,'vort','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_8,'long_name','mean absolute vorticity')
netcdf.putAtt(ncid,varid_8,'units','s^-1x10^-5')
netcdf.putAtt(ncid,varid_8,'missing_value',-9999)

varid_9 = netcdf.defVar(ncid,'pres','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_9,'long_name','mean pressure')
netcdf.putAtt(ncid,varid_9,'units','hPa')
netcdf.putAtt(ncid,varid_9,'missing_value',-9999)

varid_10 = netcdf.defVar(ncid,'wnd_speed','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_10,'long_name','mean wind speed')
netcdf.putAtt(ncid,varid_10,'units','km/h')
netcdf.putAtt(ncid,varid_10,'missing_value',-9999)

else
    
varid_1 = netcdf.defVar(ncid,'gen_dens','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_1,'long_name','cyclogenesis_density')
netcdf.putAtt(ncid,varid_1,'units','standard deviation')
netcdf.putAtt(ncid,varid_1,'missing_value',-9999)

varid_2 = netcdf.defVar(ncid,'lys_dens','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_2,'long_name','cyclolysis_density')
netcdf.putAtt(ncid,varid_2,'units','standard deviation')
netcdf.putAtt(ncid,varid_2,'missing_value',-9999)

varid_3 = netcdf.defVar(ncid,'strm_den','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_3,'long_name','storm track density')
netcdf.putAtt(ncid,varid_3,'units','standard deviation')
netcdf.putAtt(ncid,varid_3,'missing_value',-9999)

varid_4 = netcdf.defVar(ncid,'strm_den_6','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_4,'long_name','intense_storm_density')
netcdf.putAtt(ncid,varid_4,'units','standard deviation')
netcdf.putAtt(ncid,varid_4,'missing_value',-9999)

varid_5 = netcdf.defVar(ncid,'mean_centers','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_5,'long_name','mean_centers_per_storm')
netcdf.putAtt(ncid,varid_5,'units','standard deviation')
netcdf.putAtt(ncid,varid_5,'missing_value',-9999)

varid_6 = netcdf.defVar(ncid,'strm_speed','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_6,'long_name','storm_speed')
netcdf.putAtt(ncid,varid_6,'units','standard deviation')
netcdf.putAtt(ncid,varid_6,'missing_value',-9999)

varid_7 = netcdf.defVar(ncid,'circ','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_7,'long_name','mean circulation')
netcdf.putAtt(ncid,varid_7,'units','standard deviation')
netcdf.putAtt(ncid,varid_7,'missing_value',-9999)

varid_8 = netcdf.defVar(ncid,'vort','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_8,'long_name','mean absolute vorticity')
netcdf.putAtt(ncid,varid_8,'units','standard deviation')
netcdf.putAtt(ncid,varid_8,'missing_value',-9999)

varid_9 = netcdf.defVar(ncid,'pres','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_9,'long_name','mean pressure')
netcdf.putAtt(ncid,varid_9,'units','standard deviation')
netcdf.putAtt(ncid,varid_9,'missing_value',-9999)

varid_10 = netcdf.defVar(ncid,'wnd_speed','double',[dimid_lon,dimid_lat]);
netcdf.putAtt(ncid,varid_10,'long_name','mean wind speed')
netcdf.putAtt(ncid,varid_10,'units','standard deviation')
netcdf.putAtt(ncid,varid_10,'missing_value',-9999)    
   
end


%varid_lon = netcdf.defVar(ncid,'lon','double',dimid_lon);

netcdf.endDef(ncid)

% % % Ecrire les variables
    netcdf.putVar(ncid,varid_lon,lon);
    netcdf.putVar(ncid,varid_lat,lat);
    netcdf.putVar(ncid,varid_time,time);
   
    
    if ivar ==4 
      dat=[char(strcat(rep_in,char(vari(ivar)),'/' ,char(vari(ivar)),'_storm_',mod,'_',char(month(im)),'.txt'))]
      %dat=[char(strcat('/BIG1/emmanuel/TEMPETES/OTHER_TRACKS/CanESM2_RCP45/',rep_in,char(vari(ivar)),'/' ,char(vari(ivar)),'_storm_',mod,'_',char(month(im)),'.txt'))]
    else   
      dat=[char(strcat(rep_in,char(vari(ivar)),'/' ,char(vari(ivar)),'_storm_',mod,'_',char(month(im)),'_',int2str(iy),'.txt'))]
      %dat=[char(strcat('/BIG1/emmanuel/TEMPETES/OTHER_TRACKS/CanESM2_RCP45/',rep_in,char(vari(ivar)),'/' ,char(vari(ivar)),'_storm_',mod,'_',char(month(im)),'_',int2str(iy),'.txt'))]
    end
    

     data=load(dat);
    data1=reshape(data(:,3),nx,ny,1);
    data2=reshape(data(:,4),nx,ny,1);
    data3=reshape(data(:,5),nx,ny,1);
    data4=reshape(data(:,6),nx,ny,1); 
    data5=reshape(data(:,7),nx,ny,1);
    data6=reshape(data(:,8),nx,ny,1);
    data7=reshape(data(:,9),nx,ny,1);
    data8=reshape(data(:,10),nx,ny,1);
    data9=reshape(data(:,11),nx,ny,1);
    data10=reshape(data(:,12),nx,ny,1);
%for mes=1:12    
%netcdf.putVar(ncid,varid_prec,[0 0 m-1],[nx ny 1],data(:,5));

    netcdf.putVar(ncid,varid_1,data1);
    netcdf.putVar(ncid,varid_2,data2);
    netcdf.putVar(ncid,varid_3,data3);
    netcdf.putVar(ncid,varid_4,data4);
    netcdf.putVar(ncid,varid_5,data5);
    netcdf.putVar(ncid,varid_6,data6);
    netcdf.putVar(ncid,varid_7,data7);
    netcdf.putVar(ncid,varid_8,data8);
    netcdf.putVar(ncid,varid_9,data9);
    netcdf.putVar(ncid,varid_10,data10);

%end
clear dir data data1 data2 data3 data4 data5 data6 data7 data8 data9 dir1 dir2

netcdf.close(ncid) 
%system('mv file.nc file1');
end

end

end
disp('C*est fini Manu ! Reveille toi');
