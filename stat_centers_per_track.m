clear;
% Plotting interannual variability of "Indices".
model={'narr_uv850'}
rep_in0={'STORM_OUT_GUI/'}
rep_out0={'STORM_STAT_GUI/'}

mod=char(model(1))
rep_in=char(rep_in0(1))
rep_out=char(rep_out0(1))

    
  month={'01','02','03','04','05','06','07','08','09'};
 
s_year=2019
e_year=2019
year=s_year:e_year;
ny=e_year-s_year+1


for m=1:length(month);
    m
 for iy=1:ny;

year(iy)
dat1=load(strcat('/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_OUT_GUI/NARR/','out_',mod,'_',char(month(m)),int2str(year(iy)),'_gis.txt'));  
dat2=load(strcat('/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_STAT_GUI/NARR/','stat_',mod,'_',char(month(m)),int2str(year(iy)))); 
x=storm_speed(dat1);
             
 xx=x';
 
 
 
 for i=1:max(size(dat2))
 if  dat2(i,5) < 1.
     xx(i)=NaN ;
 end 
 end   
       y=cat(2,dat2(:,1:6),xx,dat2(:,8:12))    ;  
       
       
  clear dat1 dat2     
save(strcat('/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_STAT_GUI/NARR/','stat_',mod,'_',char(month(m)),int2str(year(iy)),'.2.txt'),'y','-ASCII');
          
        end
end
    
    





