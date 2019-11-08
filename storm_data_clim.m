clear;


model={'narr_uv850'}
rep_in0={'/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_STAT_GUI/NARR/'}
rep_out0={'/BIG1/emmanuel/TEMPETES/DATA_STORM/850_TRACKING/STORM_CLIM_GUI/NARR/ASCII/'}

mod=char(model(1))
rep_in=char(rep_in0(1))
rep_out=char(rep_out0(1))

main_dir=strcat(rep_in);      
save_dir1=strcat(rep_out,'/var/');
save_dir2=strcat(rep_out,'/clim/');
save_dir3=strcat(rep_out,'/anom/');
save_dir4=strcat(rep_out,'/anom_std/');


month={'janvier','fevrier','mars','avril','mai','juin','juillet','aout','septembre'};
%month={'octobre','novembre','decembre'};
%month2={'10','11','12'};
month2={'01','02','03','04','05','06','07','08','09'};


   start_year=1981;
   end_year=2019;
   ny=end_year-start_year+1
  
   year=[start_year:end_year];
   %% period clim
   y_start=1981;
   y_end=2019;

   iys=y_start-start_year+1;
   iye=y_end-start_year+1;
   
   vc=1.5
for m=1:length(month2);

     month(m)
   
   for iy=1:ny;
  yr=int2str(year(iy));
  
  

            T1=load(strcat(main_dir,'stat_',mod,'_',char(month2(m)),yr,'.2.txt'));
               for j=1:12   
          
               statn(:,j)=T1(:,j);
               end
               clear T1;%% period clim

                 
        for j=1:12;
        com1(j)=0;
        stat1(j)=0;
        end
      
           i1=0;

 % AMNO GRID  
           for i=1:4000;
                      
               if( statn(i,1)>=22.5 && statn(i,1)<=80);
               if( statn(i,2)>=220 && statn(i,2)<=310); 
 
 % NH GRID 
 %            for i=1:4205;
                   
%               if( statn(i,1)>=20. && statn(i,1)<=90);
%               if( statn(i,2)>=0. && statn(i,2)<=360);      
                   
                   
                   i1=i1+1;
  
             for j=1:12;      
               var(m,iy,i1,j)=statn(i,j);
             end
         
               end     
               end
           end   
   end
   
     siz=size(var,3);
    
     var2(:,:,:)=var(m,:,:,:);

% reverification que le cyclone atteint le minimum de tourbillon exige cad
% 1.5

    for iy=1:ny
         for i1=1:siz
              if var(m,iy,i1,10) < vc;
                       
                       for j=3:6       
                        var2(iy,i1,j)=0;
                        var(m,iy,i1,j)=0;
                       end
                       for j=7:12       
                        var2(iy,i1,j)=NaN;
                        var(m,iy,i1,j)=NaN;
                       end
                end
 
            end
      end
      
    
     for iy=1:ny
         for i1=1:siz
            for j=7:12       
                if var2(iy,i1,5)==0. ;
                var2(iy,i1,j)=NaN;
                end
                 if var(m,iy,i1,5)==0.  ;
                var(m,iy,i1,j)=NaN;
                end
 
            end
         end
     end   
         
        
    clim1(:,:)=nanmean(var2(iys:iye,:,:),1);
    st1(:,3:12)=nanstd(var2(:,:,3:12));
    st1(:,1)=var2(1,:,1);
    st1(:,2)=var2(1,:,2);
    for j=1:12
     if j==12
      clim(:,j)=clim1(:,j)*10.;
     st(:,j)=st1(:,j)*10.;    
     elseif j==8
      clim(:,j)=clim1(:,j)*2;   
      st(:,j)=st1(:,j);   
     else
      clim(:,j)=clim1(:,j);   
      st(:,j)=st1(:,j);   
     end
    end
    
     
for iy=1:length(year);
    yr=int2str(year(iy));
   
        for j=1:12;
            
            if j == 12
            var1(:,j)=var(m,iy,:,j)*10.;  
           
            elseif j==8      
            var1(:,j)=var(m,iy,:,j)*2.;
           
            else
             var1(:,j)=var(m,iy,:,j);
            end
           
          for i1=1:siz
            if j <3
            ano(i1,j)=clim(i1,j);
            ano_s(i1,j)=clim(i1,j);
            else   
            ano(i1,j)=var1(i1,j)-clim(i1,j);
            
             if st(i1,j)==0
              ano_s(i1,j)=NaN;
              else
              ano_s(i1,j)=ano(i1,j)/st(i1,j);
             end        
            
            end
          end  
        end    
        
   var1(isnan(var1))=-9999;
   dlmwrite(strcat(save_dir1,'var_storm_',mod, '_',char(month(m)),'_',yr,'.txt'), var1, 'delimiter', '\t','precision',4);
   
   ano(isnan(ano))=-9999;
   dlmwrite(strcat(save_dir3,'anom_storm_',mod,'_',char(month(m)),'_',yr,'.txt'), ano, 'delimiter', '\t','precision',4)

   ano_s(isnan(ano_s))=-9999;
   dlmwrite(strcat(save_dir4,'anom_std_storm_',mod, '_',char(month(m)),'_',yr,'.txt'), ano_s, 'delimiter', '\t','precision',4)
 
end

    clim(isnan(clim))=-9999;
     
    dlmwrite(strcat(save_dir2,'clim_storm_',mod,'_',char(month(m)),'.txt'), clim, 'delimiter', '\t','precision',4);

clear ano ano_s var var1 var2 clim clim1 statn

end




          
          
