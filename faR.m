function [faR]=faR(s,M,Z,faR_struct)%function faR is the ratio of absorptive atomic potential to elastic potential 100fa*beta/fe in the reciprocal space
% global Ztable BKABS nthatom

Ztable=faR_struct.Ztable;
BKABS=faR_struct.BKABS;
% nthatom=faR_struct.nthatom;

faR=0;
i=1;
while(Z~=Ztable(i,2))
    i=i+1;
end
     nthatom=Ztable(i,1);
    GM=[0.05 0.15 0.30 0.70 1.30 2.00]; % labelled by i i0
    GMS2=[0 0.005 0.025 0.07 0.2 0.5 1.2 2.0 3.5 6.0];  %labelled by j j0
     
         MS2=M*s*s;

        i0=2; 
        while(M>=GM(i0))&&(i0<6)
            i0=i0+1;
        end
        
  if MS2>=2 %interpolate absorptive part using bi-cubic
         if i0<3 
             i0=3;
         end
         if i0==6 
             i0=5;
         end
         i0=i0-3;
         j0=6;
         
        RF=zeros(4); 
        FATM=zeros(4);
        
         for j=1:4
             for i=1:4
             RF(i)=BKABS(j0+j,i0+i,nthatom);
             end
           FATM(j)=CINT(GM(i0+1:i0+4),M,RF);
         end
     
        faR=CINT(GMS2(j0+1:j0+4),MS2,FATM);  
        
  end
        
  if MS2<2 % interpolate absorptive part using bi-quadratic
       COUNT=0;
       BFLAG=1;
       LFLAG=1;
       RFLAG=1;
       
        QF=zeros(4,4);% in the fortran code it is 5x4 matrix
        j0=2;
          while(MS2>=GMS2(j0))&&(j0<10)
            j0=j0+1;
           end
        i0=i0-3;
        j0=j0-3;
        
        for j=1:4
            j1=j0+j;
            if j1==0
                BFLAG=0; % reach the top of the table
            else
                for i=1:4
                   i1=i0+i;
                   if i1==0
                       LFLAG=0; % reach the left edge
                   elseif i1==7
                       RFLAG=0; % reach the right edge
                   else
                       QF(j,i)=BKABS(j1,i1,nthatom);
                   end
                end
            end
        end
            
        if (BFLAG==1)&&(LFLAG==1)
            faR=faR+BIQINT(M,MS2,GM(i0+1:i0+3),GMS2(j0+1:j0+3),QF(1:3,1:3)); % different from fortran
            COUNT=COUNT+1;
        end
        if (BFLAG==1)&&(RFLAG==1)
            faR=faR+BIQINT(M,MS2,GM(i0+2:i0+4),GMS2(j0+1:j0+3),QF(1:3,2:4));
            COUNT=COUNT+1;
        end
        if LFLAG==1
            faR=faR+BIQINT(M,MS2,GM(i0+1:i0+3),GMS2(j0+2:j0+4),QF(2:4,1:3));
            COUNT=COUNT+1;
        end
         if RFLAG==1
            faR=faR+BIQINT(M,MS2,GM(i0+2:i0+4),GMS2(j0+2:j0+4),QF(2:4,2:4));
            COUNT=COUNT+1;
         end
        faR=faR/COUNT;
  end
  
  if M<GM(1)
      faR=0;
  end
  
      
end



    