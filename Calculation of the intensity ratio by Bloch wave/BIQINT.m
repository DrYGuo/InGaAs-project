function BIQINT=BIQINT(M,MS2,GM,GMS2,F)%quadratic interpolation
   FATM=zeros(3);
   
   FATM(1)=QQINT(GM,M,F(1,:)); % GM is a 3x1 array
   FATM(2)=QQINT(GM,M,F(2,:));
   FATM(3)=QQINT(GM,M,F(3,:));
   
   BIQINT=QQINT(GMS2,MS2,FATM);  % GMS2 is a 3x1 array
end