function [Ug,Vg,F_200]=Ug(g,Ug_struct)
% Now we calculate the structure factor Vg
% calculate the length of the reciprocal lattice vector

% global Metric volume structure E0 vcratio; 

Metric=Ug_struct.Metric;
volume=Ug_struct.volume;
structure=Ug_struct.structure;
E0=Ug_struct.E0; 
vcratio=Ug_struct.vcratio;
f_s_struct=Ug_struct.f_s_struct;
faR_struct=Ug_struct.faR_struct;
numberofrow=Ug_struct.numberofrow;

F200_struct.faR_struct=faR_struct;
F200_struct.vcratio=vcratio;

E00=200*1e3; % 200kV was used in the parameterization of the relativistically corrected atomic scattering factor in Rosenauer et al(2005) 
m0=9.10943e-31;
e=1.6022e-19;
c=299792458;
% 
% Ztable=Ug_struct.faR_struct.Ztable;
% number_of_species=Ug_struct.f_s_struct.number_of_species;

g_length=sqrt(g*Metric*g');
s=g_length/2; %the length of the scattering vector

% calculate the structure factors from atomic scattering factors
% [number_of_row,~]=size(structure);
Vg=0;
structure1=structure(:,1);
structure2=structure(:,2);
structure3=structure(:,3);
structure4=structure(:,4);
structure5=structure(:,5);
structure6=structure(:,6);


if abs(g(1))+abs(g(2))+abs(g(3))==2 % if g={2 0 0} this is only for GaAs and InGaAs
    occ=structure6(2);
    Vg=(478.77647/volume)*F200(occ,F200_struct,s)/(1+e*E00/(m0*c^2));% there is a difference between the units of F200(nm) and other Fg(A), so the scaling factor here is 10X different 
%     Vg=Vg/(1+e*E0/(m0*c^2));% The Vg in the line above was not properly scaled. The convention of fs and F in the paper by Rosenauer et al.(2005)PRB was different from that in the fs and F in the independent atom model by a scaling factor.
    F_200=F200(occ,F200_struct,s);
   
end


% f_total=zeros(3);
% for nthatom=1:number_of_species
%    Z=Ztable(nthatom,2);
%    f_total(nthatom)=f_s(s,Z,f_s_struct)*(1+1i*faR(s,M,Z,faR_struct)/(100*vcratio);
% end



Ztemp=zeros(numberofrow);
Mtemp=zeros(numberofrow);
ftotal=zeros(numberofrow);
rowfinal=0;
if abs(g(1))+abs(g(2))+abs(g(3))~=2 % if g={2 0 0}
  for row=2:numberofrow
r=[structure2(row); structure3(row); structure4(row)]; % the atom positions
M=8*pi^2*structure5(row);
Z=structure1(row);
Ztemp(row)=Z;
Mtemp(row)=M;

if row==2
    ftotal(row)=f_s(s,Z,f_s_struct)*(1+1i*faR(s,M,Z,faR_struct)/(100*vcratio));
end

if (row>=3)&&(Z==Ztemp(row-1))&&(M==Mtemp(row-1))
    ftotal(row)=ftotal(row-1);
end
if (row>=3)&&(M~=Mtemp(row-1))
    ftotal(row)=f_s(s,Z,f_s_struct)*(1+1i*faR(s,M,Z,faR_struct)/(100*vcratio));
end
    
occ=structure6(row);%occupancy
Vg=Vg+(47.877647/volume)*occ*ftotal(row)*exp(-2i*pi*g*r)*exp(-M*s^2); %(47.877647/volume)is the constant that converts Fg(Born) to Vg
   
% Vg=Vg+(47.877647/volume)*occ*f_s(s,Z,f_s_struct)*exp(-2i*pi*g*r)*exp(-M*s^2)*(1+1i*faR(s,M,Z,faR_struct)/(100*vcratio)); %(47.877647/volume)is the constant that converts Fg(Born) to Vg
  rowfinal=rowfinal+1;
  end
  
end

  Ug=0.006648403*(1+1.956951*1e-6*E0)*Vg; % relativistically "corrected"
end