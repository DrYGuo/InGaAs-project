%--------------------TILTE----------------------------------
%This Blochwave programe is written by Yueming Guo
%Most of the notations in this code have been adopted
%from the book by Zuo & Spence. The conventions for Sg,direction of zone axis follows his own thesis. Please see Fig.2.2  
% in Guo's PhD thesis(2017). In Chapter 2 of Guo's thesis,
%a more detailed derivation of Blochwave approach is given.
%For calculating elastic scattering, change the atomic displacements in the fifth column in the structure files,<u^2>(for Debye-waller factor),to 0. Do not leave it blank!
%You are welcomed to contact Yueming via his email address
%yueming.guo@monash.edu
%% -------------Section 1----------------read Kirkland's fitting parameters for f(s) and structure files, then calculate the metric tensor
% first we define some global variables and read in the Kirklands'
% parameterisation files and the structure file
clc;
clear;
% global Metric volume structure FPARAMS K0 K E0 fsp Ztable number_of_species BKABS vcratio 
%read fitting parameters of atomic scattering factors
E0=200*1e3; % accelrating voltage
vcratio=sqrt(1-261127/(261127+(E0/1e3)*((E0/1e3)+1022)));  % v/c the ratio of the vecolocity to speed of light
FPARAMS=xlsread('FPARAMS.xlsx');
t=500;% thickness in Angstrom
%read the structure file(in the format of xlsx) into a matrix called "structure"
structure=xlsread('GaAs.xlsx');% the name of the structure file It is desired to have the same atom species placed together rather than layer by layer
ABSratio=xlsread('absorptionratios.xlsx');

%% -------------Section 2----------------read Kirkland's parameters for f(s) % this section and change f_s to f_soriginal to compare the speed with f_s and f_soriginal
%read Z numbers from structure
nthatom=1;%natom labels the nth atom species, so it starts from 1
Ztable=zeros(10,2); %assuming there are less than 10 elements in the crystal
Ztable(1,2)=structure(2,1);
Ztable(1,1)= nthatom;
[m,~]=size(structure);% get the size of the structure file

for i=3:m % start from the third row
    if structure(i,1)~=structure(i-1,1)
        nthatom=nthatom+1;
        Ztable(2,:)=[nthatom structure(i,1)];
    end
end
number_of_species=nthatom;

[m,~] = size(FPARAMS);%get the size 

for nthatom=1:number_of_species 
    Z=Ztable(nthatom,2);
for i=1:m  %read in Kirklands parameters
    if FPARAMS(i,1)== Z
        fsp.a(nthatom,1)=FPARAMS(i+1,1); % Kirkland's parameterasation of f(s)
        fsp.b(nthatom,1)=FPARAMS(i+1,2);
        fsp.a(nthatom,2)=FPARAMS(i+1,3);
        fsp.b(nthatom,2)=FPARAMS(i+1,4);
        fsp.a(nthatom,3)=FPARAMS(i+2,1);
        fsp.b(nthatom,3)=FPARAMS(i+2,2);
        fsp.c(nthatom,1)=FPARAMS(i+2,3);
        fsp.d(nthatom,1)=FPARAMS(i+2,4);
        fsp.c(nthatom,2)=FPARAMS(i+3,1);
        fsp.d(nthatom,2)=FPARAMS(i+3,2);
        fsp.c(nthatom,3)=FPARAMS(i+3,3);
        fsp.d(nthatom,3)=FPARAMS(i+3,4);
    end
end
end
   BKABS=zeros(10,6,number_of_species);% read ABSratio from Bird & King (1990)
   [nmax,~]=size(ABSratio);
   for nthatom=1:number_of_species 
      Z=Ztable(nthatom,2);
      n=1;
      while(Z~=ABSratio(n,2))
         n=n+11;
      end
          for i=1:10
              for j=1:6
          BKABS(i,j,nthatom)=ABSratio(n+1,j);
              end
              n=n+1;
          end
   end
%% --------------Section 3----assign lattice parameters from the defined matrix "structure"
a= structure(1,2);
b=structure(1,3);
c= structure(1,4);
alpha=structure(1,5)*pi/180;
beta=structure(1,6)*pi/180;
gamma=structure(1,7)*pi/180;
%the volume of the unit cell
volume=a*b*c*sqrt(1-(cos(alpha))^2-(cos(beta))^2-(cos(gamma))^2+2*cos(alpha)*cos(beta)*cos(gamma));
%the lattice parameters of the reciprocal unit cells
astar=b*c*sin(alpha)/volume;
bstar=c*a*sin(beta)/volume;
cstar=a*b*sin(gamma)/volume;
cos_alpha_star=(cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma));
cos_beta_star=(cos(gamma)*cos(alpha)-cos(beta))/(sin(gamma)*sin(alpha));
cos_gamma_star=(cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta));
 
% In order to calculate the length of the reciprocal lattice vector g
% we need to first calculate the metric tensor G* (see IUCr table vol.B
% page 4)

Metric =[astar^2 astar*bstar*cos_gamma_star astar*cstar*cos_beta_star; bstar*astar*cos_gamma_star bstar^2 bstar*cstar*cos_alpha_star; cstar*astar*cos_beta_star cstar*bstar*cos_alpha_star cstar^2];
MetricR=[a^2 a*b*cos(gamma) a*c*cos(beta); b*a*cos(gamma) b^2 b*c*cos(alpha); c*a*cos(beta) c*b*cos(alpha) c^2];%the metric tensor in real space

%% -----------------Section4-------------selecting the reflections that are close to the zone axis of choice
%Definition of constants
%hp=6.62607004e-34;  %planck's constant in international unit
%m=9.10938356e-31; %the mass of an electron in international unit
%e=1.60217662e-19; %the charge of an electron in international unit

lambda=12.2643/sqrt(E0*(1+0.978476*1e-6*E0)); % the relativiscally corrected wavlength
%U0=0.006648403*(1+1.956951*1e-6*E0)*V([0 0 0]);
f_s_struct.fsp=fsp;
f_s_struct.Ztable=Ztable;
f_s_struct.number_of_species=number_of_species;



faR_struct.Ztable=Ztable;
faR_struct.BKABS=BKABS;
faR_struct.nthatom=nthatom;

Ug_struct.Metric=Metric;
Ug_struct.volume=volume;
Ug_struct.structure=structure;
Ug_struct.E0=E0;
Ug_struct.vcratio=vcratio;
Ug_struct.f_s_struct=f_s_struct;
Ug_struct.faR_struct=faR_struct;

U0=Ug([0 0 0],Ug_struct);
K0=sqrt((1/lambda)^2+U0^2); %the magnitude of the wave vector inside the crystal
ZA=[0 0 1]; %zone axis;
SN=[0 0 1];%surface normal
Kt=[0 14.859 0]; %the projection of the wave vector Kt=Kt1*astar+Kt2*bstar+Kt3*cstar Kt=[Kt1,Kt2,Kt3]

%in JEMS Kt has the opposite sign
Kt0=Kt;% prepared for latter use
%% ------------------Section5--------Express the Zone Axis in reciprocal lattice:ZA=[u v w] is usually expressed in real space,i.e. the bases are a b c
% In order to calculate the zone axis component Kz of the incident wave K0,
% we need to transform ZA from real space into reciprocal space, ZAK=Ua*+Vb*+Wc* 

% This requires some knowledge of contravariant and covariant vectors. For more
% visit Wikipedia

ZAK= -(MetricR*ZA')';%the negative sign means that the zone axis is opposite to the incident beam direction
SNK= -(MetricR*SN')';% the negative sign means that the surface normal is against the incident beam

%this ZAK=[U V W], where U V W are the coefficients for a* b* c*
ZAK_length=sqrt(ZAK*Metric*ZAK');%calculate the length of ZAK
ZAK(1)=ZAK(1)*astar;
ZAK(2)=ZAK(2)*bstar;
ZAK(3)=ZAK(3)*cstar;
% the new ZAK=[ZAK1 A-1,ZAK2 A-1,ZAK3 A-1] the new coordinate has a 1A-1
% for each direction a* b* c*
Kz=sqrt(K0^2-Kt*Metric*Kt')*ZAK/ZAK_length; %Unit error! [Kz1 A-1,Kz2 A-1,Kz3 A-1] rather than [Kz1 a*, Kz2 b*, Kz3 c*]
Kz(1)=Kz(1)/astar;
Kz(2)=Kz(2)/bstar;
Kz(3)=Kz(3)/cstar;
% Now Kz=Kz1 a* + Kz2 b* + Kz3 c*, therefore Kz has the same coordinates as
% Kt and g (both are in a* b* c*)
K=Kt+Kz;
%excitation_error=10*Sg(g);% in the unit of nm-1
Sg_struct.K0=K0;
Sg_struct.Metric=Metric;
Sg_struct.K=K;
%% --------------Section6------------- beam selection: beams are selected based on the value of |Sg|
%PAY ATTENTION TO THE COMMENT after if 

% The eigenmatrix is created 
% The surface normal is antiparallel to the beam direction
NHOLZ=4; % number of HOLZ layers
count=0; 
hmax=40;
kmax=40;
lmax=40;
Sgmax=0.01; % the maximum excitation error
for hh=-hmax:hmax
    for kk=-kmax:kmax
        for ll=-lmax:lmax
            
            g=[hh kk ll];
         
            if (ZA*g'>=0)&&(ZA*g'<=NHOLZ) && (abs(Sg(g,Sg_struct))<=Sgmax)&&(mod(hh+kk,2)==0)&&(mod(hh+ll,2)==0 && abs(hh)+abs(kk)+abs(ll)>0)%&&(mod(hh+kk,2)==0)&&(mod(hh+ll,2)==0 is to include the hkl which are all even or all odd so this condition is not neccessary in general; && abs(hh)+abs(kk)+abs(ll)>0 is to exclude the central beam(0 0 0)
                
               count=count+1;
               checklist(count,1:3)=g;
               %checklist(count,4)=Ug(g);
               checklist(count,5)=Sg(g,Sg_struct);
               checklist(count,6)=abs(checklist(count,5));% avoid calling Sg function again
               checklist(count,7)=dot(SNK.*[astar bstar cstar],(K+g).*[astar bstar cstar])/(sqrt((K+g)*Metric*(K+g)')*sqrt(SNK*Metric*SNK'));
            end
        end
    end
end
B=sortrows(checklist,6);%sort the table, "checklist", based on |Sg|

%% -----------Section7----------------------- create an eigenmatrix
N=count+1;% size of the eigenmatrix
disp(N);
A=zeros(N,N);% creating a the eigenmatix with the 2*K0*Sg as the diagnaol elements and U(g-h) as the off-diagnal elements

K00=K0;%define a local variable for the purpose of parallel computing

B1=B(:,1);
B2=B(:,2);
B3=B(:,3);
B5=B(:,5);
B7=B(:,7);

for m=2:N
   mm=m-1; 
    hm=B1(mm);
      km=B2(mm);
      lm=B3(mm);
      gm=[hm km lm];
      A(m,1)=Ug(gm,Ug_struct)/B7(mm); % assign the first column of A
      A(1,m)=Ug(-gm,Ug_struct);% assign the first row of A, this makes it run faster
%       A(m,m)=2*K00*B5(mm)/B7(mm);% assing the excitation errors
      disp(m);     
    for n=2:N  % for n=2:m only half of the diagonal reduces the time      
      nn=n-1;
      hn=-B1(nn);
      kn=-B2(nn);
      ln=-B3(nn);
      gn=[hn kn ln]; 
     % A(1,n)=Ug(gn); % assing the first row of A
      
%       if m==n%&&(m~=1)
%         A(m,n)=2*K0*B(m-1,5);
%       end
       if m~=n 
    A(m,n)=Ug(gm+gn,Ug_struct)/B7(mm);
       end
    end
end
          
%% Call LAPACK for matrix diagonalisation Or just use the eig function
%
[eigve,eigva]= eig(A); %[eigenvector,eigenvalue] eigenvector is the matrix in which each column is the eigenvector of a certain eigenvalue, eigenva is the matrix where the diagonal elements are the eigenvalues
eiggamma=eigva/(2*K0); % the eiggamma is the matrix where the diagonal elements are the eigenvalues, denoted as "gamma" in the book by Zuo,Spence
% inv_eigve=inv(eigve); % the inverse matrix of the eigenvector matrix denoted as "C-1" in the book by Zuo & Spence
%% Calculate the output wavefunction and intensity for each disc
% incident plane wave (1 0 0 0 0 ...)
Plane=zeros(N,1); 
Plane(1)=1;
% nt is the number of thicknesses, dt is the stepsize of change in
% thickness
nt=100;
dt=20;
% "Intensity" is a table which stores the intensity at this excitation
% error for different reflections(in different rows)and different
% thicknesses(in different columns). The first three columns are the h k l
% indeces


%uncomment these lines for testing purpose only
% calculate the intensity at the centre of the disc 
% Intensity=zeros(N,nt+3);
% Intensity(1,1:3)=[0 0 0];
% Intensity(2:N,1:3)=B(:,1:3);
% for t=0:dt:dt*nt
%   exp_gamma=diag(exp(2*pi*1i*diag(eiggamma)*t));% create the matrix with the diagonal elements as the exponential of the gamma
%   out_wave=eigve*exp_gamma*(eigve\Plane);% out_wave=eigve*exp_gamma*inv_eigve*Plane; the inverse matrix inv(eigve) was avoided to following matlab suggestion
%   conj_out_wave=conj(out_wave);
%   out_intensity=conj_out_wave.*out_wave;
%  Intensity(:,4+t/dt)= out_intensity;
% end
% trange=0:dt:dt*nt;
% ratio=Intensity(4,4:nt+4)./Intensity(5,4:nt+4);
% plot(trange,Intensity(4:5,4:nt+4));
% 

%% ----------------Section8-------------------------correlate each pixel in a CBED disc with its Kt
%There are three steps:
%1)assign the lowest positive hkl as g1
%2)find a g2=g1 X Kz, where Kz and g1 have the basis of a* b* c*. This crossproduct will result a g2 in a b c. Then, we turn the basis of g2 into a* b* c*
%3)we rescale g2 to have the same length as g1 by g2*|g1|/|g2|. The user of
%the code will specify the sampling size in reciprocal space as #pixels for the
%radius of the largest overlapping disc

%First, we need to find a ZOLZ reflection that has the lowest indeces
%By assigning the fourth column in B' with |h|+|k|+|l|,we rank the reflections in order of indeces 
Bprime=B;
Bprime(:,4)=abs(Bprime(:,1))+abs(Bprime(:,2)+Bprime(:,3));
Bprime=sortrows(Bprime,4);

n=1; % reset n to 1
while Bprime(n,1)*ZA(1)+Bprime(n,2)*ZA(2)+Bprime(n,3)*ZA(3)~=0
    n=n+1;
end
%after this while loop,we find the ZOLZ hkl with the lowest index numbers
%we prefer to have g1 as 200 rather than -200,so we carry out the following
%operations
if Bprime(n,4)~=Bprime(n+1,4)
    g1=[Bprime(n,1) Bprime(n,2) Bprime(n,3)];
end

if (Bprime(n,4)==Bprime(n+1,4))&&(Bprime(n,1)+Bprime(n,2)+Bprime(n,3)>Bprime(n+1,1)+Bprime(n+1,2)+Bprime(n+1,3))
    g1=[Bprime(n,1) Bprime(n,2) Bprime(n,3)];
end

if (Bprime(n,4)==Bprime(n+1,4))&&(Bprime(n,1)+Bprime(n,2)+Bprime(n,3)<=Bprime(n+1,1)+Bprime(n+1,2)+Bprime(n+1,3))   
    g1=[Bprime(n+1,1) Bprime(n+1,2) Bprime(n+1,3)];
end
h1=g1(1);
k1=g1(2);
l1=g1(3);
% Second, we find the crossproduct of g1 and Kz
g2=[k1*Kz(3)-l1*Kz(2) l1*Kz(1)-h1*Kz(3) h1*Kz(2)-k1*Kz(1)]; % now g2 is in the basis of a b c, besides a scaling factor, 1/(c*.(a*xb*)) was missed out
g2=(MetricR*g2')';
% Third, we rescale g2 to the same length as g1
g2_length=sqrt(g2*Metric*g2'); % in A-1
g1_length=sqrt(g1*Metric*g1');
g2=g2*g1_length/g2_length; % g2=[h2 k2 l2] in the basis of a* b* c*, h2 k2 l2 may not be integers
%To now, we have got two basis vectors g1 and g2, which are perpendicular
%to each other and share the same length
%Here, we can get the increment of Kt in both x and y directions
np=50;% the number of pixels for the length, |g1|, is defined by the users,which is normally between 50-500,this sets the scale bar
dKtx=g1/np; % dKtx is a vector, which is the increment of Kt by shifting one pixel to the right  
dKty=g2/np; % dKty is a vector, which is the increment of Kt by shifting one pixel to the top
dKtx_length=g1_length/np;
%% ---------------------------Section 9-----------------------Now calculate the intensity Ig(x,y,t) where (x,y) are the integers that specify the coordinate in a CBED disc
semi_conv_angle=1; % in mrad
dKmax=semi_conv_angle*1e-3/lambda;% semi-convergence angle in A-1
Np_disc=round(dKmax/dKtx_length);% number of pixels for the radius of the disc

Ig=zeros(2*Np_disc+1,2*Np_disc+1,nt,N);
for y=-Np_disc:Np_disc
    for x=-Np_disc:Np_disc
        Kt=Kt0+x*dKtx+y*dKty; 
% ZAK= -(MetricR*ZA')';%the negative sign means that the zone axis is opposite to the incident beam direction
% SNK= -(MetricR*SN')';% the negative sign means that the surface normal is against the incident beam
% 
% %this ZAK=[U V W], where U V W are the coefficients for a* b* c*
% ZAK_length=sqrt(ZAK*Metric*ZAK');%calculate the length of ZAK        
% 
% ZAK_length=sqrt(ZAK*Metric*ZAK');%calculate the length of ZAK
% ZAK(1)=ZAK(1)*astar;
% ZAK(2)=ZAK(2)*bstar;
% ZAK(3)=ZAK(3)*cstar;
% the new ZAK=[ZAK1 A-1,ZAK2 A-1,ZAK3 A-1] the new coordinate has a 1A-1
% for each direction a* b* c*
Kz=sqrt(K0^2-Kt*Metric*Kt')*ZAK/ZAK_length; %Unit error! [Kz1 A-1,Kz2 A-1,Kz3 A-1] rather than [Kz1 a*, Kz2 b*, Kz3 c*]
Kz(1)=Kz(1)/astar;
Kz(2)=Kz(2)/bstar;
Kz(3)=Kz(3)/cstar;
% Now Kz=Kz1 a* + Kz2 b* + Kz3 c*, therefore Kz has the same coordinates as
% Kt and g (both are in a* b* c*)
K=Kt+Kz;

Sg_struct.K=K;

for m=2:N
    
    mm=m-1; 
    hm=B1(mm);
      km=B2(mm);
      lm=B3(mm);
      gm=[hm km lm];    
%     hm=B(m-1,1);
%       km=B(m-1,2);
%       lm=B(m-1,3);
%       gm=[hm km lm];     
        A(m,m)=2*K0*Sg(gm,Sg_struct)/B7(mm);   
end

% Call LAPACK for matrix diagonalisation Or just use the eig function
[eigve,eigva]= eig(A); %[eigenvector,eigenvalue] eigenvector is the matrix in which each column is the eigenvector of a certain eigenvalue, eigenva is the matrix where the diagonal elements are the eigenvalues
eiggamma=eigva/(2*K0); % the eiggamma is the matrix where the diagonal elements are the eigenvalues, denoted as "gamma" in the book by Zuo,Spence
% inv_eigve=inv(eigve); % the inverse matrix of the eigenvector matrix denoted as "C-1" in the book by Zuo & Spence

invC_plane=eigve\Plane; % calculate the inv(C)*(1 0 0 0...)


% "Intensity" is a table which stores the intensity at this excitation
% error for different reflections(in different rows)and different
% thicknesses(in different columns). The first three columns are the h k l
% indeces
Intensity=zeros(N,nt+3);
Intensity(1,1:3)=[0 0 0];
Intensity(2:N,1:3)=B(:,1:3);

for n=1:nt
  t=n*dt;
  exp_gamma=diag(exp(2*pi*1i*diag(eiggamma)*t));% create the matrix with the diagonal elements as the exponential of the gamma
  out_wave= eigve*(exp_gamma*invC_plane); % 1.invC_plane is thickness independent, so it is wise to avoid calculating inverse matrix in the loop 2.make use of the associative rule (AB)C=A(BC), the smaller size of (BC) reduces the calculation time for matrix multplication
  conj_out_wave=conj(out_wave);
  out_intensity=conj_out_wave.*out_wave;
 Intensity(:,4+n)= out_intensity;
 COL=x+Np_disc+1;
 ROW=Np_disc+1-y;
 
 Ig(ROW,COL,n,:)=Intensity(:,4+n);
 
end

XY=[x y];
disp(XY);


    end
    
    
end


%% ----------------------------Section 10-------------------------------------- image display run this section alone
% 
% I_1111=Ig(:,:,100,4);
% imshow(I_1111);
% 
% I1111=Ig(:,:,100,5);
% imshow(I1111)
I1_total=zeros(nt);
I2_total=zeros(nt);
ratio=zeros(nt);
for n=1:nt
    
for row=1:2*Np_disc+1
    for col=1:2*Np_disc+1
        I1_total(n)=I1_total(n)+Ig(row,col,n,4);
        I2_total(n)=I2_total(n)+Ig(row,col,n,5);
        ratio(n)=I1_total(n)/I2_total(n);
    end
end

end
thickness=2:2:2*nt;
plot(thickness,ratio);

% 
I000=Ig(:,:,100,1);
imshow(I000);




    
    
    
    









