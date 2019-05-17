%--------------------TILTE----------------------------------
%This program is written by Yueming Guo
% This program is used for calculating the ratios of integrated
% intensisites between (-1,11,1) and (1, 11, 1) for In(x)Ga(1-x)As at
% different composition x and thicknesses. This program calls Blochwave.m
% for the computation of the intensity ratios. To adjust the speed and
% accuracy of the Blochwave calculation, Sgmax in Blochwave.m is to be changed. 
% Most of the notations in this code have been adopted
%from the book by Zuo & Spence. The conventions for Sg,direction of zone axis follows his own thesis. Please see Fig.2.2  
% in Guo's PhD thesis(2017). In Chapter 2 of Guo's thesis,
%a more detailed derivation of Blochwave approach is given.
%For calculating elastic scattering, change the atomic displacements in the fifth column in the structure files,<u^2>(for Debye-waller factor),to 0. Do not leave it blank!
%You are welcomed to contact Yueming via his email address
%yueming.guo@monash.edu
clc;
clear;

structure=xlsread('GaAs.xlsx');% the name of the structure file It is desired to have the same atom species placed together rather than layer by layer
ABSratio=xlsread('absorptionratios.xlsx');
FPARAMS=xlsread('FPARAMS.xlsx');

                                                                                       


number_of_compositions=5;
step_of_compositions=0.01;
X_In_initial=0.48;

ratios=zeros(number_of_compositions,200);

ratio2s=zeros(number_of_compositions,200);

log_ratios=zeros(number_of_compositions,200);



for X_integer=1:number_of_compositions
    X_In=X_In_initial+(X_integer-1)*step_of_compositions;
    disp(X_In);
   %changing the lattice constants 
    lattice_constants=5.6535+(6.0583-5.6535)*X_In;
    structure(1,2)=lattice_constants;
    structure(1,3)=lattice_constants;
    structure(1,4)=lattice_constants;
   %The site occupancy for In 
    structure(2,6)=X_In;
    structure(3,6)=X_In;
    structure(4,6)=X_In;
    structure(5,6)=X_In;
   % The site occupancy for Ga 
    structure(6,6)=1-X_In;
    structure(7,6)=1-X_In;
    structure(8,6)=1-X_In;
    structure(9,6)=1-X_In;
   % The site occupancy for As in GaAs
    structure(10,6)=1-X_In;
    structure(11,6)=1-X_In;
    structure(12,6)=1-X_In;
    structure(13,6)=1-X_In;
   % The site occupancy for As in InAs
    structure(14,6)=X_In;
    structure(15,6)=X_In;
    structure(16,6)=X_In;
    structure(17,6)=X_In;
    disp(structure);
    [ratio,ratio2,log_ratio]=Blochwave(FPARAMS,ABSratio,structure);
    ratio_temp=ratio(:,1);
    ratio_temp=ratio_temp';
    
    ratio2_temp=ratio2(:,1);
    ratio2_temp=ratio_temp';
    
    log_ratio_temp=log_ratio(:,1);
    log_ratio_temp=log_ratio_temp';
    
    for i=1:200
    ratios(X_integer,i)=ratio_temp(i);
    end
    
     for i=1:200
    ratio2s(X_integer,i)=ratio2_temp(i);
     end
    
    for i=1:200
        log_ratios(X_integer,i)=log_ratio_temp(i);
    end
  
end