%% Calculation of atomic scattering factors following Kirkland's
% parameterisation
function [fs]=f_s(s,Z,f_s_struct)
% global fsp Ztable number_of_species;
fsp=f_s_struct.fsp;
Ztable=f_s_struct.Ztable;
number_of_species=f_s_struct.number_of_species;

q=2*s; % for testing
fs=0;
for i=1:number_of_species
    if Z==Ztable(i,2)
    nthatom=Ztable(i,1);
    end
end  
    
 for n=1:3
     %fs=fs+a(n)/(q^2+b(n))+c(n)*exp(-d(n)*q^2);
     fs=fs+fsp.a(nthatom,n)/(q^2+fsp.b(nthatom,n))+fsp.c(nthatom,n)*exp(-fsp.d(nthatom,n)*q^2);
 end
 
end