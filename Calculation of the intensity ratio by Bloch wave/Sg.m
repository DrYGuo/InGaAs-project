function [Sg] = Sg(g,Sg_struct)
%returns the excitation error for a reflection g
%   refer to Advanced TEM by Zuo & Spence page 58 equation(3.16)
% global  K0 Metric K 
K0=Sg_struct.K0;
Metric=Sg_struct.Metric;
K=Sg_struct.K;

    Sg = (K0^2-(K+g)*Metric*(K+g)')/(2*K0);

end

