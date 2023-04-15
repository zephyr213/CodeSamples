function [knew,Cnew] = fitkc(T)
% find temperature dependent k and Cp 

% ABS data https://www.sciencedirect.com/science/article/pii/S0264127519308470?via%3Dihub#f0010
T0 = 0:50:250;
T0 = T0 + 273.15;
k0 = [.23,.25,.28,.29,.31,.33];
C0 = [780,1040,1490,1710,1865,2020];

% use linear interpolation as a test for now
knew = interp1(T0,k0,T);
Cnew = interp1(T0,C0,T);