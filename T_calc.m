%T_calc.m

clear
clc

%input
p=1015.1;   %[hPa]
T=299.28;  %[K]


%constants
Rd = 287;   %[J/K/kg]
Cpd = 1004; %[J/K/kg]
Rv = 461;   %[J/K/kg]
eps = Rd/Rv;
p0 = 1000;


%equations
%T = th*(p/p0)^((Rd/Cpd))    %[K]
T=T-273.15; %[C]
e_sat = 6.112*exp((17.67*T)/(T+243.5))  %[hPa]

%saturation
qv_sat = eps*(e_sat/(p-e_sat))*1000  %[g/kg]