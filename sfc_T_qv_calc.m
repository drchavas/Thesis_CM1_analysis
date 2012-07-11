%sfc_T_qv_calc.m

%Created: 05-05-11, Dan Chavas

clear
clc

%input
p=1015.1;   %[hPa]
SST=302.28  %[K]
dT_sfc = 2; %[K]; sfc thermal disequilibrium between air and water
RH_sfc = .78;   %sfc moisture disequilibrium between air and water

%constants
Rd = 287;   %[J/K/kg]
Cpd = 1004; %[J/K/kg]
Rv = 461;   %[J/K/kg]
eps = Rd/Rv;
p0 = 1000;


%equations
T_sfc = SST-dT_sfc; %[K]
SST_C = SST-273.15; %[C]
e_sat = 6.112*exp((17.67*SST_C)/(SST_C+243.5));  %[hPa]


%saturation
qv_sat = eps*(e_sat/(p-e_sat))*1000;  %[g/kg]

%th_sfc from T_sfc
th_sfc = T_sfc*(p0/p)^((Rd/Cpd))    %[K]

%qv_sfc from RH_sfc, e_sat
qv_sfc = RH_sfc*qv_sat