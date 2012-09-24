%vp_rce.m

%Calculates potential intensity, v_p, for RCE with constant tropospheric
%cooling.

clc
clear

c_CM1 = constants_CM1(); %c_CM1: [g rd cp cv p00 xlv]

g=c_CM1(1); %[m/s2]
Rd=c_CM1(2);  %[J/kg/K]
Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
Rv=c_CM1(4);   %[J/K/kg]
p0 = c_CM1(5); %[Pa]
Lv=c_CM1(6);   %[J/kg]

eps=Rd/Rv;

Cd = 1.5e-3;    %[-]
T_tpp=200;  %[K]
T_sst=300;  %[K]
Q_cool = 1/86400;   %[K/s]

u_sfc = 1.5 + .76;   %[m/s] includes lowest-level mean wind speed (CTRL = .521)
p_sfc = 101500;     %[Pa] from CTRL RCE sounding
p_tpp = 15104.5;    %[Pa] from CTRL RCE sounding
p0 = 100000;   %[Pa]
rho = 1.125;    %from CTRL RCE sounding

dp = p_sfc - p_tpp;

%account for fact that Q_cool is applied to potential temperature, and dtheta > dT for p<p0
th_fac = p0^(-Rd/Cpd)*(quad('x.^.2859',p_sfc,p_tpp))/(p_tpp-p_sfc)  %Rd/Cpd = .2859;

vp = sqrt(((T_sst-T_tpp)/T_tpp)*(Cpd*th_fac*Q_cool*dp/(g*rho*Cd*u_sfc)))

