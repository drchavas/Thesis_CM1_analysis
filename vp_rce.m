%vp_rce.m

%Calculates potential intensity, v_p, for RCE with constant tropospheric
%cooling.

clc
clear

g=9.81;
Rd = 287.04;
Cp = 1004;
Cd = 1.5e-3;
T_tpp=150;
T_sst=305;
Q_cool = 1/86400;

u_sfc = 3 + .56;   %[m/s] includes lowest-level mean wind speed (CTRL = .521)
p_sfc = 101510;     %[Pa] from CTRL RCE sounding
p_tpp = 14827.5;    %[Pa] from CTRL RCE sounding
p0 = 100000;   %[Pa]
rho = 1.125;    %from CTRL RCE sounding

dp = p_sfc - p_tpp;

%account for fact that Q_cool is applied to potential temperature, and dtheta > dT for p<p0
th_fac = p0^(-Rd/Cp)*(quad('x.^.2859',p_sfc,p_tpp))/(p_tpp-p_sfc)  %Rd/Cp = .2859;

vp = sqrt(((T_sst-T_tpp)/T_tpp)*(Cp*th_fac*Q_cool*dp/(g*rho*Cd*u_sfc)))

