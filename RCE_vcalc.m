%RCE_vcalc.m

%Purpose: calculate the total column mass of water vapor and of dry air,
%then calculate the near-surface wind speed (background or gustiness)
%needed to achieve RCE (i.e. sfc heat fluxes perfectly balanced by constant
%radiative cooling)

%Created: 05-24-11, Dan Chavas

clear
clc

%% Calculate mass per unit area of dry air, water vapor

%%Inputs
dz_frac = .5;   %[]; desired output vertical resolution as fraction of sounding resolution (1.25 km)
dTdt_day = -1; %[K/day]; radiative cooling rate
U=3.88;    %[m/s]; background wind speed at surface; 5.3 m/s is from TOGA COARE
SST= 299.28;   %[K]; skin temp in model
C_T = 1.5e-3;    %[]; thermal exchange coefficient (constant)
C_q = 1.5e-3;    %[]; water vapor exchange coefficient (constant)

%%Constants (ha, good one)--taken directly from CM1 model
Rd = 287.04;   %[J/kg/K];  dry air gas constant
Rv = 461.5;   %[J/K/kg]
eps = Rd/Rv;
Cpd = 1005.7; %[J/kg/K]; spec heat of dry air
Cpv = 1875; %[J/kg/K]; spec heat of water vapor
Lv = 2501000;    %[J/kg]
g = 9.81;    %[m/s2]
p00   = 100000;   %[Pa]; reference pressure

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding_rotunno_emanuel_stratisothermal';
dir_in='/Users/drchavas/Documents/Research/Thesis/CM1/v15/analysis/';

%%load initial sounding variables
[zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc dz] = snd_extract(dir_in,snd_file,dz_frac);
%%Note: all in standard units (e.g. Pa, m, s etc.)

%% Calculate surface sensible heat flux, following how it's done in CM1
rho_sfc = rho00(1); %total air density at lowest model level
pi0s = (p_sfc/p00)^(Rd/Cpd);    %non-dim background surface pressure
pisfc = pi0s;   %just background value
th_s=SST/pisfc; %NOT th_sfc from sounding! this is how CM1 calculates surface th

thflux=C_T*U*(th_s-th00(1));    %[Km/s]; potential temperature flux
                                % pisfc = pi0s(i,j) + ppi(i,j,1) -- background surface pi + pipert at lowest model level
                                % th = th0(i,j,1)-tha(i,j,1) -- pot temp at lowest model level

F_sh = Cpd*rho_sfc*thflux  %[W/m2]

%% Calculate surface latent heat flux, following how it's done in CM1
%%calculate saturation vapor pressures at the surface
SST_C = SST-273.15;
es_sfc = 6.112*exp((17.67*SST_C)./(SST_C+243.5))*100;  %[Pa]

%%calculate saturation mixing ratios at the surface
qvs_sfc = eps*(es_sfc/(p_sfc-es_sfc));  %[kg/kg]

qvflux=C_q*U*(qvs_sfc-qv00(1));  %[(kg wv/kg air)*m/s]; water vapor flux

F_lh = Lv*rho_sfc*qvflux  %[W/m2]

F_sfc = F_sh+F_lh

%% Calculate radiative heat flux (out), following how it's done in CM1

%%METHOD 1: calculate from mixing ratio and dry air density in each layer in height coordinates
%%calculate T, Tv, rho, rh
[T rho rh the] = thermo(pp00,th00,qv00);
temp=pp00(T>200)';  %pressure levels where radiational cooling will occur
cool_frac = 1-(temp(end)/p_sfc);

%%calculate density of dry air
rho_d00 = pp00./(Rd*T00);

%%calculate mass per unit area of water vapor (assume water vapor = 0
%%above available levels)
m_v = dz*sum(rho_d00.*qv00);  %[kg/m2]

%%calculate total mass per unit area from sfc pressure
m_tot = p_sfc / g;  %[kg/m2]

%%calculate mass per unit area of dry air (misses very upper part of atmos)
%m_d = 1000*dz*sum(rho_d_prof);  %[kg/m2]

%%calculate mass per unit area of dry air as difference between m_tot and m_v
m_d = m_tot - m_v;  %[kg/m2]

%%remove fraction of dry air above tropopause that will not be cooled radiatively


%%calculate radiative heating rate
dTdt = dTdt_day / (60*60*24);   %[K/s]; radiative heating rate
F_rad = (Cpd*cool_frac*m_d+Cpv*m_v)*dTdt  %[W/m2]; only approximately lowest 900hPa are cooling!

%%solve for RCE gustiness velocity U_rce
U_rce = -F_rad/(Cpd*rho_sfc*C_T*(th_s-th00(1)) + Lv*rho_sfc*C_q*(qvs_sfc-qv00(1)))

%{
%%METHOD 2: calculate from mixing ratio in each pressure (i.e. mass) layer
%%calculate mass total mass per unit are of water vapor
pp00_2 = (pp00(1:end-1) + pp00(2:end))/2;
pp00_2 = [p_sfc;pp00_2];    %linearly interpolated pressure values at vertical grid midpoints
dp = (pp00_2(1:end-1)-pp00_2(2:end));   %layer pressure thicknesses
dm = dp/g; %mass of each layer of atmosphere
m_v_2 = dm'*qv00(1:length(dm))

%%calculate radiative heating rate using effective heat capacity (Emanuel "Atmospheric Convection" p. 110)
dTdt = dTdt_day / (60*60*24);   %[K/s]; radiative heating rate
F_rad_2 = Cpd*dTdt*(dm'*(1+.85*qv00(1:length(dm))))  %[W/m2]
%}


%}
