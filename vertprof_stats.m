%vertprof_stats.m

%Created: 11 Sep 2012, Dan Chavas

%Purpose: Output various relevant statistics given an input sounding

%%TESTING
%clc
%clear
%load '/Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/temp.mat'
%%%%%
function [zz_mid dp_mid dthdz_mid N_v_mid H_tpp N_v_trop] = vertprof_stats(th00,qv00,thv00,T00,pp00,zz00,Ttpp,z_trop_bot); 

c_CM1 = constants_CM1(); %c_CM1: [g rd cp cv p00 xlv]

g=c_CM1(1); %[m/s2]
Rd=c_CM1(2);  %[J/kg/K]
Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
Rv=c_CM1(4);   %[J/K/kg]
p0 = c_CM1(5); %[Pa]
Lv=c_CM1(6);   %[J/kg]


% Note: the following is calculated at midpoints of levels, and thus n=nz-1
zz_mid = mean([zz00(2:end)' zz00(1:end-1)'],2);
dz_mid = zz00(2:end)-zz00(1:end-1);
dp_mid = pp00(1:end-1)-pp00(2:end);

%% dry static stability: dtheta/dz
dth_mid = th00(2:end)-th00(1:end-1);
dthdz_mid = dth_mid./dz_mid;

%% buoyancy (Brunt-Vaisala) frequency: N = sqrt((g/th_v)*d(th_v)/dz) (p.166 Emanuel)
dthv_mid = thv00(2:end)-thv00(1:end-1);
dthvdz_mid = dthv_mid./dz_mid;
thv_mid = mean([thv00(2:end)' thv00(1:end-1)'],2);
N_v_mid = sqrt((g./thv_mid').*dthvdz_mid);

%% Height of tropopause -- taken as first point where T<Ttpp moving up from below
i_ztpp_low = find(T00<Ttpp,1)-1;
H_tpp = zz00(i_ztpp_low)+dz_mid(i_ztpp_low)*(T00(i_ztpp_low)-Ttpp)/(T00(i_ztpp_low)-T00(i_ztpp_low+1)); %[m]

%% Tropospheric mean buoyancy frequency (above z_trop_bot (e.g. boundary layer top))
i_trop = zz_mid<=H_tpp & zz_mid>=z_trop_bot;   %indices for desired levels
N_v_trop = (dp_mid(i_trop)*sqrt((g./thv_mid(i_trop)).*dthdz_mid(i_trop)'))/sum(dp_mid(i_trop));

end