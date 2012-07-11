%snd_isothermal_calc.m

%Created: 28 Mar 2011, Dan Chavas

%Purpose: initial estimate of RCE profile, just to start the simulation at
%something near its final equilibrium state.

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc
clf

%% USER INPUT %%%%%%%%%%%%%%%%%%
T_isotherm = 200;   %[K]; temp set constant at and above height of this temp
dz = .625;  %desired vertical resolution
z_top = 30000;  %[m]
dTdz=-10.0/1000;    %[K/m]
Tv_bar=235;  %[K] -- free parameter :) can match pressures at z>~20 km
rh_strat=.01;   %relative humidity above tropopause
save_output_sounding=0; %0=don't save; o.w. = save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
Rd=287;  %[J/kg/K]
Rv=461.5;   %[J/K/kg]
Cpd=1005.7; %[J/kg/K]; spec heat of dry air
epsilon=Rd/Rv;
g=9.81; %[m/s2]

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding_rotunno_emanuel';

nz_sub = (z_top/1000)/dz;

%%load initial sounding variables
dir_full = '/Users/drchavas/Documents/Research/Thesis/CM1/v15/analysis';
[zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);
zz00(end)


%% Extend sounding to arbitrary height
%{
%%1) z
dz=zz00(2)-zz00(1); %[m]
nz=floor(z_top/dz);
zz00=(1:nz)'*dz - dz/2;

%%2) T
i_T220 = find(T00<220,1);
temp=1:(length(T00)-i_T220+1);
T00(i_T220:end)=T00(i_T220-1)+dTdz*dz*temp;
clear temp
temp=1:(nz-length(T00));
T00(end+1:nz)=T00(end)+dTdz*dz*temp;    %FIGURE ME OUT!
T00=T00';

%%3) p -- assume exponential (from hydrostatic)
H=Rd*Tv_bar/g;  %[m]
p_exp = p_sfc*exp(-zz00/H);
pp00(end+1:nz)=0;
pp00(zz00>20000)=p_exp(zz00>20000);
pp00=pp00';
pi00=(pp00/100000).^(Rd/Cpd);

%%4) th
th00=T00./pi00;

%%5) rh (just add zeroes at the top)
rh00(end+1:length(zz00))=rh00(end);
rh00=rh00';

%%6) u and v (just add zeroes at the top)
u00(end+1:length(zz00))=u00(end);
u00=u00';
v00(end+1:length(zz00))=v00(end);
v00=v00';
%}

%compare to isothermal atmospheric (i.e. exponential pressure)
H=Rd*Tv_bar/g;  %[m]
p_exp = p_sfc*exp(-zz00/H);


%% Adjust sounding to be isothermal at heights with T<T_isotherm
indices=find(T00<T_isotherm);
th00(indices)=T_isotherm./pi00(indices);
T00(indices)=T_isotherm;
rh00(indices)=rh_strat;

%%%Keep relative humidity constant
%%Updated saturation vapor pressures
T00_C = T00-273.15;
es00 = 6.112*exp((17.67*T00_C)./(T00_C+243.5))*100;  %[Pa]

%%Updated saturation mixing ratios
qvs00 = epsilon*(es00./(pp00-es00));  %[kg/kg]

%%Updated mixing ratios keeping relh constant
qv00 = rh00.*qvs00;

%% Save sounding file
if(save_output_sounding~=0)

    %%sounding header: %1) p_sfc (mb); 2) th_sfc (K); 3) qv_sfc (g/kg)
    %%sounding columns: 1) zz00 (m); 2) th00_RCE (K); 3) qv00_RCE (g/kg); 4) u00 (m/s); 5) v00 (m/s)
    sounding_all=[zz00/2 th00 1000*qv00 u00 v00];
    sounding_all=double(sounding_all);

    file_out = sprintf('input_sounding_stratcnstT%iK',T_isotherm);
%    save(file_out,'p_sfc','th_sfc','qv_sfc','sounding_all');    
    clear fid
    fid = fopen(file_out,'w');
    fprintf(fid,'   %6.2f     %6.2f     %6.2f\n',p_sfc/100,th_sfc,1000*qv_sfc);
    fprintf(fid,'   %8.5f     %8.5f     %8.5f     %8.5f    %8.5f\n',sounding_all');
    
end


%% Plot p_exp for comparison
figure(3)
plot(pp00/100,zz00/1000)
hold on
plot(p_exp/100,zz00/1000,'r')

%% Plot sounding
figure(1)
subplot(1,2,1)
hold on
plot(1000*[qv00],[zz00/1000],'r','LineWidth',2)

    set(gca,'fontweight','bold','fontsize',11)
    axis([0 20 0 z_top/1000])
    input_xlabel=sprintf('water vapor mixing ratio [g/kg]');
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_title=strrep(sprintf('water vapor mixing ratio [g/kg]'),'_','\_');
    title(input_title)

figure(1)
subplot(1,2,2)
hold on
plot([th00],[zz00/1000],'r','LineWidth',2)
plot([T00],[zz00/1000],'b','LineWidth',2)

    set(gca,'fontweight','bold','fontsize',11)
    axis([150 450 0 z_top/1000])
    input_xlabel=sprintf('T and theta [K]');
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_title=strrep(sprintf('T and theta [K]'),'_','\_');
    title(input_title)

%}