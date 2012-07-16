%snd_drc_stratcnst.m

%Created: 28 Mar 2011, Dan Chavas

%Purpose: Generate a reasonable initial sounding for a given T_sst and T_tpp

clear
clc
clf

%% USER INPUT %%%%%%%%%%%%%%%%%%
T_sfc_out = 287.50;    %[K]; T_sfc = T_sfc_out - 3
T_tpp = 200;   %[K]; temp set constant at and above height of this temp

p_sfc = 101500; %[Pa]
z_bl = 1500;    %[m]
rh_bl = .8;
rh_ft = .2;
rh_strat=.01;   %relative humidity above tropopause

dz = 625;  %[m]; desired vertical resolution
z_top = 39989;  %[m]

save_output_sounding=1; %0=don't save; o.w. = save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Constants
Rd=287;  %[J/kg/K]
Rv=461.5;   %[J/K/kg]
Lv=2.5e6;   %[J/kg]
Cpd=1005.7; %[J/kg/K]; spec heat of dry air
eps=Rd/Rv;
g=9.81; %[m/s2]
p0 = 100000; %[Pa]

nz_sub = z_top/dz + 1;

%% Calculate surface values of theta, qv
%%calculate saturation water vapor mixing ratio at surface
T_sfc = T_sfc_out - 5;  %because if the initial sounding is too moist it will be very stable and will take forever before convection can actually kick in!
T_sfc_C = T_sfc - 273.15;
es_sfc = 6.112*exp((17.67*T_sfc_C)./(T_sfc_C+243.5))*100;  %[Pa]
qvs_sfc = eps*(es_sfc./(p_sfc-es_sfc));  %[kg/kg]
qv_sfc = rh_bl*qvs_sfc

th_sfc = T_sfc*(p0/p_sfc)^(Rd/Cpd)
%the_sfc = th_sfc*exp(Lv*qv_sfc/(Cpd*T_sfc))

%% Create vertical profile given suface values
zz00 = dz/2:dz:z_top+dz;
gam_d = 9.8/1000;    %[K/m]

%% Calculate temperature profile by iterating upwards
gam = gam_d*(1+Lv*qv_sfc/(Rd*T_sfc))/(1+(eps*Lv^2*qv_sfc)/(Cpd*Rd*T_sfc^2));    %[K/km], approx moist adiabatic lapse rate
T00(1) = T_sfc + gam*(zz00(1)-0);

T00_C = T00(1) - 273.15;
es = 6.112*exp((17.67*T00_C)./(T00_C+243.5))*100;  %[Pa]
Tv = T_sfc*(1+qv_sfc/eps)/(1+qv_sfc);
p(1) = p_sfc - (p_sfc/(Rd*Tv))*(zz00(1)-0);
qvs00 = eps*(es./(p(1)-es));  %[kg/kg]
qv00(1) = rh_bl*qvs00;

th00(1) = T00(1)*(p0/p)^(Rd/Cpd);

for i=2:length(zz00)
    gam = gam_d*(1+Lv*qvs00/(Rd*T00(i-1)))/(1+(eps*Lv^2*qvs00)/(Cpd*Rd*T00(i-1)^2));    %[K/km], approx moist adiabatic lapse rate
    T00(i) = T00(i-1) - gam*(zz00(i)-zz00(i-1))
    if(T00(i)<T_tpp)
        T00(i)=T_tpp;
    end
    
    T00_C = T00(i) - 273.15;
    es = 6.112*exp((17.67*T00_C)./(T00_C+243.5))*100;  %[Pa]
    Tv = T00(i)*(1+qv00(i-1)/eps)/(1+qv00(i-1));
    p(i) = p(i-1) - g*(p(i-1)/(Rd*Tv))*(zz00(i)-zz00(i-1)) %presure
    qvs00 = eps*(es./(p(i)-es));  %[kg/kg] 
    if(zz00(i)<=z_bl)
        qv00(i) = rh_bl*qvs00;
    else
        qv00(i) = rh_ft*qvs00;
    end
    
    th00(i) = T00(i)*(p0/p(i))^(Rd/Cpd);
end

[zz00'/1000 T00' th00' qv00']

%% Save sounding file
u00 = 0*zz00;
v00 = 0*zz00;
if(save_output_sounding~=0)

    %%sounding header: %1) p_sfc (mb); 2) th_sfc (K); 3) qv_sfc (g/kg)
    %%sounding columns: 1) zz00 (m); 2) th00_RCE (K); 3) qv00_RCE (g/kg); 4) u00 (m/s); 5) v00 (m/s)
    sounding_all=[zz00' th00' 1000*qv00' u00' v00'];
    sounding_all=double(sounding_all);

    file_out = sprintf('input_soundings/stratcnst/input_sounding_stratcnstT%iK_T%i.00',T_tpp,T_sfc_out);
%    save(file_out,'p_sfc','th_sfc','qv_sfc','sounding_all');    
    clear fid
    fid = fopen(file_out,'w');
    fprintf(fid,'   %6.2f     %6.2f     %6.2f\n',p_sfc/100,th_sfc,1000*qv_sfc);
    fprintf(fid,'   %8.5f     %8.5f     %8.5f     %8.5f    %8.5f\n',sounding_all');
    
end

%% Plot sounding
figure(1)
subplot(1,2,1)
hold on
plot(1000*[qv00],[zz00/1000],'r','LineWidth',2)

    set(gca,'fontweight','bold','fontsize',11)
    axis([0 max(1000*qv00) 0 z_top/1000])
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
    axis([125 800 0 z_top/1000])
    input_xlabel=sprintf('T and theta [K]');
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_title=strrep(sprintf('T and theta [K]'),'_','\_');
    title(input_title)

%}