%Domain_size.m

%Created: 10 Jan 2012, Dan Chavas

%This file is the same as TC_stats_plot.m, except that it only produces one
%plot with the radial wind profiles of 5 simulations with different domain sizes all on a single plot (i.e. a modified version
%of figure 5)


clear
clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

tmean0_usr = 70;    %[day]
tmeanf_usr = 100;    %[day]

T_mean = 2; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
dt_equil = 30;  %[days]; how long must be quasi-steady to define equilibrium
Cd_in = 1.5e-3; %only used to calculate r0_Lilly

%% Plots on/off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_gradient = 1;  %1=make plots with gradient wind
plot_full = 0;  %1=make plots with full wind (note: can do both)

%SINGLE/MULTI SIMULATION
plot_usrprof = 1;   %0=no; plot radial wind profile for user-defined time-means    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
   
%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
subdirs = {

%%DOMAIN SIZE
'CTRLv0qrhSATqdz5000_nx3072'
%}

}; %name of sub-directory with nc files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



numruns=length(subdirs);  %total number of runs you want to plot (taken from below)

%%want data for entire simulation
t0 = 0;
tf = 100;

%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=100000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=100;  %last z grid point [0,end]

%%Define subset of points for analysis
rmin_sub = 0;  %[km]; lowest value plotted
rmax_sub = 10000;    %[km]; highest value plotted
zmin_subsub = .5;  %[km]; lowest value plotted
zmax_subsub = 1; %[km]; highest value plotted

%%Append ax or 3d to subdir names for plot
for i=1:length(subdirs)
    if(run_types(i)==1) %ax
        subdirs_load{i}=sprintf('ax%s',subdirs{i});
    else    %3d
        subdirs_load{i}=sprintf('3d%s',subdirs{i});
    end
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;
for ss=1:numruns

    %%Load data for given simulation
    load(sprintf('../../../../CM1_postproc_data/simdata_Tmean2_dt30_dynamic/%s.mat',subdirs_load{ss}));
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    rrad_equil(ss) = rrad_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    rrad_equil_g(ss) = rrad_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    
    %timescales to those values
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    rrad_tau_equil(ss) = rrad_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    rrad_tau_equil_g(ss) = rrad_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
%    Vmax_gen(ss) = Vmax_gen_sim;
%    rmax_gen(ss) = rmax_gen_sim;
%    rrad_gen(ss) = rrad_gen_sim;
%    r0_gen(ss) = r0_gen_sim;
%    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    rrad_gen_g(ss) = rrad_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    rrad_tau_max_g(ss) = rrad_tau_max_g_sim; %defined using the GRADIENT wind
    rrad_max_g(ss) = rrad_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
%    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
%    Vmax_max(ss) = Vmax_max_sim;
%    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
%    rmax_max(ss) = rmax_max_sim;
%    rrad_tau_max(ss) = rrad_tau_max_sim; %defined using the GRADIENT wind
%    rrad_max(ss) = rrad_max_sim;
%    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
%    r0_max(ss) = r0_max_sim;
%    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
%    r0Lil_max(ss) = r0Lil_max_sim;
    
    %% Save data for all simulations
%    Vmax_movave_all(:,ss)=Vmax_movave_sim;
%    rmax_movave_all(:,ss)=rmax_movave_sim;
%    rrad_movave_all(:,ss)=rrad_movave_sim;
%    r0_movave_all(:,ss)=r0_movave_sim;
%    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;
    
    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    rrad_movave_g_all(:,ss)=rrad_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
%    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;


end


%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',24,'defaultaxesfontweight','bold','defaultlinelinewidth',6)

%Single/Multi simulation: Plot user-defined radial wind profiles

figure(1)
clf(1)

hold off
for i=1:numruns
    i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
%    xvals_pl = xvals_sub_all{i}(i_xvals_pl);
    xvals_pl = 2:4:1000;

    %% CALCULATE THEORETICAL PROFILE FROM Emanuel and Rotunno 2011
    CkCd = 1;
    V_theo = Vmax_equil_g(i)*(rmax_equil_g(i)./xvals_pl).*((2*(xvals_pl./rmax_equil_g(i)).^2)./(2-CkCd+CkCd*(xvals_pl./rmax_equil_g(i)).^2)).^(1/(2-CkCd)) - .5*fcor*1000*xvals_pl
    
%    plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl),pl_clrs{i})
    hold on
    plot(xvals_pl,V_theo,'b')
end
%plot(0*(1:xvals_pl(end)),'k')
input_title = sprintf('Blue: lh=1500 m; Black = ER11 analytic');
%title(input_title)
%xlabel('Radius [km]')
xlabel('Radius [km]')
%ylabel('Azimuthal gradient wind speed [m/s], z=1 km')
ylabel('Azimuthal wind speed [m/s]')
grid on
box on

plot(mpi*ones(1,xvals_pl(end)),'--','Color',[0 .5 0])
axis([0 600 0 100])