%Domain_size.m

%Created: 10 Jan 2012, Dan Chavas

%This file is the same as TC_stats_plot.m, except that it only produces one
%plot with the radial wind profiles of 5 simulations with different domain sizes all on a single plot (i.e. a modified version
%of figure 5)

clear
clc

cd ../..

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
   
%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 1500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
subdirs = {

%%DOMAIN SIZE
'CTRLv0qrhSATqdz5000_nx192'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_nx768'
'CTRLv0qrhSATqdz5000_nx1536'
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
    load(sprintf('simdata/%s.mat',subdirs_load{ss}));
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    r12_equil(ss) = r12_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    r12_equil_g(ss) = r12_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    
    %timescales to those values
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    r12_tau_equil(ss) = r12_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    r12_tau_equil_g(ss) = r12_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
    Vmax_gen(ss) = Vmax_gen_sim;
    rmax_gen(ss) = rmax_gen_sim;
    r12_gen(ss) = r12_gen_sim;
    r0_gen(ss) = r0_gen_sim;
    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    r12_gen_g(ss) = r12_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    r12_tau_max_g(ss) = r12_tau_max_g_sim; %defined using the GRADIENT wind
    r12_max_g(ss) = r12_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
    Vmax_max(ss) = Vmax_max_sim;
    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
    rmax_max(ss) = rmax_max_sim;
    r12_tau_max(ss) = r12_tau_max_sim; %defined using the GRADIENT wind
    r12_max(ss) = r12_max_sim;
    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
    r0_max(ss) = r0_max_sim;
    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
    r0Lil_max(ss) = r0Lil_max_sim;
    
    %% Save data for all simulations
    Vmax_movave_all(:,ss)=Vmax_movave_sim;
    rmax_movave_all(:,ss)=rmax_movave_sim;
    r12_movave_all(:,ss)=r12_movave_sim;
    r0_movave_all(:,ss)=r0_movave_sim;
    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;
    
    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    r12_movave_g_all(:,ss)=r12_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;


end


%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%Single/Multi simulation: Plot user-defined radial wind profiles

h=figure(1)
clf(1)

set(h,'Position',[160 578 575 400])

ax1=axes('position',[0.15    0.15    0.75    0.75]);
axes(ax1)


hold off
for i=1:numruns
    i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
    xvals_pl = xvals_sub_all{i}(i_xvals_pl);

    plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl),pl_clrs{i})
    hold on
end
plot(0*(1:xvals_pl(end)),'k')
input_title = sprintf('Equilibrium radial profiles of gradient wind for variable domain width');
title(input_title)
xlabel('Radius [km]')
ylabel('Azimuthal wind speed [m/s]')
input_legend={'768 km','1536 km','3072 km','6144 km','12288 km'};
legend(input_legend)
set(ax1,'YTick',[-10 0 10 20 30 40 50 60 70 80],'XTick',[0 500 1000 1500])
grid on

cd Papers/RCE_equilibrium/
