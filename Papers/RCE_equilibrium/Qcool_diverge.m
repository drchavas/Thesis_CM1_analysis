%Qcool_diverge.m

%Created: 11 June 2012, Dan Chavas

%Purpose: This file plots the normalized equilibrium radial gradient wind profiles at
%the BL top for varying radiative cooling rate to demonstrate at what
%radius varying the radiative cooling rate becomes important.

clear
clc
figure(1)
clf(1)

cd ../..

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

t0 = 0;

T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
equil_dynamic = 0;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad
Cd_in = 1.5e-3; %only used to calculate r0_Lilly

%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 1500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
sim_sets_all = {'Qcool'};  %name out output subdir (within simsets_Tmean#/PLOTS/[sim_set]/) where plots will be saved

for jj = 1:length(sim_sets_all)

sim_set = sim_sets_all{jj};
    
switch sim_set

    case 'Qcool'
        CTRL_val = 1; %CTRL value of quantity varied across simulations
        units = 'K day^{-1}';
        multipliers = [-2 -1 0 1];
        subdirs_set = {
            %%Qcool
            'CTRLv0qrhSATqdz5000_nx3072_rad0.125K'
            'CTRLv0qrhSATqdz5000_nx3072_rad0.25K'
            'CTRLv0qrhSATqdz5000_nx3072'
            'CTRLv0qrhSATqdz5000_nx3072_rad1.0K'
            %'CTRLv0qrhSATqdz5000_nx3072_rad2.0K'
        }
    case 'QcoolVpcnst'
        CTRL_val = 1; %CTRL value of quantity varied across simulations
        units = 'K day^{-1}';
        multipliers = [-1 0 1 2 3 4 5];
        subdirs_set = {
            %%Qcool at constant Vp (can't get Vp ~ 92 m/s with rad0.125)
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K'
            'CTRLv0qrhSATqdz5000_nx3072'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K'
        }
    case 'QcoolVplvHcnst'
        CTRL_val = 1; %CTRL value of quantity varied across simulations
        units = 'K day^{-1}';
        %multipliers = [-1 0 1 2 3 4 5];
        multipliers = [-1 0 1 2];
        subdirs_set = {
            %%Q_rad at constant Vp, lv/H constant (can't get Vp ~ 92 m/s with rad0.125)
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K_lv121'
            'CTRLv0qrhSATqdz5000_nx3072'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K_lv81'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K_lv60'
            %'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K_lv44'
            %'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K_lv26'
            %'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K_lv19'
        }
end
        
    
%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numruns=length(subdirs_set);  %total number of runs you want to plot (taken from below)

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
for i=1:length(subdirs_set)
    if(run_types(i)==1) %ax
        subdirs_load{i}=sprintf('ax%s',subdirs_set{i});
    else    %3d
        subdirs_load{i}=sprintf('3d%s',subdirs_set{i});
    end
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;
clear Vmax_tau_equil rmax_tau_equil rmid_tau_equil r0_tau_equil r0Lil_tau_equil Vmax_tau_equil_g rmax_tau_equil_g rmid_tau_equil_g r0_tau_equil_g r0Lil_tau_equil_g
clear tau_gen Vmax_gen_g rmax_gen_g rmid_gen_g r0Lil_gen_g Vmax_tau_max_g Vmax_max_g rmax_tau_max_g rmax_max_g r0_tau_max_g r0_max_g r0Lil_tau_max_g r0Lil_max_g
clear Vmax_equil rmax_equil rmid_equil r0_equil r0Lil_equil Vmax_equil_g rmax_equil_g rmid_equil_g r0_equil_g r0Lil_equil_g
clear Vmax_movave_g_all rmax_movave_g_all r0_movave_g_all r0Lil_movave_g_all
clear xvals_sub_all data_tmean_usr_g_all
clear mpi_all
for ss=1:numruns

    %%Load data for given simulation
    if(equil_dynamic==1)
        if(wrad_const == 1)
            load(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst/%s.mat',T_mean,dt_final_dynamic,subdirs_load{ss}));
        else
            load(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic/%s.mat',T_mean,dt_final_dynamic,subdirs_load{ss}));
        end
    else
        if(wrad_const == 1)
            load(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst/%s.mat',T_mean,tf-dt_final,tf,subdirs_load{ss}));
        else
            load(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/%s.mat',T_mean,tf-dt_final,tf,subdirs_load{ss}));
        end
    end
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    rmid_equil(ss) = rmid_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    rmid_equil_g(ss) = rmid_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    
    %timescales to those values
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    rmid_tau_equil(ss) = rmid_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    rmid_tau_equil_g(ss) = rmid_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
%    Vmax_gen(ss) = Vmax_gen_sim;
%    rmax_gen(ss) = rmax_gen_sim;
%    rmid_gen(ss) = rmid_gen_sim;
%    r0_gen(ss) = r0_gen_sim;
%    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    rmid_gen_g(ss) = rmid_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    rmid_tau_max_g(ss) = rmid_tau_max_g_sim; %defined using the GRADIENT wind
    rmid_max_g(ss) = rmid_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
%    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
%    Vmax_max(ss) = Vmax_max_sim;
%    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
%    rmax_max(ss) = rmax_max_sim;
%    rmid_tau_max(ss) = rmid_tau_max_sim; %defined using the GRADIENT wind
%    rmid_max(ss) = rmid_max_sim;
%    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
%    r0_max(ss) = r0_max_sim;
%    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
%    r0Lil_max(ss) = r0Lil_max_sim;
    
    %% Save data for all simulations
%    Vmax_movave_all(:,ss)=Vmax_movave_sim;
%    rmax_movave_all(:,ss)=rmax_movave_sim;
%    rmid_movave_all(:,ss)=rmid_movave_sim;
%    r0_movave_all(:,ss)=r0_movave_sim;
%    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;

    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    rmid_movave_g_all(:,ss)=rmid_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
%    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;

    %% Extract and keep mpi %%%%%
    mpi_all(ss) = mpi;
end

%% PLOTTING %%%%%%%%%%%%%%%%
%Single/Multi simulation: Plot user-defined radial wind profiles

hold off
for i=1:numruns
    i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
    xvals_pl = xvals_sub_all{i}(i_xvals_pl);

    %%DRC 01/31/12
%        angmom = 1000*xvals_pl'.*data_tmean_usr_all{i}(i_xvals_pl) + .5*fcor*(1000*xvals_pl').^2;
%        dangmomdr = (angmom(3:end)-angmom(1:end-2))/(2*dx);
    %plot(xvals_pl,angmom,pl_clrs{i},'LineWidth',1)
%        plot(xvals_pl(2:end-1),dangmomdr,pl_clrs{i},'LineWidth',1)
%        angmomloss(i) = sum(dangmomdr)
    %%DRC

    r_nondim = xvals_pl/(mpi_all(i)/fcor/1000);
    V_nondim = data_tmean_usr_g_all{i}(i_xvals_pl)/mpi_all(i);
    plot(r_nondim,V_nondim,pl_clrs{i},'LineWidth',1)
    hold on
    mpi_all(i)
    
    %%identify first radial point beyond rmax with V_nondim < .1
    clear indices
    indices = find(xvals_pl>rmax_equil_g(i));
    r_div(i) = r_nondim(find(V_nondim(indices)<=.1,1)+indices(1)-1);
    
end
plot(mean(r_div),.1,'o','MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize',12)
input_title = sprintf('Time-mean radial GRADIENT wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
title(input_title)
xlabel('r / (V_p/f)')
ylabel('V_g / V_p')
legend({'0.25 K day^{-1}' '0.5 K day^{-1}' '1 K day^{-1}' '2 K day^{-1}'})
grid on


end


cd Papers/RCE_equilibrium/
