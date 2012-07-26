%Qcool_diverge.m

%Created: 11 June 2012, Dan Chavas

%Purpose: This file plots the normalized equilibrium radial gradient wind profiles at
%the BL top for varying radiative cooling rate to demonstrate at what
%radius varying the radiative cooling rate becomes important.

clear
clc
close all
figure(1)
clf(1)

cd ../..

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

t0 = 0;

T_mean = 2; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
equil_dynamic = 0;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 1500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
sim_sets = {'Qcool'};  %name out output subdir (within simsets_Tmean#/PLOTS/[sim_set]/) where plots will be saved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

if(wrad_const == 1)
    wrad_str = 'ctrl';
else
    wrad_str = 'rce';
end

xvals_pl = nan(1000,2);
data_pl = nan(1000,2);

for m=1:length(sim_sets)

    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out,sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    numruns = length(subdirs_set);
    
    clear fcor_all lh_all mpi_all
    for ss=1:numruns

        subdir = subdirs_set{ss};
        %%Determine output subdirectory pathname for given sim_set
        if(equil_dynamic == 1)
            if(wrad_const == 1)
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
            else
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
            end
        else
            if(wrad_const == 1)
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
            else
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
            end
        end

        if(wrad_const == 1)
            wrad_str = 'ctrl';
        else
            wrad_str = 'rce';
        end
        %%Load data for given simulation
        load(sprintf('%s/ax%s.mat',subdir_out2,subdir));

        %% Other
        fcor_all(ss) = fcor;
        mpi_all(ss) = mpi;
        lh_all(ss) = lh;

    end
    
end

%% PLOTTING %%%%%%%%%%%%%%%%
%Single/Multi simulation: Plot user-defined radial wind profiles

hold off
for i=1:numruns-1
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
if(equil_dynamic == 0)
    input_title = sprintf('Time-mean radial GRADIENT wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
else
    input_title = sprintf('Time-mean radial GRADIENT wind profiles: dynamic; z=%5.2f %s',zvals(i_zvals),zunits);
end
title(input_title)
xlabel('r / (V_p/f)')
ylabel('V_g / V_p')
%legend({'0.25 K day^{-1}' '0.5 K day^{-1}' '1 K day^{-1}' '2 K day^{-1}' '4 K day^{-1}'})
legend({'0.25 K day^{-1}' '0.5 K day^{-1}' '1 K day^{-1}' '2 K day^{-1}'})
grid on


cd Papers/RCE_equilibrium/
