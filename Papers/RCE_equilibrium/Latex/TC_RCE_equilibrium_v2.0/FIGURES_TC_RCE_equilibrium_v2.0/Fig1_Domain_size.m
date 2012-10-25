%Fig1_Domain_size.m

%Created: 10 Jan 2012, Dan Chavas

%This file is the same as TC_stats_plot.m, except that it only produces one
%plot with the radial wind profiles of 5 simulations with different domain sizes all on a single plot (i.e. a modified version
%of figure 5)

clear
clc
close all
figure(1)
clf(1)

cd '/Users/drchavas/Documents/Research/Thesis/CM1/v15/'

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
sim_sets = {'Lx_drag'};  %name out output subdir (within simsets_Tmean#/PLOTS/[sim_set]/) where plots will be saved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out = sprintf('CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out = sprintf('CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out = sprintf('CM1_postproc_data/simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out = sprintf('CM1_postproc_data/simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
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
                subdir_out2 = sprintf('CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
            else
                subdir_out2 = sprintf('CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
            end
        else
            if(wrad_const == 1)
                subdir_out2 = sprintf('CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
            else
                subdir_out2 = sprintf('CM1_postproc_data/simdata_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
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
for i=1:numruns
    i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
    xvals_pl = xvals_sub_all{i}(i_xvals_pl);

    plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl),pl_clrs{i},'LineWidth',1)
    hold on
     
end


%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%Single/Multi simulation: Plot user-defined radial wind profiles

h=figure(1)
clf(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 20 16];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.15    0.15    0.8    0.8]);
axes(ax1)

hold off
for i=1:numruns
    i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
    xvals_pl = xvals_sub_all{i}(i_xvals_pl);

    plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl),pl_clrs{i})
    hold on
end
plot(0*(1:xvals_pl(end)),'k')
%input_title = sprintf('Equilibrium radial profiles of gradient wind for variable domain width');
%title(input_title)
xlabel('radius [km]','FontSize',18)
ylabel('azimuthal wind speed [m/s]','FontSize',18)
input_legend={'768 km','1536 km','3072 km','6144 km','12288 km'};
legend(input_legend)
set(ax1,'YTick',[-10 0 10 20 30 40 50 60 70 80],'XTick',[0 500 1000 1500])
grid on

cd 'Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0'

%save PDF
print -dpdf -r300 Fig1_Domain_size.pdf
