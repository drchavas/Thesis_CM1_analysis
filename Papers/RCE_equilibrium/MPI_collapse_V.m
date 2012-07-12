%MPI_collapse_V.m

%Created: 14 Feb 2012, Dan Chavas

%This file is the same as TC_stats_plotdim.m, except that it plots with color-coding the 4 input variables that modulate MPI
%and displays them all togther based upon their respective mpi values.

clear all
%close all
clear
clc
figure(1)
clf(1)

cd ../..

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%clf(1)

%%variables of interest (sim_set name): 'dx' 'dz' 'domain' 'lh' 'lv' 'H' 'Qrad' 'Vpot' 'cor' 'qro' 'ro' 'rodrmax'
sim_sets = {'Tsst' 'Ttpp' 'Qcool' 'usfc'}
%sim_sets = {'Ttpp'}
T_mean = 5; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad -- DOESN"T MATTER FOR Vmax!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat_max=0;
dat_min=0;

%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

if(wrad_const == 1)
    wrad_str = 'ctrl';
else
    wrad_str = 'rce';
end

xvals_pl_all = [];
data_pl_all = [];

for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out2,sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    pl_shapes={'x' '*' 's' 'd' '+'};
    
    %%Adjust MPI for u_sfc runs only
%{
    if(strcmp('usfc',sim_set))
    VmVp = .7790;   %=Vmax_equil_g_CTRL/mpi_CTRL
    mpi_all = Vmax_equil_g/VmVp; %%u_sfc adjustment DRC 07 Jun 2012
    end
%}
    
    i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    mpi_ctrl = mpi_all(i_ctrl);
    Vmax_equil_g_ctrl = Vmax_equil_g(i_ctrl);

    [junk i_sort] = sort(mpi_all);
    clear junk
    multipliers = log2(mpi_all(i_sort)/mpi_ctrl);
    
    xvals_pl = multipliers;    %values defined by user at top
    
    figure(1)
%    subplot(3,1,1)
    if(m==1)
    ax1=axes('position',[0.15    0.35    0.70    0.56]);
    end
    axes(ax1)
    data_temp = Vmax_equil_g(i_sort);
    data_pl = log2(data_temp./Vmax_equil_g_ctrl);
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    scatter(xvals_pl,data_pl,'Marker',pl_shapes{m},'MarkerEdgeColor',pl_clrs{m})
    hold on
    
    %%need to accumulate all points into single vector for xvals and data
    xvals_pl_all = [xvals_pl_all xvals_pl];
    data_pl_all = [data_pl_all data_pl];
   
end

%% Plot a best-fit line to the data
%options = fitoptions('Method','Smooth','SmoothingParam',0.3)
%f = fit(xvals_pl_all', data_pl_all', 'smooth',options)
f = fit(xvals_pl_all', data_pl_all', 'poly1')
% plot(f, xvals_pl_all, data_pl_all)
%%Linear model: f(x) = p1*x + p2
 xmin_pl = min(xvals_pl_all);
 xmax_pl = max(xvals_pl_all);
 xdiff_pl = xmax_pl-xmin_pl;
 xfit = xmin_pl-xdiff_pl/10:xdiff_pl/20:xmax_pl+xdiff_pl/10
 yfit = f.p1.*xfit + f.p2;
 plot(xfit,yfit,'r--')

 
%% Define the size of the figure
h=figure(1)
set(h,'Position',[360 278 375 350])


%% Plot aesthetics
%subplot(3,1,1)
axes(ax1)
axis([-2 2 -2 2])
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(V_p/V^*_p)')})
input_title1=sprintf('$V_m$');
text1=text(-1.7,1.7,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
%set(text1,'EdgeColor','k')
set(ax1,'YTick',[-2 -1 0 1 2],'XTick',[-2 -1 0 1 2])
grid on
box on

h=legend({'T_{sst}' 'T_{tpp}' 'Q_{cool}' 'u_{sfc}'},'Orientation','horizontal','Position',[0.15    0.15    0.70    0.02],'EdgeColor','white')
grid on

if(equil_dynamic == 1)
    input_title = sprintf('Equilibrium: dynamic %i day; T_{mean} = %i; wrad: %s',dt_equil,T_mean,wrad_str);
else
    input_title = sprintf('Equilibrium: days %i-%i ; T_{mean} = %i; wrad: %s',tf-dt_final,tf,T_mean,wrad_str);
end
title(input_title)

cd Papers/RCE_equilibrium/

