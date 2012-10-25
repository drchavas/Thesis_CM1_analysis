%Fig11_Qcool_Mcnst.m

%Created: 23 Oct 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%

subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'QcoolVpcnst_usfc_drag'};    %nondim=single-parm; nondim2=extremes
T_mean = 2; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

Ck = .0015;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
h=figure(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 15 15];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.15    0.15    0.80    0.70]);

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

xvals_pl = nan(2,1000);
data_pl = nan(2,1000);
xvals_all = [];
Vm_all = [];
rm_all = [];
r0_all = [];

sim_set=sim_sets{1};
load(sprintf('%s/%s',subdir_out,sim_set))

for m=1:length(sim_sets)
%{
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
%}    
%%TESTING%%%%
%mpi_all = [261.1 200.5 133.1 mpi_all(4) 63.6 44.7 31.6];
%Cd_all = .0015*[0.125 .25 .501 1 2 4 8];   %THIS SHOULD BE INCLUDED AUTOMATICALLY
%%%%%%%%%%%%%
CkCd = Ck./Cd_all;



%% PLOTTING %%%%%%%%%%%%%%%%
xvals_pl(1:numruns,m) = multipliers;        %values defined by user at top
xvals_pl = xvals_pl(1:numruns,m);

if(m==length(sim_sets))
    xvals_all = [xvals_all;xvals_pl(1:numruns,m)];
end

dat_max=0;
dat_min=0;

%%%%TESTING%%%%%%
%Vmax_equil_g_ts_max = 1.2*Vmax_equil_g;
%Vmax_equil_g_ts_min = .8*Vmax_equil_g;
%%%%%%%%%%%%%%%%%%%%
    
%    Vmax_equil_g_ts_max_ctrl = Vmax_equil_g_ts_max(i_ctrl);
%    Vmax_equil_g_ts_min_ctrl = Vmax_equil_g_ts_min(i_ctrl);

    data_pl = Vmax_equil_g.*rmax_equil_g;
%    L = data_pl - Vm_adj.*Vmax_equil_g_ts_min./mpi_all;
%    U = Vm_adj.*Vmax_equil_g_ts_max./mpi_all - data_pl;

%% Create grid for angular momentum contours over relevant domain of V and r
rgrid_min = .9*min(rmax_equil_g);
rgrid_max = 1.1*max(rmax_equil_g);
rgrid = rgrid_min:(rgrid_max-rgrid_min)/1000:rgrid_max;

Vgrid_min = .9*min(Vmax_equil_g);
Vgrid_max = 1.1*max(Vmax_equil_g);
Vgrid = Vgrid_min:(Vgrid_max-Vgrid_min)/1000:Vgrid_max;

Mrel = Vgrid'*rgrid; %Matrix of angular momentum values

%% Plot M contour grid
contour(rgrid,Vgrid,Mrel,20)
colorbar
hold on

axes(ax1)

clr1 = [0 0 1];

%% Plot error bars
%h1 = errorbar(x,data_pl,L,U,'LineStyle','none','Color',[.5 .5 .5])
%hold on

%h2 = plot(multipliers,data_pl,'x','MarkerEdgeColor',clr1,'MarkerSize',10,'LineWidth',2)
for i=1:6
    if(multipliers(i)==0)
        h2 = plot(rmax_equil_g(i),Vmax_equil_g(i),'.','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',10*2^(.5*i),'LineWidth',2)
    else
        h2 = plot(rmax_equil_g(i),Vmax_equil_g(i),'.','MarkerEdgeColor','k','MarkerSize',10*2^(.5*i),'LineWidth',2)
    end
    hold on
    
end
%xlabel('$\\log_2(Q_{cool}/Q*_{cool})$','Interpreter','Latex','FontSize',18)
xlabel('$r_m$ [$km$]','Interpreter','Latex','FontSize',18)
ylabel(sprintf('$V_m$ [$ms^{-1}$]'),'Interpreter','Latex','FontSize',18)
axis([rgrid_min rgrid_max Vgrid_min Vgrid_max])
title('M(r_m): Larger marker = 2x Qcool; gray = 1 K/day')
grid on
hold on

%set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on        
    

end

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig11_Qcool_Mcnst.pdf