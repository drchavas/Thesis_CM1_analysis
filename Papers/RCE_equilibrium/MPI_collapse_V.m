%MPI_collapse_V_poster.m

%Created: 24 Jul 2012, Dan Chavas

%Purpose: Create poster-ready MPI_collapse plot for Vmax

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'Tsst' 'usfc' 'Qcool' 'Ttpp'}
sim_sets_str = {'T_{sst}' 'u_{sfc}' 'Q_{cool}' 'T_{tpp}'};  %make sure this matches!

T_mean = 2; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',36,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
h=figure(1)
set(h,'Position',[0 0 575 575])
ax1=axes('position',[0.20    0.25    0.70    0.70]);

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
dat_max=0;
dat_min=0;

for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out2,sim_set));
    pl_clrs={'b' 'c' 'g' 'r' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    pl_shapes={'o' 's' 'd' 'v'};
    
    %%Adjust MPI for u_sfc runs only
%{
    if(strcmp('usfc',sim_set))
    VmVp = .7790;   %=Vmax_equil_g_CTRL/mpi_CTRL
    mpi_all = Vmax_equil_g/VmVp; %%u_sfc adjustment DRC 07 Jun 2012
    end
%}
    
    i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    if(strcmp(sim_set,'usfc_drag')==1)
        i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072_drag')==1,1);
    end
    mpi_ctrl = mpi_all(i_ctrl);
    Vmax_equil_g_ctrl = Vmax_equil_g(i_ctrl);

    [junk i_sort] = sort(mpi_all);
    clear junk
    multipliers = log2(mpi_all(i_sort)/mpi_ctrl);
    
    xvals_pl = multipliers;    %values defined by user at top
    
    figure(1)
%    subplot(3,1,1)
    axes(ax1)
    data_temp = Vmax_equil_g(i_sort);
    data_pl = log2(data_temp./Vmax_equil_g_ctrl);
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,pl_shapes{m},'MarkerFaceColor',pl_clrs{m},'MarkerEdgeColor','k','MarkerSize',20)
    hold on
    
    %%need to accumulate all points into single vector for xvals and data
    xvals_pl_all = [xvals_pl_all xvals_pl];
    data_pl_all = [data_pl_all data_pl];
   
end
    
pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 2]);

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
 plot(xfit,yfit,'--','Color',[.5 0 0],'LineWidth',3)


input_title1=sprintf('$V_m$');
text1=text(-1.9,1.9,input_title1,'FontSize',60);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('$\\log_2(V_m/V_m^*)$','Interpreter','Latex')
xlabel('$\\log_2(V_p/V_p^*)$','Interpreter','Latex')
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
grid on
set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on
h=legend(sim_sets_str,'Orientation','vertical','Position',[0.85    0.5    0.02    0.20],'EdgeColor','white')
grid on

cd Papers/RCE_equilibrium/

