%Fig5half_LR_nocollapse_r0ctrlwcool.m

%Created: 11 Sep 2012, Dan Chavas

%This file is the same as MPI_collapse_r0ctrlwcool.m, except sorted by L_R = (N_v*H)/f
%instead of V_p.  The purpose is to show how the scaling does not work well


clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%

%%NOTE!! L_R_all IS SET TO RANDOM VALUES BELOW FOR TESTING!

subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'Tsst' 'usfc_drag' 'Qcool' 'Ttpp'}
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

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
h=figure(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 20 20];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.1    0.2    0.85    0.75]);

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
   
    i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    if(strcmp(sim_set,'usfc_drag')==1)
        i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072_drag')==1,1);
    end
    
    %%Rossby deformation radius  L_R = (N_v*H)/f

%%%%% RANDOM VALUES!!! CHANGE ME! %%%%%%%%%%%%
    L_R_all = 3000*1000*rand(1,length(multipliers))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    L_R_ctrl = L_R_all(i_ctrl);
    r0Lil_Lilctrl_equil_g_ctrl = r0Lil_Lilctrl_equil_g(i_ctrl);

    [junk i_sort] = sort(L_R_all);
    clear junk
    multipliers = log2(L_R_all(i_sort)/L_R_ctrl);
    
    xvals_pl = multipliers;    %values defined by user at top
    
    figure(1)
%    subplot(3,1,1)
    axes(ax1)
    data_temp = r0Lil_Lilctrl_equil_g(i_sort);
    data_pl = log2(data_temp./r0Lil_Lilctrl_equil_g_ctrl);
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,pl_shapes{m},'MarkerFaceColor',pl_clrs{m},'MarkerEdgeColor','k','MarkerSize',10)
    hold on
    
    %%need to accumulate all points into single vector for xvals and data
    xvals_pl_all = [xvals_pl_all xvals_pl];
    data_pl_all = [data_pl_all data_pl];
   
end
    
pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 1.5]);

%% Plot a best-fit line to the data
%{
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
%}

input_title1=sprintf('$r_{0''}$');
text1=text(-pl_edge+.1,pl_edge-.1,input_title1,'FontSize',30);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('$\\log_2(r_{0''}/r_{0''}^*)$','Interpreter','Latex','FontSize',18)
xlabel('$\\log_2(L_R/L_R^*)$','Interpreter','Latex','FontSize',18)
%xlabel('','Interpreter','Latex','FontSize',18)
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
grid on
set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on
%h=legend(sim_sets_str,'Orientation','horizontal','Position',[1 0.1  0.2    0.02],'EdgeColor','black','FontSize',18)
h=legend(sim_sets_str,'Orientation','horizontal','Location','SouthOutside','EdgeColor','black','FontSize',18)
set(h,'Position',get(h,'Position')+[0 -.18 0 0])
grid on

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig5half_LR_nocollapse_r0ctrlwcool.pdf
