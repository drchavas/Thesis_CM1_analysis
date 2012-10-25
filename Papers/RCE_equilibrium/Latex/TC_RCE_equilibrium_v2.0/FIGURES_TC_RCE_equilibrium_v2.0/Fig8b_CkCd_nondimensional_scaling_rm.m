%Fig8b_CkCd_nondimensional_scaling_rm.m

%Created: 19 Sep 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%

subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'Cd_drag'};    %nondim=single-parm; nondim2=extremes
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

ax1=axes('position',[0.15    0.15    0.80    0.80]);

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

%% Calculate non-dim number, C
C = CkCd;
i_ctrl = find(CkCd==1,1);
C_ctrl = C(i_ctrl);
multipliers = log2(C/C_ctrl);

%% PLOTTING %%%%%%%%%%%%%%%%

xvals_pl(m,1:numruns) = multipliers;        %values defined by user at top
xvals_pl = xvals_pl(m,1:numruns);

if(m==length(sim_sets))
    xvals_all = [xvals_all;xvals_pl(m,1:numruns)];
end

dat_max=0;
dat_min=0;

pl_edge = max([abs(floor(min(xvals_pl(m,1:numruns)))) abs(ceil(max(xvals_pl(m,1:numruns)))) 3]);

%%%%TESTING%%%%%%
%rmax_equil_g_ts_max = 1.2*rmax_equil_g;
%rmax_equil_g_ts_min = .8*rmax_equil_g;
%%%%%%%%%%%%%%%%%%%%
    
    rmax_equil_g_ts_max_ctrl = rmax_equil_g_ts_max(i_ctrl);
    rmax_equil_g_ts_min_ctrl = rmax_equil_g_ts_min(i_ctrl);
        
%% rmax
%subplot(3,1,1)
axes(ax1)
data_temp = rmax_equil_g./mpi_all;
L_temp = rmax_equil_g_ts_min./mpi_all;
U_temp = rmax_equil_g_ts_max./mpi_all;

data_pl(m,1:numruns) = log2(data_temp./data_temp(i_ctrl));
data_pl = data_pl(m,1:numruns);
L = data_pl - log2(L_temp./data_temp(i_ctrl));
U = data_pl - log2(U_temp./data_temp(i_ctrl));

dat_max = max(dat_max,max(data_pl(m,1:numruns)));
dat_min = min(dat_min,min(data_pl(m,1:numruns)));

stats = regstats(data_pl(m,1:numruns),xvals_pl(m,1:numruns),'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p1 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
%rsq1 = stats.adjrsquare -- not useful for a simple linear fit
exp_rmVpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;

%clr1 = [0 0 1];
clr1 = [0 0 1];
clr2 = [0 1 0];
clr3 = [1 0 0];
%clr3 = [0 0 1];
%clr2 = clr3 + .4*[1 1 0];
%clr1 = clr2 + .4*[1 1 0];
switch m
    case 1,
        %% Plot error bars
        h2 = errorbar(xvals_pl,data_pl,L,U,'LineStyle','none','Color',[.5 .5 .5])
        hold on
    
        %% Plot data        
        h1(m) = plot(xvals_pl(m,1:numruns),data_pl(m,1:numruns),'x','MarkerEdgeColor',clr1,'MarkerSize',10,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color','k','LineWidth',2)
        
        text1=text(1.5,2.5,sprintf('$\\alpha_3 = $ %5.2f',exp_rmVpf),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
        text2=text(1.5,2.0,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

    case 2,
        plot(xvals_pl(m,1:numruns),data_pl(m,1:numruns),'x','MarkerEdgeColor',clr2,'MarkerSize',20,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color',clr2,'LineWidth',2)

        text1=text(2,3.1,sprintf('$\\alpha_2 = $ %5.2f',exp_rmVpf),'fontweight','bold','FontSize',24,'Color',clr2)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');    
        text2=text(1.5,2.6,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

    case 3,
        plot(xvals_pl(m,1:numruns),data_pl(m,1:numruns),'x','MarkerFaceColor',clr3,'MarkerEdgeColor',clr3,'MarkerSize',20,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color',clr3,'LineWidth',2)
        
        text1=text(2,3.7,sprintf('$\\alpha_1 = $ %5.2f',exp_rmVpf),'fontweight','bold','FontSize',24,'Color',clr3)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');    
        text2=text(1.5,3.2,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

end

input_title1=sprintf('$\\frac{r_m}{V_p/f}$');
text1=text(-2.75,-2.75,input_title1,'FontSize',30);
set(text1,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('$\\log_2(\widetilde{r}_m/\widetilde{r}_m^*)$','Interpreter','Latex','FontSize',18)
xlabel('$\\log_2(\frac{C_k}{C_d}/{\frac{C_k}{C_d}}^*)$','Interpreter','Latex','FontSize',18)
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
grid on
hold on

set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

end

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig8b_CkCd_nondimensional_scaling_rm.pdf
