%Fig7b_Nondimensional_scaling_r0.m

%Created: 7 Sep 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%
sim_sets = {'nondim_all_drag'};    %nondim=single-parm; nondim2=extremes
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

xvals_pl = nan(1000,2);
data_pl = nan(1000,2);
xvals_all = [];
Vm_all = [];
rm_all = [];
r0_all = [];

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
    
%% Calculate non-dim number, C
C = mpi_all./(fcor_all.*lh_all);
i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072_drag')==1,1);
C_ctrl = C(i_ctrl);
multipliers = log2(C/C_ctrl);

%% PLOTTING %%%%%%%%%%%%%%%%

xvals_pl(1:numruns,m) = multipliers;        %values defined by user at top
[junk isort] = sort(xvals_pl(1:numruns,m));
clear junk
list=(1:length(xvals_pl(1:numruns,m)))';

if(m==length(sim_sets))
    xvals_all = [xvals_all;xvals_pl(1:numruns,m)];
end

dat_max=0;
dat_min=0;

pl_edge = max([abs(floor(min(xvals_pl(1:numruns,m)))) abs(ceil(max(xvals_pl(1:numruns,m)))) 5]);


%% R0
%subplot(3,1,1)
axes(ax1)
data_temp = r0Lil_equil_g./(mpi_all./fcor_all);
data_pl(1:numruns,m) = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl(1:numruns,m)));
dat_min = min(dat_min,min(data_pl(1:numruns,m)));

r0_all = [r0_all;data_pl(1:numruns,m)];

stats = regstats(data_pl(1:numruns,m),xvals_pl(1:numruns,m),'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p1 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
%rsq1 = stats.adjrsquare -- not useful for a simple linear fit
exp_r0Vpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;

%clr1 = [0 0 1];
clr1 = [0 .75 0];
clr2 = [0 1 0];
clr3 = [1 0 0];
%clr3 = [0 0 1];
%clr2 = clr3 + .4*[1 1 0];
%clr1 = clr2 + .4*[1 1 0];
switch m
    case 1,
        plot(xvals_pl(1:numruns,m),data_pl(1:numruns,m),'x','MarkerEdgeColor',clr1,'MarkerSize',20,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color',clr1,'LineWidth',2)
        
        text1=text(1.5,2.5,sprintf('$\\alpha_3 = $ %5.2f',exp_r0Vpf),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
        text2=text(1.5,2.0,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

    case 2,
        plot(xvals_pl(1:numruns,m),data_pl(1:numruns,m),'x','MarkerEdgeColor',clr2,'MarkerSize',20,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color',clr2,'LineWidth',2)

        text1=text(1.5,3.1,sprintf('$\\alpha_2 = $ %5.2f',exp_r0Vpf),'fontweight','bold','FontSize',24,'Color',clr2)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');    
        text2=text(1.5,2.6,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

    case 3,
        plot(xvals_pl(1:numruns,m),data_pl(1:numruns,m),'x','MarkerFaceColor',clr3,'MarkerEdgeColor',clr3,'MarkerSize',20,'LineWidth',2)
        hold on
        plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color',clr3,'LineWidth',2)
        
        text1=text(1.5,3.7,sprintf('$\\alpha_1 = $ %5.2f',exp_r0Vpf),'fontweight','bold','FontSize',24,'Color',clr3)
        set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');    
        text2=text(1.5,3.2,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
        set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');

end

input_title1=sprintf('$\\frac{r_0}{V_p/f}$');
text1=text(-4.75,-4.75,input_title1,'FontSize',30);
set(text1,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('$\\log_2(\widetilde{r}_0/\widetilde{r}_0^*)$','Interpreter','Latex','FontSize',18)
xlabel('$\\log_2(C/C^*)$','Interpreter','Latex','FontSize',18)
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
grid on
hold on

set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

VmVpf=[list list(isort) xvals_pl(isort,m) data_pl(isort,m) mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

end

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig7c_Nondimensional_scaling_r0.pdf
