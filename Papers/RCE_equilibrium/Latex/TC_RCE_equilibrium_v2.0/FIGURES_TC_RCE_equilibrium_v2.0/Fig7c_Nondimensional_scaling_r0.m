%Fig7c_Nondimensional_scaling_r0.m

%Created: 6 Sep 2012, Dan Chavas

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

xvals_pl = multipliers;        %values defined by user at top
[junk isort] = sort(xvals_pl);
clear junk
list=(1:length(xvals_pl))';

if(m==length(sim_sets))
    xvals_all = [xvals_all;xvals_pl];
end

dat_max=0;
dat_min=0;

pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 5]);

%% RMAX
%subplot(3,1,1)
axes(ax1)
data_temp = r0Lil_equil_g./(mpi_all./fcor_all);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));


%% ERROR BAR DATA
%%%%TESTING%%%%%%
%r0Lil_equil_g_ts_max = 1.2*r0Lil_equil_g;
%r0Lil_equil_g_ts_min = .8*r0Lil_equil_g;
%%%%%%%%%%%%%%%%%%%%
    
data_temp_max = r0Lil_equil_g_ts_max./(mpi_all./fcor_all);
data_pl_max = log2(data_temp_max./data_temp(i_ctrl));

data_temp_min = r0Lil_equil_g_ts_min./(mpi_all./fcor_all);
data_pl_min = log2(data_temp_min./data_temp(i_ctrl));

L = data_pl - data_pl_min;
U = data_pl_max - data_pl;

rm_all = [rm_all;data_pl];


stats = regstats(data_pl,xvals_pl,'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p1 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
%rsq1 = stats.adjrsquare -- not useful for a simple linear fit
exp_r0Vpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;

clr1 = [0 .75 0];


%% Plot error bars
errorbar(xvals_pl,data_pl,L,U,'LineStyle','none','Color',[.5 .5 .5])
hold on

%% Plot values + best fit line
plot(xvals_pl,data_pl,'x','MarkerEdgeColor',clr1,'MarkerSize',10,'LineWidth',2)
hold on

%%%%% TEST DRC: Vp effect %%
%plot(xvals_pl(mpi_all>92),data_pl(mpi_all>92),'x','MarkerEdgeColor','r','MarkerSize',10,'LineWidth',2)
%hold on
%plot(xvals_pl(mpi_all<90),data_pl(mpi_all<90),'x','MarkerEdgeColor','g','MarkerSize',10,'LineWidth',2)
%title('red = high MPI, green = low MPI, blue = ctrl mpi')
%%%%% TEST DRC %%

plot(xfit,coefs(1)+coefs(2).*xfit,'--','Color','k','LineWidth',2)

%% Print best fit slope and p-value
text1=text(1.5,3.5,sprintf('$\\alpha_3 = $ %5.2f',exp_r0Vpf),'fontweight','bold','FontSize',24,'Color',clr1)
set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
text2=text(1.5,2.5,sprintf('(p = %5.2f)',p1(2)),'fontweight','bold','FontSize',24,'Color',clr1)
set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');


input_title1=sprintf('$\\frac{r_0}{V_p/f}$');
text1=text(-4.75,-4.75,input_title1,'FontSize',30);
set(text1,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('$\\log_2(\widetilde{r}_0/\widetilde{r}_0^*)$','Interpreter','Latex','FontSize',18)
%xlabel('$\\log_2(C/C^*)$','Interpreter','Latex','FontSize',18)
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
grid on
hold on

set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

VmVp=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

end

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig7c_Nondimensional_scaling_r0.pdf
