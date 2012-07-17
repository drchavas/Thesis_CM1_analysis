%Nondimensional_scaling.m

%Created: 12 Jan 2012, Dan Chavas

%This file is the same as TC_stats_plotnondim.m, except that it only produces
%three of the subplots, corresponding to the scalings for Vm/Vp, rm/(Vp/f), and r0/(Vp/f).

clear
clc
figure(1)
clf(1)

cd ../..

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'nondim'}
T_mean = 5; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dat_max=0;
dat_min=0;

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


for m=1:length(sim_sets)

    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out,sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    numruns = length(subdirs_set);
    
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

%% Calculate non-dim number, C
C = mpi_all./(fcor_all.*lh_all);
i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
C_ctrl = C(i_ctrl);
multipliers = log2(C/C_ctrl);

%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

xvals_pl = multipliers;        %values defined by user at top
[junk isort] = sort(xvals_pl);
clear junk
list=(1:length(xvals_pl))';

h=figure(1)

set(h,'Position',[360 278 400 700])

dat_max=0;
dat_min = 0;

pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 3]);

%subplot(3,1,1)
%ax1=axes('position',[0.15    0.70    0.70    0.27]);
ax1=axes('position',[0.1    0.6    0.35    0.35]);
axes(ax1)
data_temp = Vmax_equil_g./mpi_all;
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'g*')
stats = regstats(data_pl,xvals_pl,'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p1 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
%rsq1 = stats.adjrsquare -- not useful for a simple linear fit
exp_VmVp = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,coefs(1)+coefs(2).*xfit,'r--')
input_title1=sprintf('$\\frac{V_m}{V_p}$');
text1=text(-4.65,4.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
grid on
hold on
text1=text(1.8,2.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_VmVp,p1(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

VmVp=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

%subplot(3,1,1)
%ax2=axes('position',[0.15    0.40    0.70    0.27]);
ax2=axes('position',[0.55    0.6    0.35    0.35]);
axes(ax2)
data_temp = rmax_equil_g./(mpi_all./fcor_all);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'bx')
stats = regstats(data_pl,xvals_pl,'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p2 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
exp_rmVpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,coefs(1)+coefs(2).*xfit,'r--')
input_title1=sprintf('$\\frac{r_m}{V_p/f}$');
text1=text(-4.65,4.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text2=text(1.8,2.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_rmVpf,p2(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax2,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

rmVpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];


%subplot(3,1,1)
%ax3=axes('position',[0.15    0.10    0.70    0.27]);
ax3=axes('position',[0.1    0.2    0.35    0.35]);
axes(ax3)
data_temp = r0ER11_equil_g./(mpi_all./fcor_all);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'m*')
stats = regstats(data_pl,xvals_pl,'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p3 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
exp_r0Vpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,coefs(1)+coefs(2).*xfit,'r--')
input_title1=sprintf('$\\frac{r_{0 ER11}}{V_p/f}$');
text1=text(-4.65,4.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(C/C*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text3=text(1.8,2.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_r0Vpf,p3(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text3,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax3,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

r0ER11Vpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];


%subplot(3,1,1)
%ax4=axes('position',[0.15    0.10    0.70    0.27]);
ax4=axes('position',[0.55    0.2    0.35    0.35]);
axes(ax4)
data_temp = r0Lil_equil_g./(mpi_all./fcor_all);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'m*')
stats = regstats(data_pl,xvals_pl,'linear'); %fit to y=beta1+beta2*x
coefs = stats.beta %coefs(1) = beta1, coefs(2) = beta2
p3 = stats.tstat.pval   %two-sided t-test: p(1) = p-value of beta1; p(2) = p-value of beta2; p<.025 = value is significantly different from zero at 95% CI (i.e. can reject null hypothesis of value = 0)
exp_r0Vpf = coefs(2);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,coefs(1)+coefs(2).*xfit,'r--')
input_title1=sprintf('$\\frac{r_0}{V_p/f}$');
text1=text(-4.65,4.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(C/C*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text3=text(1.8,2.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_r0Vpf,p3(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text3,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax4,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

r0Vpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

%input_title1=sprintf('Non-dim variables = G(C), C = V_p / (f l_h)');
%input_title2=sprintf('X* = C* = %5.2f [%s]; V_p* =%5.0f [ms-1]; f* =%5.1f [*10^5 s-1]; l_h* =%5.1f [km]',C_ctrl,units,mpi_all(i_ctrl),10^5*fcor_all(i_ctrl),.001*lh_all(i_ctrl));
%suptitle({input_title1,input_title2})
%grid on

cd Papers/RCE_equilibrium/

