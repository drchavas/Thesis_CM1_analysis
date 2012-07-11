%Cd_scaling.m

%Created: 28 Mar 2012, Dan Chavas

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

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
   
%%Simulations
subdirs = {
%%DRAG COEFFICIENT
'CTRLv0qrhSATqdz5000_nx3072_Cddiv8'
'CTRLv0qrhSATqdz5000_nx3072_Cddiv4'
'CTRLv0qrhSATqdz5000_nx3072_Cddiv2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Cdx2'
'CTRLv0qrhSATqdz5000_nx3072_Cdx4'
'CTRLv0qrhSATqdz5000_nx3072_Cdx8'

}; %name of sub-directory with nc files
subdirs_set = subdirs;
Cds = .0015*[1/8 1/4 1/2 1 2 4 8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numruns=length(subdirs);  %total number of runs you want to plot (taken from below)

%%Append ax or 3d to subdir names for plot
for i=1:length(subdirs)
    if(run_types(i)==1) %ax
        subdirs_load{i}=sprintf('ax%s',subdirs{i});
    else    %3d
        subdirs_load{i}=sprintf('3d%s',subdirs{i});
    end
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;

for ss=1:numruns

    %%Load data for given simulation
    file_in=sprintf('simdata/%s.mat',subdirs_load{ss});
    load(file_in);

    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    r12_equil(ss) = r12_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    r12_equil_g(ss) = r12_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    
    %timescales to those values
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    r12_tau_equil(ss) = r12_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    r12_tau_equil_g(ss) = r12_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
    Vmax_gen(ss) = Vmax_gen_sim;
    rmax_gen(ss) = rmax_gen_sim;
    r12_gen(ss) = r12_gen_sim;
    r0_gen(ss) = r0_gen_sim;
    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    r12_gen_g(ss) = r12_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    r12_tau_max_g(ss) = r12_tau_max_g_sim; %defined using the GRADIENT wind
    r12_max_g(ss) = r12_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
    Vmax_max(ss) = Vmax_max_sim;
    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
    rmax_max(ss) = rmax_max_sim;
    r12_tau_max(ss) = r12_tau_max_sim; %defined using the GRADIENT wind
    r12_max(ss) = r12_max_sim;
    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
    r0_max(ss) = r0_max_sim;
    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
    r0Lil_max(ss) = r0Lil_max_sim;
    
    %% Save data for all simulations
    Vmax_movave_all(:,ss)=Vmax_movave_sim;
    rmax_movave_all(:,ss)=rmax_movave_sim;
    r12_movave_all(:,ss)=r12_movave_sim;
    r0_movave_all(:,ss)=r0_movave_sim;
    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;
    
    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    r12_movave_g_all(:,ss)=r12_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;
    
    %% Other
    fcor_all(ss) = fcor;
    mpi_all(ss) = mpi;
    lh_all(ss) = lh;

end

%% Calculate non-dim number, C
C = Cds;
i_ctrl = find(Cds==.0015,1);
C_ctrl = C(i_ctrl);
multipliers = log2(C/C_ctrl);

%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

xvals_pl = multipliers;        %values defined by user at top
[junk isort] = sort(xvals_pl);
clear junk
list=(1:length(xvals_pl))';

h=figure(1)

set(h,'Position',[360 278 375 700])

dat_max=0;
dat_min = 0;

pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 3]);

%subplot(3,1,1)
ax1=axes('position',[0.15    0.70    0.70    0.27]);
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
text1=text(-2.65,2.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
grid on
hold on
text1=text(.8,1.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_VmVp,p1(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text1,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax1,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

VmVp=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

%subplot(3,1,1)
ax2=axes('position',[0.15    0.40    0.70    0.27]);
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
text1=text(-2.65,2.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text2=text(.8,1.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_rmVpf,p2(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text2,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax2,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

rmVpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

%subplot(3,1,1)
ax3=axes('position',[0.15    0.10    0.70    0.27]);
axes(ax3)
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
text1=text(-2.65,2.4,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(C_d/C_d^*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text3=text(.8,1.5,sprintf('$\\alpha = $ %5.2f (p = %5.2f)',exp_r0Vpf,p3(2)),'fontweight','bold','FontSize',12,'Color','r')
set(text3,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','Latex');
set(ax3,'YTick',[-4 -3 -2 -1 0 1 2 3 4],'XTick',[-4 -3 -2 -1 0 1 2 3 4])
box on

r0Vpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' mpi_all(isort)' fcor_all(isort)'*10^5 lh_all(isort)'/1000];

%input_title1=sprintf('Non-dim variables = G(C), C = V_p / (f l_h)');
%input_title2=sprintf('X* = C* = %5.2f [%s]; V_p* =%5.0f [ms-1]; f* =%5.1f [*10^5 s-1]; l_h* =%5.1f [km]',C_ctrl,units,mpi_all(i_ctrl),10^5*fcor_all(i_ctrl),.001*lh_all(i_ctrl));
%suptitle({input_title1,input_title2})
%grid on

cd Papers/RCE_equilibrium/


