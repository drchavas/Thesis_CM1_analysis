%TC_stats.m

%Created: 10-05-11, Dan Chavas

%Purpose: Plot sensitivity of storm structure to individual scale changes
%1) 

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc


%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

tmean0_usr = 70;    %[day]
tmeanf_usr = 100;    %[day]

T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
dt_equil = 30;  %[days]; how long must be quasi-steady to define equilibrium
Cd_in = 1.5e-3; %only used to calculate r0_Lilly

%SINGLE SIMULATION
plot_ts = 1;    %0=no; plot time series of vmax, rmax, r12, r0, r0Lil
plot_proftevol = 1; %0=no; plot radial wind profile at multiple times

%SINGLE/MULTI SIMULATION
plot_usrprof = 0;   %0=no; plot radial wind profile for user-defined time-means

%MULTI SIMULATION
plot_ts_multi = 0;    %0=no; plot time series of vmax, rmax, r12, r0
plot_stats = 0; %0=no; plot comparisons of basic statistics: Vmax, Rmax, R0, tau_equil for each; can be done at "equil time", genesis time, time of first V~Vmax, final (day 70-100 ave)
    save_file = 1;  %0=don't save; 1=save file with file name [out_file]
    i_ctrl = 4;   %index of control simulation value (i.e. place in list where 'CTRLv0qrhSATqdz5000_nx3072' is found)
    units = '-';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
   
%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 2500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
sim_set = 'nondim';  %name out output subdir (within simsets/PLOTS/[sim_set]/) where plots will be saved
Vpots = [49 93 93 93 49 93 93 49 40 40 121 121 80 80 160 160 160 71 74 74 115 80 80 40 160 160 40 84 117 94 80 120 125];
lhs = [1500 1500 3000 1500 1500 1500 750 750 1500 1500 1500 1500 1500 1500 1500 1500 1500 1500 1500 750 1500 3000 750 3000 3000 750 3000 1500 1500 1500 1500 1500 1500];
subdirs = {

%%NON-DIM NUMBER SCALING Vpot/(f*lh)
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'    %49
'CTRLv0qrhSATqdz5000_nx3072_fx2'    %93
'CTRLv0qrhSATqdz5000_nx3072_lh3000' %93
'CTRLv0qrhSATqdz5000_nx3072'    %93
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2'  %49
'CTRLv0qrhSATqdz5000_nx3072_fdiv2'  %93
'CTRLv0qrhSATqdz5000_nx3072_lh750'  %93
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh750'    %49
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2'   %40
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2' %40
'CTRLv0qrhSATqdz5000_nx3072_usfc1_fdiv2'    %121
'CTRLv0qrhSATqdz5000_nx3072_usfc1_fx2'  %121
'CTRLv0qrhSATqdz5000_nx3072_usfc5_fdiv2'    %80
'CTRLv0qrhSATqdz5000_nx3072_usfc5_fx2'  %80
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_fx2'  %160
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_fdiv2'    %160
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1'  %160
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_Tthresh250K_usfc1'   %71
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_usfc5'   %74
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_lh750' %74
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_Tthresh150K'   %115
'CTRLv0qrhSATqdz5000_nx3072_usfc5_lh3000'   %80
'CTRLv0qrhSATqdz5000_nx3072_usfc5_lh750'    %80
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_lh3000'    %40
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_lh3000'   %160
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_lh750'    %160
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2_lh3000'    %40
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K' %84
'CTRLv0qrhSATqdz5000_nx3072_rad1K' %117
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K' %94
'CTRLv0qrhSATqdz5000_nx3072_usfc5' %80 m/s
'CTRLv0qrhSATqdz5000_nx3072_usfc1' %120 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'    %125 m/s
}; %name of sub-directory with nc files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



numruns=length(subdirs);  %total number of runs you want to plot (taken from below)

%%want data for entire simulation
t0 = 0;
tf = 100;

%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=100000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=100;  %last z grid point [0,end]

%%Define subset of points for analysis
rmin_sub = 0;  %[km]; lowest value plotted
rmax_sub = 10000;    %[km]; highest value plotted
zmin_subsub = .5;  %[km]; lowest value plotted
zmax_subsub = 1; %[km]; highest value plotted

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
    
    istr = strfind(file_in,'SST');
    if(isempty(istr))
        sstk(ss) = 300.00;
    else
        sstk(ss) = str2num(file_in(istr+3:istr+8));
    end
    
    istr = strfind(file_in,'Tthresh');
    if(isempty(istr))
        Ttpp(ss) = 200;
    else
        Ttpp(ss) = str2num(file_in(istr+7:istr+9));
    end
    
    istr = strfind(file_in,'usfc');
    if(isempty(istr))
        usfc(ss) = 3;
    else
        usfc(ss) = str2num(file_in(istr+4:istr+4));
    end
    
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
    fcors(ss) = fcor;
    %Vpots(ss) = Vpot;
    %lhs(ss) = lh;

end

cd simsets/PLOTS/
mkdir(sim_set)
cd ../..


%% Calculate non-dim number, C
C = Vpots./(fcors.*lhs);
CTRL_val = C(i_ctrl);
multipliers = log2(C/CTRL_val);

%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',8,'defaultaxesfontweight','bold')

xvals_pl = multipliers;        %values defined by user at top
[junk isort] = sort(xvals_pl);
clear junk
list=(1:length(xvals_pl))';

figure(77)
clf(77)

hold off
dat_max=0;
dat_min = 0;

pl_edge = max([abs(floor(min(xvals_pl))) abs(ceil(max(xvals_pl))) 3]);

subplot(2,2,1)
data_temp = Vmax_equil_g./Vpots;
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'g*')
p = polyfit(xvals_pl,data_pl,1);
exp_VmVp = p(1);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,p(1).*xfit+p(2),'r--')
title('V_m / V_p')
axis([-pl_edge pl_edge -pl_edge pl_edge])
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
grid on
hold on
text(-1.8,-2.5,sprintf('exp=%5.2f',exp_VmVp),'fontweight','bold')

VmVp=[list list(isort) xvals_pl(isort)' data_pl(isort)' Vpots(isort)' fcors(isort)'*10^5 lhs(isort)'/1000]

subplot(2,2,2)
data_temp = rmax_equil_g./r0Lil_equil_g;
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'b+')
p = polyfit(xvals_pl,data_pl,1);
exp_rmr0 = p(1);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,p(1).*xfit+p(2),'r--')
title('r_m / r_0')
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text(-1.8,-2.5,sprintf('exp=%5.2f',exp_rmr0),'fontweight','bold')

rmr0=[list list(isort) xvals_pl(isort)' data_pl(isort)' Vpots(isort)' fcors(isort)'*10^5 lhs(isort)'/1000]

subplot(2,2,3)
data_temp = rmax_equil_g./(Vpots./fcors);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'bx')
p = polyfit(xvals_pl,data_pl,1);
exp_rmVpf = p(1);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,p(1).*xfit+p(2),'r--')
title('r_m / (V_p/ f )')
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text(-1.8,-2.5,sprintf('exp=%5.2f',exp_rmVpf),'fontweight','bold')

rmVpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' Vpots(isort)' fcors(isort)'*10^5 lhs(isort)'/1000]

subplot(2,2,4)
data_temp = r0Lil_equil_g./(Vpots./fcors);
data_pl = log2(data_temp./data_temp(i_ctrl));
dat_max = max(dat_max,max(data_pl));
dat_min = min(dat_min,min(data_pl));
scatter(xvals_pl,data_pl,'m*')
p = polyfit(xvals_pl,data_pl,1);
exp_r0Vpf = p(1);
hold on
xfit = -pl_edge:.1:pl_edge;
plot(xfit,p(1).*xfit+p(2),'r--')
title('r_0 / (V_p/ f )')
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
axis([-pl_edge pl_edge -pl_edge pl_edge])
grid on
hold on
text(-1.8,-2.5,sprintf('exp=%5.2f',exp_r0Vpf),'fontweight','bold')


r0Vpf=[list list(isort) xvals_pl(isort)' data_pl(isort)' Vpots(isort)' fcors(isort)'*10^5 lhs(isort)'/1000]

input_title1=sprintf('Non-dim variables = G(C), C = V_p / (f l_h)');
input_title2=sprintf('X* = C* = %5.2f [%s]; V_p* =%5.0f [ms-1]; f* =%5.1f [*10^5 s-1]; l_h* =%5.1f [km]',CTRL_val,units,Vpots(i_ctrl),10^5*fcors(i_ctrl),.001*lhs(i_ctrl));
suptitle({input_title1,input_title2})
grid on

%axis([-2 2 -2 2])
%ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})

%legend({'V_m/V_p','r_m/r_0','r_m/(V_p/f)','r_0/(V_p/f)'},'Location','SouthEast')
%input_title1=sprintf('GRADIENT Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_equil,tf,sim_set,units);
%input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_equil_g(i_ctrl),rmax_equil_g(i_ctrl),r0Lil_equil_g(i_ctrl));
%title({input_title1,input_title2})
%grid on

cd(sprintf('simsets/PLOTS/%s/',sim_set))
saveas(gcf,'nondim','jpeg')
cd ../../..
