%time_series_compare.m

%Created: 28 Mar 2012, Dan Chavas

%This file is the same as TC_stats_plot.m, but prettier

clear
clc

cd ../..

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

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
subdirs = {

%{
'CTRLv12.5qrhSATqroVpfdiv128qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv64qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv32qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv16qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv8qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv4qdz5000_nx3072'
%'CTRLv12.5qrhSATqroVpfdiv2qdz5000_nx3072'
%'CTRLv12.5qrhSATqroVpfqdz5000_nx3072'
%}

'CTRLv0qrhSATqroVpfdiv128qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv64qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv32qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv16qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv8qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv4qdz5000_nx3072'
%'CTRLv0qrhSATqroVpfdiv2qdz5000_nx3072'
%'CTRLv0qrhSATqroVpfqdz5000_nx3072'
%}
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
    load(sprintf('simdata/%s.mat',subdirs_load{ss}));
    
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


end


%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
%Single simulation: Plot time-series of Vmax, rmax, r12, r0 for both V and Vg

%Plot time series
h=figure(1)
clf(1)

set(h,'Position',[160 578 575 400])

ax1=axes('position',[0.15    0.15    0.75    0.75]);
axes(ax1)

for i=1:numruns

    subplot(3,1,1)
    Vm_plot = Vmax_movave_g_all(:,i)/Vmax_equil_g(i);
    plot(t_day,Vm_plot,pl_clrs{i})
    hold on
%    if(~isnan(Vmax_tau_equil_g(i)))
%        h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i)/max(Vmax_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
%    end


    subplot(3,1,2)
    rm_plot = rmax_movave_g_all(:,i)/rmax_equil_g(i);
    plot(t_day(Vm_plot>.7),rm_plot(Vm_plot>.7),pl_clrs{i})
    hold on
%    if(~isnan(rmax_tau_equil_g(i)))
%        h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i)/max(rmax_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
%    end


    subplot(3,1,3)
    r0_plot = r0Lil_movave_g_all(:,i)/r0Lil_equil_g(i);
    plot(t_day(Vm_plot>.7),r0_plot(Vm_plot>.7),pl_clrs{i})
    hold on
%    if(~isnan(r0Lil_tau_equil_g(i)))
%        h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i)/max(r0Lil_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',.5*[.49 1 .63]);
%    end

%    plot(Vmax_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{1},'MarkerSize',12)
%    plot(rmax_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{2},'MarkerSize',12)
%    plot(r0Lil_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{3},'MarkerSize',12)
end

input_legend = {'1/128','1/64','1/32','1/16','1/8','1/4','1/2','1'};
h_legend=legend(input_legend,'Location','SouthEast')
set(h_legend,'FontSize',9);


subplot(3,1,1)
plot(t_day,.9*ones(length(t_day),1),'k--','LineWidth',1)
plot(t_day,1.1*ones(length(t_day),1),'k--','LineWidth',1)
subplot(3,1,2)
plot(t_day,.9*ones(length(t_day),1),'k--','LineWidth',1)
plot(t_day,1.1*ones(length(t_day),1),'k--','LineWidth',1)
subplot(3,1,3)
plot(t_day,.9*ones(length(t_day),1),'k--','LineWidth',1)
plot(t_day,1.1*ones(length(t_day),1),'k--','LineWidth',1)


subplot(3,1,1)
input_title1=sprintf('Evolution of V_m/V_m^*');
%input_title2=sprintf('Max values: V_m=%3.0f ms^{-1}; r_m=%3.0f km; r_0=%4.0f km',max(Vmax_movave_g_all(:,i)),max(rmax_movave_g_all(:,i)),max(r0Lil_movave_g_all(:,i)));
%title({input_title1,input_title2})
%title(input_title1)
%ylabel('Fraction of equilibrium value')
ymax = 1.5;%1*max(Vm_plot)
axis([t0 tf/3 .2 ymax])
grid on

subplot(3,1,2)
input_title1=sprintf('Evolution of r_m/r_m^*');
%input_title2=sprintf('Max values: V_m=%3.0f ms^{-1}; r_m=%3.0f km; r_0=%4.0f km',max(Vmax_movave_g_all(:,i)),max(rmax_movave_g_all(:,i)),max(r0Lil_movave_g_all(:,i)));
%title({input_title1,input_title2})
%title(input_title1)
%ylabel('Fraction of equilibrium value')
ymax = 2;%1*max(rm_plot(Vm_plot>.7))
axis([t0 tf/3 .2 ymax])
grid on

subplot(3,1,3)
input_title1=sprintf('Evolution of r_0/r_0^*');
%input_title2=sprintf('Max values: V_m=%3.0f ms^{-1}; r_m=%3.0f km; r_0=%4.0f km',max(Vmax_movave_g_all(:,i)),max(rmax_movave_g_all(:,i)),max(r0Lil_movave_g_all(:,i)));
%title({input_title1,input_title2})
%title(input_title1)
%ylabel('Fraction of equilibrium value')
ymax = 2;%1*max(r0_plot(Vm_plot>.7))
axis([t0 tf/3 .2 ymax])
grid on
xlabel('Time [day]')
%set(ax1,'YTick',[0 .2 .4 .6 .8 1 1.2 1.4 1.6 1.8 2],'XTick',[0 20 40 60 80 100])

cd Papers/RCE_transient/

