%Fig2_Control_run.m

%Created: 22 May 2012, Dan Chavas

%This file plots and saves the timeseries of Vm, rm, and r0 for each
%simulation

clear
clc
close all

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_type=1; %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

T_mean = 2; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pl_clrs={'b','r',[0 .75 0],'m'};
   
%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 2500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
subdir = 'CTRLv0qrhSATqdz5000_nx3072'; %name of sub-directory with nc files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

if(wrad_const == 1)
    wrad_str = 'ctrl';
else
    wrad_str = 'rce';
end

%%Append ax or 3d to subdir names for plot
if(run_type==1) %ax
    subdir_load=sprintf('ax%s',subdir);
else    %3d
    subdir_load=sprintf('3d%s',subdir);
end

%% CALCULATE STEADY STATE RADIAL PROFILE

%%Load data for given simulation
if(equil_dynamic==1)
    load(sprintf('%s/%s.mat',subdir_out,subdir_load));
else
    load(sprintf('%s/%s.mat',subdir_out,subdir_load));
end

%% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
%variable values
Vmax_equil = Vmax_equil_sim;
rmax_equil = rmax_equil_sim;
rrad_equil = rrad_equil_sim;
r0_equil = r0_equil_sim;
r0Lil_equil = r0Lil_equil_sim;
Vmax_equil_g = Vmax_equil_g_sim;
rmax_equil_g = rmax_equil_g_sim;
rrad_equil_g = rrad_equil_g_sim;
r0_equil_g = r0_equil_g_sim;
r0Lil_equil_g = r0Lil_equil_g_sim;
r0Lil_Lilctrl_equil_g = r0Lil_Lilctrl_equil_g_sim;
r0ER11_equil_g = r0ER11_equil_g_sim;

%timescales to those values
Vmax_tau_equil = Vmax_tau_equil_sim;
rmax_tau_equil = rmax_tau_equil_sim;
rrad_tau_equil = rrad_tau_equil_sim;
r0_tau_equil = r0_tau_equil_sim;
r0Lil_tau_equil = r0Lil_tau_equil_sim;
Vmax_tau_equil_g = Vmax_tau_equil_g_sim;
rmax_tau_equil_g = rmax_tau_equil_g_sim;
rrad_tau_equil_g = rrad_tau_equil_g_sim;
r0_tau_equil_g = r0_tau_equil_g_sim;
r0Lil_tau_equil_g = r0Lil_tau_equil_g_sim;
r0Lil_Lilctrl_tau_equil_g = r0Lil_Lilctrl_tau_equil_g_sim;
r0ER11_tau_equil_g = r0ER11_tau_equil_g_sim;

%% Transient data %%%%%
tau_gen = tau_gen_sim; %defined using the GRADIENT wind
%    Vmax_gen = Vmax_gen_sim;
%    rmax_gen = rmax_gen_sim;
%    rrad_gen = rrad_gen_sim;
%    r0_gen = r0_gen_sim;
%    r0Lil_gen = r0Lil_gen_sim;
Vmax_gen_g = Vmax_gen_g_sim;
rmax_gen_g = rmax_gen_g_sim;
rrad_gen_g = rrad_gen_g_sim;
r0_gen_g = r0_gen_g_sim;
r0Lil_gen_g = r0Lil_gen_g_sim;
r0Lil_Lilctrl_gen_g = r0Lil_Lilctrl_gen_g_sim;
r0ER11_gen_g = r0ER11_gen_g_sim;

Vmax_tau_max_g = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
Vmax_max_g = Vmax_max_g_sim;
rmax_tau_max_g = rmax_tau_max_g_sim; %defined using the GRADIENT wind
rmax_max_g = rmax_max_g_sim;
rrad_tau_max_g = rrad_tau_max_g_sim; %defined using the GRADIENT wind
rrad_max_g = rrad_max_g_sim;
r0_tau_max_g = r0_tau_max_g_sim; %defined using the GRADIENT wind
r0_max_g = r0_max_g_sim;
r0Lil_tau_max_g = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
r0Lil_max_g = r0Lil_max_g_sim;
r0Lil_Lilctrl_tau_max_g = r0Lil_Lilctrl_tau_max_g_sim; %defined using the GRADIENT wind
r0Lil_Lilctrl_max_g = r0Lil_Lilctrl_max_g_sim;
r0ER11_tau_max_g = r0ER11_tau_max_g_sim; %defined using the GRADIENT wind
r0ER11_max_g = r0ER11_max_g_sim;

%    Vmax_tau_max = Vmax_tau_max_sim; %defined using the GRADIENT wind
%    Vmax_max = Vmax_max_sim;
%    rmax_tau_max = rmax_tau_max_sim; %defined using the GRADIENT wind
%    rmax_max = rmax_max_sim;
%    rrad_tau_max = rrad_tau_max_sim; %defined using the GRADIENT wind
%    rrad_max = rrad_max_sim;
%    r0_tau_max = r0_tau_max_sim; %defined using the GRADIENT wind
%    r0_max = r0_max_sim;
%    r0Lil_tau_max = r0Lil_tau_max_sim; %defined using the GRADIENT wind
%    r0Lil_max = r0Lil_max_sim;

%% Save data for all simulations
%    Vmax_movave_all=Vmax_movave_sim;
%    rmax_movave_all=rmax_movave_sim;
%    rrad_movave_all=rrad_movave_sim;
%    r0_movave_all=r0_movave_sim;
%    r0Lil_movave_all=r0Lil_movave_sim;

Vmax_movave_g_all=Vmax_movave_g_sim;
rmax_movave_g_all=rmax_movave_g_sim;
rrad_movave_g_all=rrad_movave_g_sim;
r0_movave_g_all=r0_movave_g_sim;
r0Lil_movave_g_all=r0Lil_movave_g_sim;
r0Lil_Lilctrl_movave_g_all=r0Lil_Lilctrl_movave_g_sim;
r0ER11_movave_g_all=r0ER11_movave_g_sim;

%% User profile %%%%%%%%%%%%%%%%%%%%%%%%
xvals_sub_all = xvals_sub_sim;
%    data_tmean_usr_all = data_tmean_usr_sim;
data_tmean_usr_g_all = data_tmean_usr_g_sim;

%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
%Single simulation: Plot time-series of Vmax, rmax, rrad, r0 for both V and Vg

%Plot time series
h = figure(1);
clf(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 20 16];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));


ax1=axes('position',[0.15    0.15    0.75    0.75]);

Vm_plot = Vmax_movave_g_all/Vmax_equil_g;
plot(ax1,t_day,Vm_plot,'Color',pl_clrs{1})
hold on
%    if(~isnan(Vmax_tau_equil_g(i)))
%        h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i)/max(Vmax_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
%    end

rm_plot = rmax_movave_g_all/rmax_equil_g;
plot(ax1,t_day(Vm_plot>.5),rm_plot(Vm_plot>.5),'Color',pl_clrs{2})
%    if(~isnan(rmax_tau_equil_g(i)))
%        h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i)/max(rmax_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
%    end

%r0Lil_Lilctrl_plot = r0Lil_Lilctrl_movave_g_all/r0Lil_Lilctrl_equil_g;
%plot(ax1,t_day(Vm_plot>.5),r0Lil_Lilctrl_plot(Vm_plot>.5),'Color',pl_clrs{4})
%    if(~isnan(r0Lil_tau_equil_g(i)))
%        h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i)/max(r0Lil_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',.5*[.49 1 .63]);
%    end

r0_plot = r0Lil_movave_g_all/r0Lil_equil_g;
plot(ax1,t_day(Vm_plot>.5),r0_plot(Vm_plot>.5),'Color',pl_clrs{3})
%    if(~isnan(r0Lil_tau_equil_g(i)))
%        h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i)/max(r0Lil_movave_g_all(:,i)),'d');
%        set(h,'markersize',10,'MarkerFaceColor',.5*[.49 1 .63]);
%    end

plot(ax1,t_day,.9*ones(length(t_day),1),'k--','LineWidth',1)
plot(ax1,t_day,1.1*ones(length(t_day),1),'k--','LineWidth',1)
plot(ax1,Vmax_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{1},'MarkerSize',12)
plot(ax1,rmax_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{2},'MarkerSize',12)
%plot(ax1,r0Lil_Lilctrl_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{4},'MarkerSize',12)
plot(ax1,r0Lil_tau_equil_g,0,'x-.','MarkerEdgeColor',pl_clrs{3},'MarkerSize',12)

if(equil_dynamic == 1)
    t_equil = t0_equil:dt/60/60/24:tf_equil;
    plot(ax1,t_equil,ones(length(t_equil),1),'m','LineWidth',2)
else
    t_equil = tf-dt_final:dt/60/60/24:tf;
    plot(ax1,t_equil,ones(length(t_equil),1),'m','LineWidth',2)
end

input_legend = {'$V_m$','$r_m$','$r_0$'};
legend(input_legend,'Location','SouthEast','Interpreter','Latex')
input_title1=sprintf('Evolution of control storm structure');
input_title2=sprintf('%s',subdir);
%title({input_title1},'Linestyle','none','Interpreter','none')
xlabel('Time [day]','FontSize',18)
ylabel('Non-dimensional units','FontSize',18)
%ymax = 1*max([max(Vm_plot) max(rm_plot(Vm_plot>.7)) max(r0_plot(Vm_plot>.7))])
ymax = 2;   %just make plot go [0,2] in y axis for simplicity
axis([t0 tf 0 min([ymax 2])])
set(ax1,'YTick',[0 .2 .4 .6 .8 1 1.2 1.4 1.6 1.8 2],'XTick',[0 20 40 60 80 100 120 140])
grid on
text1=text(.73*tf,.97*min([ymax 2]),sprintf('$V^*_m = $ %5.1f $ms^{-1}$',Vmax_equil_g),'fontweight','bold','FontSize',14,'Color','k');
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
text2=text(.73*tf,.89*min([ymax 2]),sprintf('$r^*_m = $ %5.1f $km$',rmax_equil_g),'fontweight','bold','FontSize',14,'Color','k');
set(text2,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');
text3=text(.73*tf,.82*min([ymax 2]),sprintf('$r^*_0 = $ %5.1f $km$',r0Lil_equil_g),'fontweight','bold','FontSize',14,'Color','k');
set(text3,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','BackgroundColor','white');

%text5=text(-.1*tf,1.13*min([ymax 2]),sprintf('T_{mean} = %i days',T_mean),'fontweight','bold','FontSize',12,'Color','k');
%set(text5,'HorizontalAlignment','left','VerticalAlignment','top');

%text6=text(.9*tf,1.10*min([ymax 2]),sprintf('V_p = %3.0f ms^{-1}',mpi),'fontweight','bold','FontSize',12,'Color','k');
%set(text6,'HorizontalAlignment','left','VerticalAlignment','top');

f_exp = floor(log10(fcor));
f_coef = fcor/(10^f_exp);
%text7=text(.9*tf,1.06*min([ymax 2]),sprintf('f = %1.0f*10^{%2.0f} s^{-1}',f_coef,f_exp),'fontweight','bold','FontSize',12,'Color','k');
%set(text7,'HorizontalAlignment','left','VerticalAlignment','top');

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0


%save PDF
print -dpdf -r300 Fig2_Control_run.pdf
