%TC_stats_set.m

%Created: 11-29-11, Dan Chavas

%Purpose: Plot sensitivity of each of Vmax, rmax, r0 to ALL parameters in
%single plot

clear
clc
figure(1)
clf(1)
figure(2)
clf(2)

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold')

%clf(1)
%clf(2)
%clf(3)

%%variables of interest (sim_set name): 'dx' 'dz' 'domain' 'lh' 'lv' 'H' 'Qrad' 'Vpot' 'cor' 'qro' 'ro' 'rodrmax'
%sim_sets = {'dx' 'dz' 'domain' 'lh' 'lv' 'H' 'Qrad' 'Vpot' 'cor' 'qro' 'ro' 'rodrmax'};
sim_sets = {'lh' 'lv' 'cor' 'qro' 'ro' 'mpi'}
%sim_sets = {'cor'}
%sim_sets = {'test'}

dat_max=0;
dat_min = 0;


for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('simsets/%s.mat',sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    
    %xvals_pl = CTRL_val*multipliers;    %values defined by user at top
    xvals_pl = multipliers;    %values defined by user at top
    i_ctrl = find(multipliers==0,1);

    figure(1)
    subplot(2,2,1)
    data_temp = Vmax_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,2)
    data_temp = rmax_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,3)
    data_temp = r0Lil_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on
    
    subplot(2,2,4)
    data_temp = Vmax_equil_g./rmax_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    
    figure(2)
    subplot(2,2,1)
    data_temp = Vmax_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,2)
    data_temp = rmax_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,3)
    data_temp = r0Lil_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    
%{
    figure(2)
    for pp=1:length(data_tmean_usr_g_all)
        data_temp(pp) = xvals(find(smooth(smooth(data_tmean_usr_g_all{pp},200),200)<0,1));
    end
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on
%}
    
%    cd simsets
%    saveas(gcf,'Vmax','jpeg')
%    cd ../../..


end

figure(1)
subplot(2,2,1)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = Vmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Equilibrium storm: scaling of Vmax');
input_title2=sprintf('Vmax* = %5.2f ms^{-1}',Vmax_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,2)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = rmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Equilibrium storm: scaling of rmax');
input_title2=sprintf('rmax* = %5.2f km',rmax_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,3)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = r0')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Equilibrium storm: scaling of r0');
input_title2=sprintf('r0* = %5.2f km',r0Lil_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,4)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = dV/dr in eye')
xlabel({sprintf('log_2(X/X*)')})

legend(sim_sets,'Location','SouthEast')
input_title1=sprintf('GRADIENT Equilibrium storm: scaling of dV/dr in eye');
input_title2=sprintf('(dV/dr)* = %5.2f ms^{-1} km^{-1}',Vmax_equil(i_ctrl)/rmax_equil(i_ctrl));
title({input_title1,input_title2})
grid on


figure(2)
subplot(2,2,1)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = Vmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Equilibrium storm: scaling of Vmax timescale');
input_title2=sprintf('Vmax* = %5.2f ms^{-1}',Vmax_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,2)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = rmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Equilibrium storm: scaling of rmax timescale');
input_title2=sprintf('rmax* = %5.2f km',rmax_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,3)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = r0')
xlabel({sprintf('log_2(X/X*)')})

legend(sim_sets,'Location','SouthEast')
input_title1=sprintf('GRADIENT Equilibrium storm: scaling of r0 timescale');
input_title2=sprintf('r0* = %5.2f km',r0Lil_equil_g(i_ctrl));
title({input_title1,input_title2})
grid on










dat_max=0;
dat_min = 0;


for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('simsets/%s.mat',sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    
    %xvals_pl = CTRL_val*multipliers;    %values defined by user at top
    xvals_pl = multipliers;    %values defined by user at top
    i_ctrl = find(multipliers==0,1);

    figure(3)
    subplot(2,2,1)
    data_temp = Vmax_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,2)
    data_temp = rmax_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    subplot(2,2,3)
    data_temp = r0Lil_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on
    
    subplot(2,2,4)
    data_temp = Vmax_gen_g./rmax_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on

    
  
%{
    figure(2)
    for pp=1:length(data_tmean_usr_g_all)
        data_temp(pp) = xvals(find(smooth(smooth(data_tmean_usr_g_all{pp},200),200)<0,1));
    end
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    hold on
%}
    
%    cd simsets
%    saveas(gcf,'Vmax','jpeg')
%    cd ../../..


end

figure(3)
subplot(2,2,1)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = Vmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Genesis storm: scaling of Vmax');
input_title2=sprintf('Vmax* = %5.2f ms^{-1}',Vmax_gen_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,2)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = rmax')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Genesis storm: scaling of rmax');
input_title2=sprintf('rmax* = %5.2f km',rmax_gen_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,3)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = r0')
xlabel({sprintf('log_2(X/X*)')})

input_title1=sprintf('GRADIENT Genesis storm: scaling of r0');
input_title2=sprintf('r0* = %5.2f km',r0Lil_gen_g(i_ctrl));
title({input_title1,input_title2})
grid on

subplot(2,2,4)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)  Y* = dV/dr in eye')
xlabel({sprintf('log_2(X/X*)')})

legend(sim_sets,'Location','SouthEast')
input_title1=sprintf('GRADIENT Genesis storm: scaling of dV/dr in eye');
input_title2=sprintf('(dV/dr)* = %5.2f ms^{-1} km^{-1}',Vmax_gen(i_ctrl)/rmax_gen(i_ctrl));
title({input_title1,input_title2})
grid on

