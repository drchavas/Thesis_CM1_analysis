%TC_stats_plot.m

%Created: 10-05-11, Dan Chavas
%Updated: 12-02-11, Dan Chavas

%Purpose: This file makes lots of plots of the stats data for a desired set
%of simulations whose individual files are located in stats_files/.
%Output plots are saved in simsets_Tmean#_#_#/PLOTS/[sim_set]/.

%%TESTING
%{
run_type = 1;
T_mean = 5;
equil_dynamic = 0;
dt_final = 50;
tf = 150;
dt_final_dynamic = 30;
rmin_plot = 0
rmax_plot = 1500;
CTRL_val = 1; %CTRL value of quantity varied across simulations
units = '-';
subdirs_set = {
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4_lh375'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh187.5'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv4_lh750'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4_lh750'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh187.5'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv4'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_lh187.5'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv4_lh3000'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fx4_lh187.5'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4_lh3000'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh187.5'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4_lh6000'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh375'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fx4_lh750'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh750'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fx4'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh3000'
    'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh6000'
    };
sim_set = 'nondim';
multipliers = ones(length(subdirs_set),1);

%}
function [junk] = TC_stats_plot(run_type,T_mean,equil_dynamic,dt_final,tf,dt_final_dynamic,rmin_plot,rmax_plot,CTRL_val,units,multipliers,subdirs_set,sim_set);
junk='junk';
%clear
%clc        

%% USER INPUT %%%%%%%%%%%%%%%%%%
wrad_const = 0; %1 = use CTRL value for wrad

%% Plots on/off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_gradient = 1;  %1=make plots with gradient wind
plot_full = 0;  %1=make plots with full wind (note: can do both)

%SINGLE/MULTI SIMULATION
plot_usrprof = 1;   %0=no; plot radial wind profile for user-defined time-means

%MULTI SIMULATION
plot_ts_multi = 1;    %0=no; plot time series of vmax, rmax, rmid, r0
plot_stats = 1; %0=no; plot comparisons of basic statistics: Vmax, Rmax, R0, tau_equil for each; can be done at "equil time", genesis time, time of first V~Vmax, final (day 70-100 ave)

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    
numruns=length(subdirs_set);  %total number of runs you want to plot (taken from below)

%%Append ax or 3d to subdir names for plot
for i=1:length(subdirs_set)
    if(run_type==1) %ax
        subdirs_load{i}=sprintf('ax%s',subdirs_set{i});
    else    %3d
        subdirs_load{i}=sprintf('3d%s',subdirs_set{i});
    end
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;
clear Vmax_tau_equil rmax_tau_equil rmid_tau_equil r0_tau_equil r0Lil_tau_equil Vmax_tau_equil_g rmax_tau_equil_g rmid_tau_equil_g r0_tau_equil_g r0Lil_tau_equil_g
clear tau_gen Vmax_gen_g rmax_gen_g rmid_gen_g r0Lil_gen_g Vmax_tau_max_g Vmax_max_g rmax_tau_max_g rmax_max_g r0_tau_max_g r0_max_g r0Lil_tau_max_g r0Lil_max_g
clear Vmax_equil rmax_equil rmid_equil r0_equil r0Lil_equil Vmax_equil_g rmax_equil_g rmid_equil_g r0_equil_g r0Lil_equil_g
clear Vmax_movave_g_all rmax_movave_g_all r0_movave_g_all r0Lil_movave_g_all
clear xvals_sub_all data_tmean_usr_g_all
clear mpi_all fcor_all
for ss=1:numruns

    %%Load data for given simulation
    if(equil_dynamic==1)
        if(wrad_const == 1)
            load(sprintf('simdata_Tmean%i_dt%i_dynamic_wradconst/%s.mat',T_mean,dt_final_dynamic,subdirs_load{ss}));
        else
            load(sprintf('simdata_Tmean%i_dt%i_dynamic/%s.mat',T_mean,dt_final_dynamic,subdirs_load{ss}));
        end
    else
        if(wrad_const == 1)
            load(sprintf('simdata_Tmean%i_%i_%i_wradconst/%s.mat',T_mean,tf-dt_final,tf,subdirs_load{ss}));
        else
            load(sprintf('simdata_Tmean%i_%i_%i/%s.mat',T_mean,tf-dt_final,tf,subdirs_load{ss}));
        end
    end
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    rmid_equil(ss) = rmid_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    rmid_equil_g(ss) = rmid_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    
    %timescales to those values
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    rmid_tau_equil(ss) = rmid_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    rmid_tau_equil_g(ss) = rmid_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
%    Vmax_gen(ss) = Vmax_gen_sim;
%    rmax_gen(ss) = rmax_gen_sim;
%    rmid_gen(ss) = rmid_gen_sim;
%    r0_gen(ss) = r0_gen_sim;
%    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    rmid_gen_g(ss) = rmid_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    rmid_tau_max_g(ss) = rmid_tau_max_g_sim; %defined using the GRADIENT wind
    rmid_max_g(ss) = rmid_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
%    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
%    Vmax_max(ss) = Vmax_max_sim;
%    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
%    rmax_max(ss) = rmax_max_sim;
%    rmid_tau_max(ss) = rmid_tau_max_sim; %defined using the GRADIENT wind
%    rmid_max(ss) = rmid_max_sim;
%    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
%    r0_max(ss) = r0_max_sim;
%    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
%    r0Lil_max(ss) = r0Lil_max_sim;
    
    %% Save data for all simulations
%    Vmax_movave_all(:,ss)=Vmax_movave_sim;
%    rmax_movave_all(:,ss)=rmax_movave_sim;
%    rmid_movave_all(:,ss)=rmid_movave_sim;
%    r0_movave_all(:,ss)=r0_movave_sim;
%    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;

    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    rmid_movave_g_all(:,ss)=rmid_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
%    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;

    %% Extract and keep mpi and fcor%%%%%
    mpi_all(ss) = mpi;
    fcor_all(ss) = fcor;
end

%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out = sprintf('simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out = sprintf('simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out = sprintf('simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out = sprintf('simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

%% SAVE SIMULATION SET DATA AS simsets_Tmean#_#_#/[sim_set].mat %%%
save temp.mat
movefile('temp.mat',sprintf('%s/%s.mat',subdir_out,sim_set))

%% PLOTTING %%%%%%%%%%%%%%%%
mkdir(sprintf('%s/PLOTS/%s',subdir_out,sim_set))

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
%Single/Multi simulation: Plot user-defined radial wind profiles
if(plot_usrprof==1)

    if(plot_full == 1)
    figure(4)
    set(gcf, 'Visible', 'off') 
    clf(4)
    
    hold off
    for i=1:numruns
        i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
        xvals_pl = xvals_sub_all{i}(i_xvals_pl);

        plot(xvals_pl,data_tmean_usr_all{i}(i_xvals_pl),pl_clrs{i},'LineWidth',1)
        hold on
    end
    plot(0*(1:xvals_pl(end)),'k')
    input_title = sprintf('Time-mean radial FULL wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
    title(input_title)
    xlabel('Radius [km]')
    ylabel('Azimuthal wind speed [m/s]')
    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
    legend(input_legend)
    grid on
    end

    if(plot_gradient == 1)
    figure(5)
    set(gcf, 'Visible', 'off') 
    clf(5)
    
    hold off
    for i=1:numruns
        i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
        xvals_pl = xvals_sub_all{i}(i_xvals_pl);

        %%DRC 01/31/12
%        angmom = 1000*xvals_pl'.*data_tmean_usr_all{i}(i_xvals_pl) + .5*fcor*(1000*xvals_pl').^2;
%        dangmomdr = (angmom(3:end)-angmom(1:end-2))/(2*dx);
        %plot(xvals_pl,angmom,pl_clrs{i},'LineWidth',1)
%        plot(xvals_pl(2:end-1),dangmomdr,pl_clrs{i},'LineWidth',1)
%        angmomloss(i) = sum(dangmomdr)
        %%DRC
        
        
        plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl'),pl_clrs{i},'LineWidth',1)
        hold on
        
    end
    plot(0*(1:xvals_pl(end)),'k')
    input_title = sprintf('Time-mean radial GRADIENT wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
    title(input_title)
    xlabel('Radius [km]')
    ylabel('Azimuthal gradient wind speed [m/s]')
    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
    legend(input_legend)
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'radprof','jpeg')
    cd ../../..
    
    %%SAME PLOT BUT WITH Vm SCALED BY Vp AND RADII SCALED BY Vp/f
    figure(51)
    set(gcf, 'Visible', 'off') 
    clf(51)
    
    hold off
    for i=1:numruns
        i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
        xvals_pl = xvals_sub_all{i}(i_xvals_pl);

        %%DRC 01/31/12
%        angmom = 1000*xvals_pl'.*data_tmean_usr_all{i}(i_xvals_pl) + .5*fcor*(1000*xvals_pl').^2;
%        dangmomdr = (angmom(3:end)-angmom(1:end-2))/(2*dx);
        %plot(xvals_pl,angmom,pl_clrs{i},'LineWidth',1)
%        plot(xvals_pl(2:end-1),dangmomdr,pl_clrs{i},'LineWidth',1)
%        angmomloss(i) = sum(dangmomdr)
        %%DRC
        
        
        plot(xvals_pl/(mpi_all(i)/fcor),data_tmean_usr_g_all{i}(i_xvals_pl)/mpi_all(i),pl_clrs{i},'LineWidth',1)
        hold on
        mpi_all(i)
    end
    plot(0*(.01:.01:xvals_pl(end)/(min(mpi_all(i))/fcor)),'k')
    input_title = sprintf('Time-mean radial GRADIENT wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
    title(input_title)
    xlabel('r/(V_p/f)')
    ylabel('V_g/V_p')
    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
    legend(input_legend)
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'radprof_nondim','jpeg')
    cd ../../..
    
    
    end
    
end

%Multi simulation: Plot comparison of various statistical quantities
if(plot_stats==1)

    %xvals_pl = CTRL_val*multipliers;    %values defined by user at top
    xvals_pl = multipliers;    %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    if(isempty(i_ctrl))
        i_ctrl = 1;
    end

    if(strcmp(sim_set,'mpi'))
        i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
        CTRL_val = mpi_all(i_ctrl)
        multipliers = log2(mpi_all/CTRL_val);  %log2(1.5) = x1.5;
        xvals_pl = multipliers;
    end
    
    %%Equilibrium storm characteristics
    if(plot_full == 1)
    figure(6)
    set(gcf, 'Visible', 'off') 
    clf(6)
    
%    subplot(2,1,2)
    hold off
    dat_max=0;
    dat_min=1000;

    data_temp = Vmax_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl)); 
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

%    data_temp = rmid_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = r0Lil_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on

    %axis([min(xvals_pl) max(xvals_pl) min(.9*dat_min,1.1*dat_min) 1.1*dat_max])
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})

    legend({'Vmax','rmax','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_final,tf,sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_equil(i_ctrl),rmax_equil(i_ctrl),r0Lil_equil(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    figure(61)
    set(gcf, 'Visible', 'off') 
    clf(61)
    
%    subplot(2,1,2)
    hold off
    dat_max=0;
    dat_min=1000;
    
    data_temp = Vmax_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-.')
    hold on
    
    data_temp = rmax_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-.')
    hold on

%    data_temp = rmid_tau_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-.')
%    hold on

%    data_temp = r0_tau_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-.')
%    hold on
    
    data_temp = r0Lil_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-.')
    hold on
    
    %axis([min(xvals_pl) max(xvals_pl) min(.9*dat_min,1.1*dat_min) 1.1*dat_max])
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})

    legend({'tauvmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_final,tf,sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_equil(i_ctrl),rmax_tau_equil(i_ctrl),r0Lil_tau_equil(i_ctrl));
    title({input_title1,input_title3})
    grid on
    end

    if(plot_gradient == 1)
    figure(7)
    set(gcf, 'Visible', 'off') 
    clf(7)
    
%    subplot(2,1,2)
    hold off
    dat_max=0;
    dat_min = 0;

    data_temp = Vmax_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

%    data_temp = rmid_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = rmid_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on

    data_temp = r0Lil_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'gx-')
    hold on

    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
    legend({'Vmax','rmax','rmid','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_final,tf,sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; rmid* = %5.2f km; r0Lil* = %5.2f km',Vmax_equil_g(i_ctrl),rmax_equil_g(i_ctrl),rmid_equil_g(i_ctrl),r0Lil_equil_g(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'equil','jpeg')
    cd ../../..
    
    %%Timescales to equilibrium
    figure(71)
    set(gcf, 'Visible', 'off') 
    clf(71)
    
%    subplot(2,1,2)
    hold off
    dat_max=0;
    dat_min = 0;

    data_temp = Vmax_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-.')
    hold on
    
    data_temp = rmax_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-.')
    hold on

%    data_temp = rmid_tau_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-.')
%    hold on

%    data_temp = r0_tau_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-.')
%    hold on
    
    data_temp = rmid_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-.')
    hold on
    
    data_temp = r0Lil_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'gx-.')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
    legend({'tauvmax','taurmax','taurmid','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_final,tf,sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taurmid* = %5.2f d; taur0* = %5.2f d',Vmax_tau_equil_g(i_ctrl),rmax_tau_equil_g(i_ctrl),rmid_tau_equil_g(i_ctrl),r0Lil_tau_equil_g(i_ctrl));
    title({input_title1,input_title3})
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'tau_equil','jpeg')
    cd ../../..

    end
    
    %%Transient timescales
    if(plot_full == 1)
    figure(8)
    set(gcf, 'Visible', 'off') 
    clf(8)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;

    data_temp = Vmax_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

%    data_temp = rmid_tau_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_tau_equil;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = r0Lil_tau_equil;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on

    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})

    legend({'Vmax','rmax','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Equilibration timescales: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f day; rmax* = %5.2f day; r0Lil* = %5.2f day',Vmax_tau_equil(i_ctrl),rmax_tau_equil(i_ctrl),r0Lil_tau_equil(i_ctrl));
    title({input_title1,input_title2})
    grid on
    end
    
    %%Genesis storm characeristics
    if(plot_full == 1)
    figure(10)
    set(gcf, 'Visible', 'off') 
    clf(10)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;

    data_temp = Vmax_gen;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_gen;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

%    data_temp = rmid_gen;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_gen;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = r0Lil_gen;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on
    
%    data_temp = r0Lil_gen;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'mx-')
%    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
	xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax','rmax','rmid','r0','r0Lil'},'Location','SouthEast')
    legend({'Vmax','rmax','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL storm at Genesis: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_gen(i_ctrl),rmax_gen(i_ctrl),r0Lil_gen(i_ctrl));
    input_title3=sprintf('Y*: taugen* = %5.2f',tau_gen(i_ctrl));
    title({input_title1,input_title2,input_title3})
    grid on
    end
    
    if(plot_gradient == 1)
    figure(11)
    set(gcf, 'Visible', 'off') 
    clf(11)
    
%    subplot(2,1,2)
    hold off
    dat_max=0;
    dat_min = 1000;

    data_temp = Vmax_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

%    data_temp = rmid_gen_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_gen_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = rmid_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on

    data_temp = r0Lil_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'gx-')
    hold on
    
%    data_temp = r0Lil_gen_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'mx-')
%    hold on

    data_temp = tau_gen;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-.')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'Vmax','rmax','rmid','r0Lil','taugen'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT storm at Genesis: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; rmid* = %5.2f km; r0Lil* = %5.2f km',Vmax_gen_g(i_ctrl),rmax_gen_g(i_ctrl),rmid_gen_g(i_ctrl),r0Lil_gen_g(i_ctrl));
    input_title3=sprintf('Y*: taugen* = %5.2f',tau_gen(i_ctrl));
    title({input_title1,input_title2,input_title3})
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'gen','jpeg')
    cd ../../..

    end
    
    %%Peak values of Vmax, rmax, r0; time scales to each
    if(plot_full == 1)
    figure(12)
    set(gcf, 'Visible', 'off') 
    clf(12)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;
    
    data_temp = Vmax_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

    data_temp = r0Lil_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'kx-')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'Vmaxmax','rmaxmax','r0Lilmax'},'Location','SouthEast')
    input_title1=sprintf('FULL Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_max(i_ctrl),rmax_max(i_ctrl),r0Lil_max(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    figure(121)
    set(gcf, 'Visible', 'off') 
    clf(121)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;
    
    data_temp = Vmax_tau_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-.')
    hold on
    
    data_temp = rmax_tau_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-.')
    hold on
    
    data_temp = r0Lil_tau_max;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-.')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'tauVmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_max(i_ctrl),rmax_tau_max(i_ctrl),r0Lil_tau_max(i_ctrl));
    title({input_title1,input_title3})
    grid on
    end
    
    if(plot_gradient == 1)
    figure(13)
    set(gcf, 'Visible', 'off') 
    clf(13)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;
    
    data_temp = Vmax_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-')
    hold on
    
    data_temp = rmax_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-')
    hold on

    data_temp = rmid_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on
    
    data_temp = r0Lil_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'gx-')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'Vmaxmax','rmaxmax','rmidmax','r0Lilmax'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; rmid* = %5.2f km; r0Lil* = %5.2f km',Vmax_max_g(i_ctrl),rmax_max_g(i_ctrl),rmid_max_g(i_ctrl),r0Lil_max_g(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'max','jpeg')
    cd ../../..

    figure(131)
    set(gcf, 'Visible', 'off') 
    clf(131)
    
%    subplot(2,1,1)
    hold off
    dat_max = 0;
    dat_min = 1000;

    data_temp = Vmax_tau_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'rx-.')
    hold on
    
    data_temp = rmax_tau_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'bx-.')
    hold on
    
    data_temp = rmid_tau_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-.')
    hold on
    
    data_temp = r0Lil_tau_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'gx-.')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'tauVmax','taurmax','taurmid','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taurmid* = %5.2f d; taur0* = %5.2f d',Vmax_tau_max_g(i_ctrl),rmax_tau_max_g(i_ctrl),rmid_tau_max_g(i_ctrl),r0Lil_tau_max_g(i_ctrl));
    title({input_title1,input_title3})
    grid on
    
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'tau_max','jpeg')
    cd ../../..

    end
    
end


%Multi simulation: plot time-series of each variable
if(plot_ts_multi==1)  

    if(plot_full == 1)
%% FULL WIND %%%%%%%%%%%    
    %Plot time series
    figure(14)
    set(gcf, 'Visible', 'off') 
    clf(14)
    
    
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(Vmax_movave_all))])

    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),Vmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max(i)))
            h=plot(t_day(Vmax_tau_max(i)/(dt/60/60/24)),Vmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil(i)))
            h=plot(t_day(Vmax_tau_equil(i)/(dt/60/60/24)),Vmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(rmax_movave_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),rmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max(i)))
            h=plot(t_day(rmax_tau_max(i)/(dt/60/60/24)),rmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil(i)))
            h=plot(t_day(rmax_tau_equil(i)/(dt/60/60/24)),rmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
%{
    %rmid
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day,rmid_movave_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(rmid_tau_equil(i)))
            h=plot(t_day(rmid_tau_equil(i)/(dt/60/60/24)),rmid_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('rmid [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(rmid_movave_all(100:end,:)))])
%}
%{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day,r0_movave_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil(i)))
            h=plot(t_day(r0_tau_equil(i)/(dt/60/60/24)),r0_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('r0 [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(r0_movave_all(100:end,:)))])
%}
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day,r0Lil_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(r0Lil_movave_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),r0Lil_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max(i)))
            h=plot(t_day(r0Lil_tau_max(i)/(dt/60/60/24)),r0Lil_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil(i)/(dt/60/60/24)),r0Lil_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
        
    hold on
    input_title1=sprintf('FULL storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    
%% Zoom in on genesis period; plot from 0:tau_gen
    figure(16)
    set(gcf, 'Visible', 'off') 
    clf(16)
    
    
    for i=1:numruns
        i_gen(i) = find(t_day == tau_gen(i));
    end
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day(1:i_gen(i))/tau_gen(i),Vmax_movave_all(1:i_gen(i),i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
%    xlabel('t/taugen')
    ylabel(input_title)
    axis([0 1 0 1.1*max(max(Vmax_movave_all))])

    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day(1:i_gen(i))/tau_gen(i),rmax_movave_all(1:i_gen(i),i),pl_clrs{i})
        hold on

    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
%    xlabel('t/taugen')
    ylabel(input_title)
    axis([0 1 0 1.1*max(max(rmax_movave_all(20:end,:)))])

    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day(1:i_gen(i))/tau_gen(i),r0Lil_movave_all(1:i_gen(i),i),pl_clrs{i})
        hold on

    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    axis([0 1 0 1.1*max(max(r0Lil_movave_all(20:end,:)))])
    
    hold on
    input_title1=sprintf('FULL pre-genesis evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    xlabel('t/taugen')
    ylabel(input_title)
    suptitle({input_title1,input_title2})
    grid on
    end
    
    if(plot_gradient == 1)
%% GRADIENT WIND %%%%%%%%

%%Plot time series
    figure(15)
    set(gcf, 'Visible', 'off') 
    clf(15)
    
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_g_all(:,i),pl_clrs{i})
        hold on
        
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(Vmax_movave_g_all))])

    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),Vmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max_g(i)))
            h=plot(t_day(Vmax_tau_max_g(i)/(dt/60/60/24)),Vmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_g_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(rmax_movave_g_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),rmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max_g(i)))
            h=plot(t_day(rmax_tau_max_g(i)/(dt/60/60/24)),rmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil_g(i)))
            h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
%{
    %rmid
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day,rmid_movave_g_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(rmid_tau_equil_g(i)))
            h=plot(t_day(rmid_tau_equil_g(i)/(dt/60/60/24)),rmid_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('rmid [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(rmid_movave_g_all(100:end,:)))])
%}
%{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day,r0_movave_g_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil_g(i)))
            h=plot(t_day(r0_tau_equil_g(i)/(dt/60/60/24)),r0_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('r0 [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(r0_movave_g_all(100:end,:)))])
%}
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day,r0Lil_movave_g_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(r0Lil_movave_g_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),r0Lil_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max_g(i)))
            h=plot(t_day(r0Lil_tau_max_g(i)/(dt/60/60/24)),r0Lil_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
        
    hold on
    input_title1=sprintf('GRADIENT storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    set(gcf, 'Visible', 'off') 
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'ts','jpeg')
    cd ../../..
    
%%Plot non-dimensional time series
    figure(151)
    set(gcf, 'Visible', 'off') 
    clf(151)
    
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_g_all(:,i)/mpi_all(i),pl_clrs{i})
        hold on
        
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('V_m / V_p');
    %title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(Vmax_movave_g_all/mpi_all(i)))])

    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),Vmax_gen_g(i)/mpi_all(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max_g(i)))
            h=plot(t_day(Vmax_tau_max_g(i)/(dt/60/60/24)),Vmax_max_g(i)/mpi_all(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i)/mpi_all(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r_m / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
%    axis([0 tf 0 1.1*max(max(rmax_movave_g_all(20:end,:)/(mpi_all(i)/fcor_all(i)/1000)))])
    axis([0 tf 0 .07])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),rmax_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max_g(i)))
            h=plot(t_day(rmax_tau_max_g(i)/(dt/60/60/24)),rmax_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil_g(i)))
            h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
%{
    %rmid
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day,rmid_movave_g_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(rmid_tau_equil_g(i)))
            h=plot(t_day(rmid_tau_equil_g(i)/(dt/60/60/24)),rmid_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('rmid [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(rmid_movave_g_all(100:end,:)))])
%}
%{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day,r0_movave_g_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil_g(i)))
            h=plot(t_day(r0_tau_equil_g(i)/(dt/60/60/24)),r0_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('r0 [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(r0_movave_g_all(100:end,:)))])
%}
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day,r0Lil_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r_0 / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
%    axis([0 tf 0 1.1*max(max(r0Lil_movave_g_all(20:end,:)/(mpi_all(i)/fcor_all(i)/1000)))])
    axis([0 tf 0 1.1])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),r0Lil_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max_g(i)))
            h=plot(t_day(r0Lil_tau_max_g(i)/(dt/60/60/24)),r0Lil_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
        
    hold on
    input_title1=sprintf('GRADIENT non-dimensional storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    set(gcf, 'Visible', 'off') 
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'ts_nondim','jpeg')
    cd ../../..
    
%% Zoom in on genesis period; plot from 0:tau_gen
    figure(17)
    set(gcf, 'Visible', 'off') 
    clf(17)
    
    for i=1:numruns
        i_gen(i) = find(t_day == tau_gen(i));
    end
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_load(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    dat_max = 0;
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day(1:i_gen(i))/tau_gen(i),Vmax_movave_g_all(1:i_gen(i),i),pl_clrs{i})
        dat_max = max(dat_max,max(Vmax_movave_g_all(1:i_gen(i),i)));
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
%    xlabel('t/taugen')
    ylabel(input_title)
    axis([0 1 0 1.1*dat_max])

    %rmax
    dat_max = 0;
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day(1:i_gen(i))/tau_gen(i),rmax_movave_g_all(1:i_gen(i),i),pl_clrs{i})
        dat_max = max(dat_max,max(rmax_movave_g_all(1:i_gen(i),i)));
        hold on

    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
%    xlabel('t/taugen')
    ylabel(input_title)
    axis([0 1 0 1.1*dat_max])

    %r0Lil
    dat_max = 0;
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day(1:i_gen(i))/tau_gen(i),r0Lil_movave_g_all(1:i_gen(i),i),pl_clrs{i})
        dat_max = max(dat_max,max(r0Lil_movave_g_all(1:i_gen(i),i)));
        hold on

    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    axis([0 1 0 1.1*dat_max])
    
    hold on
    input_title1=sprintf('GRADIENT pre-genesis evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    xlabel('t/taugen')
    ylabel(input_title)
    suptitle({input_title1,input_title2,input_title3})
    set(gcf, 'Visible', 'off') 
    cd(sprintf('%s/PLOTS/%s',subdir_out,sim_set))
    saveas(gcf,'pregen','jpeg')
    cd ../../..

    end
end

end
