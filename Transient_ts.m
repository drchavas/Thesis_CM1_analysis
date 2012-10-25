%Transient_ts.m

%Created: 2 Oct 2012, Dan Chavas

%Purpose: To plot time series, both dimensional and non-dimensional, of
%structural variables for various simulations



clear
clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
run_type = 1;
wrad_const = 0;
T_mean = 2;
equil_dynamic = 1;
dt_final = 30;
tf = 100;
dt_final_dynamic = 10;
CTRL_val = 1; %CTRL value of quantity varied across simulations
units = '-';    %for simset
subdirs_set = {
    'CTRLv0qrhSATqdz5000_nx3072_drag'   %150 days
    'CTRLv0qrhSATqro110000qdz5000_nx3072_Tthresh250K_lh825_drag'    %150 days
    'CTRLv0qrhSATqro400000qdz5000_nx3072_fdiv2_lh3000_drag' %100 days
    'CTRLv0qrhSATqro100000qdz5000_nx3072_fx2_lh750_drag'    %100 days
    'CTRLv0qrhSATqro296000qdz5000_nx3072_usfc.5_lh2220_drag'    %150 days
    };
legend_info = {'CTRL' 'qro110Tthresh250lh825' 'qro400fdiv2lh3000' 'qro100fx2lh750' 'qro296usfc.5lh2220'}
sim_set = 'transient';
multipliers = ones(length(subdirs_set),1);
Vmax_pl = 100;
rmmax_pl = 150;
r0max_pl = 3500;
Vmax_nd_pl = 1.5;
rmmax_nd_pl = .1;
r0max_nd_pl = 1.5;
%}

time_rescale = [1 .55 1/2 1 1]; %rescale simulation time by these factors
%time_rescale = ones(length(subdirs_set),1); %rescale simulation time by these factors

numday_smooth = 0;  %CHECK ME!! [day]; 0 = no smoothing;
%%%%%%%%%%%%%%%%%%%%%%%%

dir_in_dat = sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/CM1_postproc_data/simdata_Tmean%i_%i_%i/',T_mean,tf-dt_final,dt_final);
dir_in_dat_dyn = sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic/',T_mean,dt_final_dynamic);

%% Plots on/off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_gradient = 1;  %1=make plots with gradient wind
plot_full = 0;  %1=make plots with full wind (note: can do both)

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


t_day_min = 0;
t_day_max = 0;
clear Vmax_tau_equil rmax_tau_equil rrad_tau_equil r0_tau_equil r0ER11_tau_equil r0Lil_tau_equil Vmax_tau_equil_g rmax_tau_equil_g rrad_tau_equil_g r0_tau_equil_g r0ER11_tau_equil_g r0Lil_tau_equil_g
clear tau_gen Vmax_gen_g rmax_gen_g rrad_gen_g r0ER11_gen_g r0Lil_gen_g Vmax_tau_max_g Vmax_max_g rmax_tau_max_g rmax_max_g r0_tau_max_g r0_max_g r0ER11_tau_max_g r0ER11_max_g r0Lil_tau_max_g r0Lil_max_g
clear Vmax_equil rmax_equil rrad_equil r0_equil r0ER11_equil r0Lil_equil Vmax_equil_g rmax_equil_g rrad_equil_g r0_equil_g r0ER11_equil_g r0Lil_equil_g
clear Vmax_movave_g_all rmax_movave_g_all rrad_movave_g_all r0_movave_g_all r0ER11_movave_g_all r0Lil_movave_g_all
clear xvals_sub_all data_tmean_usr_g_all
clear mpi_all fcor_all
for ss=1:numruns

    %%Load data for given simulation
    if(equil_dynamic==1)
        if(wrad_const == 1)
            load(sprintf('%s_wradconst/%s.mat',dir_in_dat_dyn,subdirs_load{ss}));
        else
            load(sprintf('%s/%s.mat',dir_in_dat_dyn,subdirs_load{ss}));
        end
    else
        if(wrad_const == 1)
            load(sprintf('%s_wradconst/%s.mat',dir_in_dat,subdirs_load{ss}));
        else
            load(sprintf('%s/%s.mat',dir_in_dat,subdirs_load{ss}));
        end
    end
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(ss) = Vmax_equil_sim;
    rmax_equil(ss) = rmax_equil_sim;
    rrad_equil(ss) = rrad_equil_sim;
    r0_equil(ss) = r0_equil_sim;
    r0Lil_equil(ss) = r0Lil_equil_sim;
    Vmax_equil_g(ss) = Vmax_equil_g_sim;
    rmax_equil_g(ss) = rmax_equil_g_sim;
    rrad_equil_g(ss) = rrad_equil_g_sim;
    r0_equil_g(ss) = r0_equil_g_sim;
    r0ER11_equil_g(ss) = r0ER11_equil_g_sim;
    r0Lil_equil_g(ss) = r0Lil_equil_g_sim;
    r0Lil_Lilctrl_equil_g(ss) = r0Lil_Lilctrl_equil_g_sim;
    
    if(equil_dynamic==1)
        %equilibrium values based on independent equilibration of each variable
        Vmax_equil_g_dynamic_all(ss) = Vmax_equil_g_dynamic;
        rmax_equil_g_dynamic_all(ss) = rmax_equil_g_dynamic;
        rrad_equil_g_dynamic_all(ss) = rrad_equil_g_dynamic;
        r0_equil_g_dynamic_all(ss) = r0_equil_g_dynamic;
        r0Lil_equil_g_dynamic_all(ss) = r0Lil_equil_g_dynamic;
        r0Lil_Lilctrl_equil_g_dynamic_all(ss) = r0Lil_Lilctrl_equil_g_dynamic;
    
        %equilibrium start times associated with above independent equilibration values
        t0_equil_Vmax_all(ss) = t0_equil_Vmax;
        t0_equil_rmax_all(ss) = t0_equil_rmax;
        t0_equil_rrad_all(ss) = t0_equil_rrad;
        t0_equil_r0_all(ss) = t0_equil_r0;
        t0_equil_r0Lil_all(ss) = t0_equil_r0Lil;
        t0_equil_r0Lil_Lilctrl_all(ss) = t0_equil_r0Lil_Lilctrl;
    end
    
    %Equilibration time-scales
    Vmax_tau_equil(ss) = Vmax_tau_equil_sim;
    rmax_tau_equil(ss) = rmax_tau_equil_sim;
    rrad_tau_equil(ss) = rrad_tau_equil_sim;
    r0_tau_equil(ss) = r0_tau_equil_sim;
    r0Lil_tau_equil(ss) = r0Lil_tau_equil_sim;
    Vmax_tau_equil_g(ss) = Vmax_tau_equil_g_sim;
    rmax_tau_equil_g(ss) = rmax_tau_equil_g_sim;
    rrad_tau_equil_g(ss) = rrad_tau_equil_g_sim;
    r0_tau_equil_g(ss) = r0_tau_equil_g_sim;
    r0ER11_tau_equil_g(ss) = r0ER11_tau_equil_g_sim;
    r0Lil_tau_equil_g(ss) = r0Lil_tau_equil_g_sim;
    r0Lil_Lilctrl_tau_equil_g(ss) = r0Lil_Lilctrl_tau_equil_g_sim;
    
    %% Transient data %%%%%
    tau_gen(ss) = tau_gen_sim; %defined using the GRADIENT wind
%    Vmax_gen(ss) = Vmax_gen_sim;
%    rmax_gen(ss) = rmax_gen_sim;
%    rrad_gen(ss) = rrad_gen_sim;
%    r0_gen(ss) = r0_gen_sim;
%    r0Lil_gen(ss) = r0Lil_gen_sim;
    Vmax_gen_g(ss) = Vmax_gen_g_sim;
    rmax_gen_g(ss) = rmax_gen_g_sim;
    rrad_gen_g(ss) = rrad_gen_g_sim;
    r0_gen_g(ss) = r0_gen_g_sim;
    r0ER11_gen_g(ss) = r0ER11_gen_g_sim;
    r0Lil_gen_g(ss) = r0Lil_gen_g_sim;
    r0Lil_Lilctrl_gen_g(ss) = r0Lil_Lilctrl_gen_g_sim;
    
    Vmax_tau_max_g(ss) = Vmax_tau_max_g_sim; %defined using the GRADIENT wind
    Vmax_max_g(ss) = Vmax_max_g_sim;
    rmax_tau_max_g(ss) = rmax_tau_max_g_sim; %defined using the GRADIENT wind
    rmax_max_g(ss) = rmax_max_g_sim;
    rrad_tau_max_g(ss) = rrad_tau_max_g_sim; %defined using the GRADIENT wind
    rrad_max_g(ss) = rrad_max_g_sim;
    r0_tau_max_g(ss) = r0_tau_max_g_sim; %defined using the GRADIENT wind
    r0_max_g(ss) = r0_max_g_sim;
    r0ER11_tau_max_g(ss) = r0ER11_tau_max_g_sim; %defined using the GRADIENT wind
    r0ER11_max_g(ss) = r0ER11_max_g_sim;
    r0Lil_tau_max_g(ss) = r0Lil_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_max_g(ss) = r0Lil_max_g_sim;
    r0Lil_Lilctrl_tau_max_g(ss) = r0Lil_Lilctrl_tau_max_g_sim; %defined using the GRADIENT wind
    r0Lil_Lilctrl_max_g(ss) = r0Lil_Lilctrl_max_g_sim;
%    Vmax_tau_max(ss) = Vmax_tau_max_sim; %defined using the GRADIENT wind
%    Vmax_max(ss) = Vmax_max_sim;
%    rmax_tau_max(ss) = rmax_tau_max_sim; %defined using the GRADIENT wind
%    rmax_max(ss) = rmax_max_sim;
%    rrad_tau_max(ss) = rrad_tau_max_sim; %defined using the GRADIENT wind
%    rrad_max(ss) = rrad_max_sim;
%    r0_tau_max(ss) = r0_tau_max_sim; %defined using the GRADIENT wind
%    r0_max(ss) = r0_max_sim;
%    r0Lil_tau_max(ss) = r0Lil_tau_max_sim; %defined using the GRADIENT wind
%    r0Lil_max(ss) = r0Lil_max_sim;
    
    %Variability during equilibration period about those values
    if(equil_dynamic==1)
        Vmax_equil_g_ts_max(ss) = max(Vmax_equil_g_ts);
        Vmax_equil_g_ts_min(ss) = min(Vmax_equil_g_ts);
        Vmax_equil_g_ts_var(ss) = nanvar(Vmax_equil_g_ts);
        rmax_equil_g_ts_max(ss) = max(rmax_equil_g_ts);
        rmax_equil_g_ts_min(ss) = min(rmax_equil_g_ts);
        rmax_equil_g_ts_var(ss) = nanvar(rmax_equil_g_ts);
        rrad_equil_g_ts_max(ss) = max(rrad_equil_g_ts);
        rrad_equil_g_ts_min(ss) = min(rrad_equil_g_ts);
        rrad_equil_g_ts_var(ss) = nanvar(rrad_equil_g_ts);
        r0_equil_g_ts_max(ss) = max(r0_equil_g_ts);
        r0_equil_g_ts_min(ss) = min(r0_equil_g_ts);
        r0_equil_g_ts_var(ss) = nanvar(r0_equil_g_ts);
        r0Lil_equil_g_ts_max(ss) = max(r0Lil_equil_g_ts);
        r0Lil_equil_g_ts_min(ss) = min(r0Lil_equil_g_ts);
        r0Lil_equil_g_ts_var(ss) = nanvar(r0Lil_equil_g_ts);
        r0Lil_Lilctrl_equil_g_ts_max(ss) = max(r0Lil_Lilctrl_equil_g_ts);
        r0Lil_Lilctrl_equil_g_ts_min(ss) = min(r0Lil_Lilctrl_equil_g_ts);
        r0Lil_Lilctrl_equil_g_ts_var(ss) = nanvar(r0Lil_Lilctrl_equil_g_ts);
    end    

    %% Save data for all simulations
%    Vmax_movave_all(:,ss)=Vmax_movave_sim;
%    rmax_movave_all(:,ss)=rmax_movave_sim;
%    rrad_movave_all(:,ss)=rrad_movave_sim;
%    r0_movave_all(:,ss)=r0_movave_sim;
%    r0Lil_movave_all(:,ss)=r0Lil_movave_sim;

    Vmax_movave_g_all(:,ss)=Vmax_movave_g_sim;
    rmax_movave_g_all(:,ss)=rmax_movave_g_sim;
    rrad_movave_g_all(:,ss)=rrad_movave_g_sim;
    r0_movave_g_all(:,ss)=r0_movave_g_sim;
    r0ER11_movave_g_all(:,ss)=r0ER11_movave_g_sim;
    r0Lil_movave_g_all(:,ss)=r0Lil_movave_g_sim;
    r0Lil_Lilctrl_movave_g_all(:,ss)=r0Lil_Lilctrl_movave_g_sim;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{ss} = xvals_sub_sim;
%    data_tmean_usr_all{ss} = data_tmean_usr_sim;
    data_tmean_usr_g_all{ss} = data_tmean_usr_g_sim;

    %% Extract and keep mpi, fcor, Cd, L_R%%%%%
    mpi_all(ss) = mpi;
    fcor_all(ss) = fcor;
    Cd_all(ss) = Cd_in;
    L_R_all(ss) = L_R;

    %% Extract H, tropospheric depth from RCE sounding
    Ttpp = str2num(subdir(strfind(subdir,'Tthresh')+7:strfind(subdir,'Tthresh')+9));
    if(isempty(Ttpp))
        Ttpp = 200;
    end
    Htpp = 1000*zz00(find(T00<Ttpp,1)-1)/1000+dz*(T00(find(T00<Ttpp,1)-1)-Ttpp)/(T00(find(T00<Ttpp,1)-1)-T00(find(T00<Ttpp,1)));

    %% Extract mass-weighted N = sqrt((g/th0)*(dtheta/dz)) from RCE sounding
    th_trop = th00(zz00<=Htpp);
    zz_trop = zz00(zz00<=Htpp);
    pp_trop = pp00(zz00<=Htpp);
    dthdz_trop = (th_trop(2:end)-th_trop(1:end-1))./(zz_trop(2:end)-zz_trop(1:end-1));
    dp_trop = pp_trop(2:end)-pp_trop(1:end-1);  %mass of each layer
    dthdz_mean = sum((dthdz_trop.*dp_trop))/sum(dp_trop);
    
end



%% PLOTTING %%%%%%%%%%%%%%%%

if(numday_smooth == 0)
    numpt_smooth = 1;
else
    numpt_smooth = numday_smooth*86400/dt;
end

set(0,'defaultaxesfontsize',8,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%% RESCALE TIME
t_day_normal = t_day;
for i=1:numruns
    t_day(i,:) = t_day_normal*time_rescale(i);
end

%Multi simulation: plot time-series of each variable
%{
if(plot_full == 1)
    %% FULL WIND %%%%%%%%%%%
    %Plot time series
    figure(14)
    
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
        plot(t_day(i,:),Vmax_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    %    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 Vmax_pl])
    
    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),Vmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max(i)))
            h=plot(t_day(i,round(Vmax_tau_max(i)/(dt/60/60/24))),Vmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil(i)))
            h=plot(t_day(i,round(Vmax_tau_equil(i)/(dt/60/60/24))),Vmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day(i,:),rmax_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 rmmax_pl])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),rmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max(i)))
            h=plot(t_day(i,round(rmax_tau_max(i)/(dt/60/60/24))),rmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil(i)))
            h=plot(t_day(i,round(rmax_tau_equil(i)/(dt/60/60/24))),rmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %{
%%WHOOPS NOT GOOD TO USE IN A DIMENSIONAL PLOT
    %rrad
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day(i,:),rrad_movave_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(rrad_tau_equil(i)))
            h=plot(t_day(i,round(rrad_tau_equil(i)/(dt/60/60/24))),rrad_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    input_title=sprintf('rrad [km]');
    title(input_title)
    axis([0 tf 0 1.1*max(max(rrad_movave_all(100:end,:)))])
    %}
    %{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day(i,:),r0_movave_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil(i)))
            h=plot(t_day(i,round(r0_tau_equil(i)/(dt/60/60/24))),r0_equil(i),'s');
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
        plot(t_day(i,:),r0Lil_movave_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 tf 0 r0max_pl])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),r0Lil_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max(i)))
            h=plot(t_day(i,round(r0Lil_tau_max(i)/(dt/60/60/24))),r0Lil_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(i,round(r0Lil_tau_equil(i)/(dt/60/60/24))),r0Lil_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    hold on
    input_title1=sprintf('BLUE TIME HALVED FULL storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    
end
%}
if(plot_gradient == 1)
    %% GRADIENT WIND %%%%%%%%
    
    %%Plot time series
    figure(15)
    hold off
    
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
        plot(t_day(i,:),smooth(Vmax_movave_g_all(:,i),numpt_smooth),pl_clrs{i})
        hold on
        
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 Vmax_pl])
    
    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),Vmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max_g(i)))
            h=plot(t_day(i,round(Vmax_tau_max_g(i)/(dt/60/60/24))),Vmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(i,round(Vmax_tau_equil_g(i)/(dt/60/60/24))),Vmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        rmax_movave_g_all(rmax_movave_g_all(:,i)>rmax_max_g(i),i) = NaN;    %get rid of very high values at beginning
        subplot(3,1,2)
        plot(t_day(i,:),smooth(rmax_movave_g_all(:,i),numpt_smooth),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 rmmax_pl])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),rmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max_g(i)))
            h=plot(t_day(i,round(rmax_tau_max_g(i)/(dt/60/60/24))),rmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil_g(i)))
            h=plot(t_day(i,round(rmax_tau_equil_g(i)/(dt/60/60/24))),rmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %{
%%WHOOPS NOT GOOD TO USE IN A DIMENSIONAL PLOT
    %rrad
    for i=1:numruns
        
        subplot(4,1,3)
        plot(t_day(i,:),rrad_movave_g_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rrad [km]');
    title(input_title)
    ylabel(input_title);
    axis([0 tf 0 1.1*max(max(rrad_movave_g_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),rrad_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rrad_tau_max_g(i)))
            h=plot(t_day(i,round(rrad_tau_max_g(i)/(dt/60/60/24))),rrad_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rrad_tau_equil(i)))
            h=plot(t_day(i,round(rrad_tau_equil_g(i)/(dt/60/60/24))),rrad_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    %}
    %r0Lil_Lilctrl
%{
    for i=1:numruns
        
        subplot(4,1,3)
        plot(t_day(i,:),r0Lil_Lilctrl_movave_g_all(:,i),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil_{ctrl} [km]');
    title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 r0max_pl])
%}    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),r0Lil_Lilctrl_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_Lilctrl_tau_max_g(i)))
            h=plot(t_day(i,round(r0Lil_Lilctrl_tau_max_g(i)/(dt/60/60/24))),r0Lil_Lilctrl_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_Lilctrl_tau_equil_g(i)))
            h=plot(t_day(i,round(r0Lil_Lilctrl_tau_equil_g(i)/(dt/60/60/24))),r0Lil_Lilctrl_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %r0Lil
    for i=1:numruns
        r0Lil_movave_g_all(r0Lil_movave_g_all(:,i)>r0Lil_max_g(i),i) = NaN;    %get rid of very high values at beginning
        subplot(3,1,3)
        plot(t_day(i,:),smooth(r0Lil_movave_g_all(:,i),numpt_smooth),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    legend(legend_info,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 tf 0 r0max_pl])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),r0Lil_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max_g(i)))
            h=plot(t_day(i,round(r0Lil_tau_max_g(i)/(dt/60/60/24))),r0Lil_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil_g(i)))
            h=plot(t_day(i,round(r0Lil_tau_equil_g(i)/(dt/60/60/24))),r0Lil_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    hold on
    input_title1=sprintf('GRADIENT storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    
    runnames=subdirs_set;
    text1=text(-4.75,-4.75,runnames,'FontSize',30);
    set(text1,'HorizontalAlignment','left','VerticalAlignment','bottom','Interpreter','Latex','BackgroundColor','white','EdgeColor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot non-dimensional time series %%%%%%%%%%%%%%%%
    figure(151)
    hold off
    
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
        plot(t_day(i,:),smooth(Vmax_movave_g_all(:,i)/mpi_all(i),numpt_smooth),pl_clrs{i})
        hold on
        
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('V_m / V_p');
    %title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 Vmax_nd_pl])
    
    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),Vmax_gen_g(i)/mpi_all(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_max_g(i)))
            h=plot(t_day(i,round(Vmax_tau_max_g(i)/(dt/60/60/24))),Vmax_max_g(i)/mpi_all(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(i,round(Vmax_tau_equil_g(i)/(dt/60/60/24))),Vmax_equil_g(i)/mpi_all(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day(i,:),smooth(rmax_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),numpt_smooth),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r_m / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 rmmax_nd_pl])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),rmax_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_max_g(i)))
            h=plot(t_day(i,round(rmax_tau_max_g(i)/(dt/60/60/24))),rmax_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rmax_tau_equil_g(i)))
            h=plot(t_day(i,round(rmax_tau_equil_g(i)/(dt/60/60/24))),rmax_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    %{
%%WHOOPS NOT GOOD TO USE IN A DIMENSIONAL PLOT
    %rrad
    for i=1:numruns
        
        subplot(4,1,3)
        plot(t_day(i,:),rrad_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r_{rad} / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
    %xlabel('time [days]')
    axis([0 tf 0 1.1*max(max(rrad_movave_g_all(20:end,:)/(mpi_all(i)/fcor_all(i)/1000)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),rrad_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rrad_tau_max_g(i)))
            h=plot(t_day(i,round(rrad_tau_max_g(i)/(dt/60/60/24))),rrad_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(rrad_tau_equil_g(i)))
            h=plot(t_day(i,round(rrad_tau_equil_g(i)/(dt/60/60/24))),rrad_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    %}
    %r0Lil_Lilctrl
%{
    for i=1:numruns
        
        subplot(4,1,3)
        plot(t_day(i,:),r0Lil_Lilctrl_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r_{0 ctrl} / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
    %    xlabel('time [days]')
    axis([0 tf 0 r0max_nd_pl])
%}    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),r0Lil_Lilctrl_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_Lilctrl_tau_max_g(i)))
            h=plot(t_day(i,round(r0Lil_Lilctrl_tau_max_g(i)/(dt/60/60/24))),r0Lil_Lilctrl_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_Lilctrl_tau_equil_g(i)))
            h=plot(t_day(i,round(r0Lil_Lilctrl_tau_equil_g(i)/(dt/60/60/24))),r0Lil_Lilctrl_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day(i,:),smooth(r0Lil_movave_g_all(:,i)/(mpi_all(i)/fcor_all(i)/1000),numpt_smooth),pl_clrs{i})
        hold on
    end
    grid on
%    legend(input_legend,'Location','EastOutside')
    legend(legend_info,'Location','EastOutside')
    input_title=sprintf('r_0 / (V_p / f)');
    %title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    %axis([0 tf 0 1.1*max(max(r0Lil_movave_g_all(20:end,:)/(mpi_all(i)/fcor_all(i)/1000)))])
    axis([0 tf 0 1.5])
    
    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(i,round(tau_gen(i)/(dt/60/60/24))),r0Lil_gen_g(i)/(mpi_all(i)/fcor_all(i)/1000),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_max_g(i)))
            h=plot(t_day(i,round(r0Lil_tau_max_g(i)/(dt/60/60/24))),r0Lil_max_g(i)/(mpi_all(i)/fcor_all(i)/1000),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
        if(~isnan(r0Lil_tau_equil_g(i)))
            h=plot(t_day(i,round(r0Lil_tau_equil_g(i)/(dt/60/60/24))),r0Lil_equil_g(i)/(mpi_all(i)/fcor_all(i)/1000),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i}(1));
        end
    end
    
    hold on
    input_title1=sprintf('TIME-MOD GRADIENT non-dimensional storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    
    
end

