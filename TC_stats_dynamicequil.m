%TC_stats_dynamicequil.m

%Created: 04 June 12, Dan Chavas

%Purpose: This file copies the file output by TC_stats, opens it,
%recalculates the equilibrium values for each structural variable,
%and then resaves the file in a new sub-directory.
%The equilibrium is calculated as the most stable 30 day period anytime
%after day 70.

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

function [junk] = TC_stats_dynamicequil(tf,T_mean,dt_final,dt_final_dynamic,subdir,dir_home,dir_in_dat,dir_in_dat_dyn,t_start_dynequil);

%%Write out to screen whats going on
sprintf('TC_stats_dynamicequil for: %s',subdir)

junk='junk';
%clear
%clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
%T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
%dt_final = 50;  %[days]; original length of period over which equilibrium is calculated
%tf = 150;   %[days]; original end day of equilibrium calculation
%dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
save_file_dynamicequil = 1;

stop_overwrite = 0;
if(save_file_dynamicequil==1)
%check if new subdirectory exists
subdir_dyn = sprintf('%s',dir_in_dat_dyn);

%if(exist(subdir_dyn)==7)
assert(exist(subdir_dyn)==7)    %don't let the program go farther if this sub-directory doesn't already exist

%check to make sure you're not overwriting original TC_stats output
if(exist(sprintf('%s/NOTANORIGINALTCSTATS.mat',subdir_dyn))==2)   %TRUE if directory exists but does NOT contain this file
else
    sprintf('You are trying to overwrite original TC_stats output!  Bad.')
    stop_overwrite = 1;
end

    
%else
            
%    mkdir(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
%    mkdir(sprintf('../CM1_postproc_data/simplots_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
%    mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
%    mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic/PLOTS',T_mean,dt_final_dynamic))

%    cd(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
%    save NOTANORIGINALTCSTATS.mat T_mean tf dt_final dt_final_dynamic
%    cd(dir_home)
    
%end
%end

if(stop_overwrite == 0)

    file_in_dyn = sprintf('%s/ax%s.mat',dir_in_dat,subdir);
    
    if(exist(file_in_dyn)==2)
        sprintf('Updating equilibrium values for ax%s',subdir)

        %%Load data for given simulation to fix things within if needed
        load(file_in_dyn);

        %%Find most stable [dt_final_dynamic]-day period after day 60 in the timeseries for Vmax
        i_start = t_start_dynequil/(dt/60/60/24) + 1;    %index of day [t_start_dynequil]
        t_end = min(tf,(numfiles-1)*(dt/60/60/24));    %[days]
        i_end = t_end/(dt/60/60/24) + 1;    %index of day [t_end]
        ni_final_dynamic = dt_final_dynamic/(dt/60/60/24)+1;    %number of points needed to cover period
        
        %%Extract data for [dt_final_dynamic]-day period and calculate its mean
        for ii = 1:((i_end-ni_final_dynamic)-i_start+2)
            %Note: ii corresponds to index for FIRST timestep in period (relative to i_start)
            
            indices = i_start+ii-1:i_start+ii-1+(ni_final_dynamic-1);
            
            Vmax_movave_g_temp = Vmax_movave_g(indices);
            rmax_movave_g_temp = rmax_movave_g(indices);
            rrad_movave_g_temp = rrad_movave_g(indices);
            r0_movave_g_temp = r0_movave_g(indices);
            r0Lil_movave_g_temp = r0Lil_movave_g(indices);
            r0Lil_Lilctrl_movave_g_temp = r0Lil_Lilctrl_movave_g(indices);
            clear indices

            %%Save all [dt_equil]-day running mean values
            Vmax_equil_g_ts(ii) = nanmean(Vmax_movave_g_temp); 
            rmax_equil_g_ts(ii) = nanmean(rmax_movave_g_temp); 
            rrad_equil_g_ts(ii) = nanmean(rrad_movave_g_temp); 
            r0_equil_g_ts(ii) = nanmean(r0_movave_g_temp); 
            r0Lil_equil_g_ts(ii) = nanmean(r0Lil_movave_g_temp); 
            r0Lil_Lilctrl_equil_g_ts(ii) = nanmean(r0Lil_Lilctrl_movave_g_temp); 

            %%Calculate and save timeseries of variances in Vmax
            Vmax_variance_g_temp(ii) = nanvar(Vmax_movave_g_temp);
            rmax_variance_g_temp(ii) = nanvar(rmax_movave_g_temp);
            rrad_variance_g_temp(ii) = nanvar(rrad_movave_g_temp);
            r0_variance_g_temp(ii) = nanvar(r0_movave_g_temp);
            r0Lil_variance_g_temp(ii) = nanvar(r0Lil_movave_g_temp);
            r0Lil_Lilctrl_variance_g_temp(ii) = nanvar(r0Lil_Lilctrl_movave_g_temp);

        end
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(Vmax_variance_g_temp==min(Vmax_variance_g_temp),1);
        indices_equil_Vmax = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        Vmax_equil_g_dynamic = Vmax_equil_g_ts(i_start_equil);
        t0_equil_Vmax = t_day(indices_equil_Vmax(1));
        tf_equil_Vmax = t_day(indices_equil_Vmax(end));
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(rmax_variance_g_temp==min(rmax_variance_g_temp),1);
        indices_equil_rmax = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        rmax_equil_g_dynamic = rmax_equil_g_ts(i_start_equil);
        t0_equil_rmax = t_day(indices_equil_rmax(1));
        tf_equil_rmax = t_day(indices_equil_rmax(end));

        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(rrad_variance_g_temp==min(rrad_variance_g_temp),1);
        indices_equil_rrad = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        rrad_equil_g_dynamic = rrad_equil_g_ts(i_start_equil);
        t0_equil_rrad = t_day(indices_equil_rrad(1));
        tf_equil_rrad = t_day(indices_equil_rrad(end));
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(r0_variance_g_temp==min(r0_variance_g_temp),1);
        indices_equil_r0 = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        r0_equil_g_dynamic = r0_equil_g_ts(i_start_equil);
        t0_equil_r0 = t_day(indices_equil_r0(1));
        tf_equil_r0 = t_day(indices_equil_r0(end));
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(r0Lil_variance_g_temp==min(r0Lil_variance_g_temp),1);
        indices_equil_r0Lil = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        r0Lil_equil_g_dynamic = r0Lil_equil_g_ts(i_start_equil);
        t0_equil_r0Lil = t_day(indices_equil_r0Lil(1));
        tf_equil_r0Lil = t_day(indices_equil_r0Lil(end));
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(r0Lil_Lilctrl_variance_g_temp==min(r0Lil_Lilctrl_variance_g_temp),1);
        indices_equil_r0Lil_Lilctrl = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        r0Lil_Lilctrl_equil_g_dynamic = r0Lil_Lilctrl_equil_g_ts(i_start_equil);
        t0_equil_r0Lil_Lilctrl = t_day(indices_equil_r0Lil_Lilctrl(1));
        tf_equil_r0Lil_Lilctrl = t_day(indices_equil_r0Lil_Lilctrl(end));
        
        %Used for the equilibration period for the radial wind profile
        t0_equil = t0_equil_Vmax;
        tf_equil = tf_equil_Vmax;
        indices_equil = indices_equil_Vmax;
        
        %%Calculate mean radial wind profile during equilibrium period
        data_tmean_usr_g_sim = mean(v_r_all(:,indices_equil),2);
        
        numvars = 11;     %Vmax, rmax, rrad, r0, r0Lil, Vmax_g, rmax_g, rrad_g, r0_g, r0Lil_g
        for l=6:numvars     %SKIP THE FULL WIND ONES FOR NOW

            switch l
                case 1
                    var_movave = Vmax_movave;
                case 2
                    var_movave = rmax_movave;
                case 3
                    var_movave = rrad_movave;
                case 4
                    var_movave = r0_movave;
                case 5
                    var_movave = r0Lil_movave;
                case 6
                    var_movave = Vmax_movave_g;
                    %sprintf('Vmax_g')
                case 7
                    var_movave = rmax_movave_g;
                    %sprintf('rmax_g')
                case 8
                    var_movave = rrad_movave_g;
                    %sprintf('rrad_g')
                case 9
                    var_movave = r0_movave_g;
                    %sprintf('r0_g')
                case 10
                    var_movave = r0Lil_movave_g;
                    %sprintf('r0Lil_g')
                case 11
                    var_movave = r0Lil_Lilctrl_movave_g;
                    %sprintf('r0Lil_Lilctrl_g')
                case 12
                    var_movave = r0ER11_movave_g;
                    %sprintf('r0ER11_g')
            end

            %%Extract data for (dt_final)-day period at end of simulation and calculate its mean
            var_dat_end = var_movave(indices_equil);
            var_mean_end = nanmean(var_dat_end);
            var_equil(l) = var_mean_end;
            tdat_end = t_day(indices_equil);

            %%Identify timescale to equilibration for dynamic equilibrium
            %%Check slope of linear trend to determine if it is at statistical equilibrium
            p = polyfit(tdat_end-mean(tdat_end),var_dat_end'-mean(var_dat_end)',1);
            slope=p(1); %[ms-1/day]
            slopes(l)=slope;
            clear p

            slope_equil = .01*var_mean_end;    %[1 percent/day]; slopes below this are considered at statistical equilibrium
            slopes_equil(l) = slope_equil;
            var_equil(l) = var_mean_end;
        %        if(abs(slope)<abs(slope_equil))   %IS at statistical equilibrium

            %%Find tau_equil by repeating same procedure forward in time
            %%from the start until first period where slope<slope_c is satisfied 
            %%and mean is <10% of equilibrium value; then repeat recursively
            %%within this subset with an averaging period half the length
            %%until T = T_mean
            disequil_frac = .1;    %define disequilibrium as having mean that deviates from equilibrium value by greater than 10%
            istep = 1;  %steps backwards in time to find disequilibrium point
            t_day_temp = t_day;
            var_movave_temp = var_movave;
            dt_equil_temp = dt_equil;

            while(dt_equil_temp > T_mean)
                tdat_sub = t_day_temp(istep:ceil(dt_equil_temp/(dt/60/60/24))+istep);
                var_dat_sub = var_movave_temp(istep:ceil(dt_equil_temp/(dt/60/60/24))+istep);

                %%DRC 06 Jul 2012 make sure there are sufficient good data!
                nanfrac = (length(var_dat_sub)-sum(isnan(var_dat_sub)))/length(var_dat_sub);    %fraction of datapoints that are NaN
                if(nanfrac>.9)  %if there are sufficient good data, then check equilibration, otherwise move to averaging period up one timestep
                    indices = ~isnan(var_dat_sub);
                    tdat_sub = tdat_sub(indices);
                    var_dat_sub = var_dat_sub(indices);
                    clear indices
                    var_mean_sub = mean(var_dat_sub);

                    %%Check slope of linear trend to determine if it is at statistical equilibrium
                    p = polyfit(tdat_sub-mean(tdat_sub),var_dat_sub'-mean(var_dat_sub)',1);
                    slope=p(1); %[ms-1/day]
                    clear p

                    if((abs(var_mean_sub-var_equil(l))<disequil_frac*var_equil(l) && abs(slope)<abs(slope_equil)) || ceil(dt_equil_temp/(dt/60/60/24))+istep == length(t_day_temp)) %reaches statistical equilibrium
                        t_day_equil0 = tdat_sub(1);
                        t_day_equilf = tdat_sub(end);
                        t_day_temp = tdat_sub;
                        var_movave_temp = var_dat_sub;
                        dt_equil_temp = dt_equil_temp/2;   %halve this distance for each equilibrium-check iteration
                        istep = 0;
                    end
                end

                istep = istep + 1;    %shift up 1 timestep at a time until reach equilibrium
            end

            var_tau_equil(l) = t_day_equilf;

        end
%}

        %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
        %variable values
        Vmax_equil_g_sim = var_equil(6);
        rmax_equil_g_sim = var_equil(7);
        rrad_equil_g_sim = var_equil(8);
        r0_equil_g_sim = var_equil(9);
        r0Lil_equil_g_sim = var_equil(10);
        r0Lil_Lilctrl_equil_g_sim = var_equil(11);
        r0ER11_equil_g_sim = var_equil(12);
        
        %timescales to those values
        Vmax_tau_equil_g_sim = var_tau_equil(6);
        rmax_tau_equil_g_sim = var_tau_equil(7);
        rrad_tau_equil_g_sim = var_tau_equil(8);
        r0_tau_equil_g_sim = var_tau_equil(9);
        r0Lil_tau_equil_g_sim = var_tau_equil(10);
        r0Lil_Lilctrl_tau_equil_g_sim = var_tau_equil(11);
        r0ER11_tau_equil_g_sim = var_tau_equil(12);

        %% Save data to file
        if(save_file_dynamicequil == 1)
            save tempstuff.mat l ii
            clear l i ii
            save temp.mat
            load tempstuff.mat
            
            movefile('temp.mat',sprintf('%s/ax%s.mat',dir_in_dat_dyn,subdir))

            delete('tempstuff.mat')
        end
        
    end


end

end
