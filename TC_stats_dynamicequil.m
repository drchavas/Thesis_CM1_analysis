%TC_stats_dynamicequil.m

%Created: 04 June 12, Dan Chavas

%Purpose: This file copies the file output by TC_stats, opens it,
%recalculates the equilibrium values for each structural variable,
%and then resaves the file in a new sub-directory.
%The equilibrium is calculated as the most stable 30 day period anytime
%after day 70.

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

function [junk] = TC_stats_dynamicequil(tf,T_mean,dt_final,dt_final_dynamic,subdir,dir_home);
junk='junk';
%clear
%clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
%T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
%dt_final = 50;  %[days]; original length of period over which equilibrium is calculated
%tf = 150;   %[days]; original end day of equilibrium calculation
%dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad
save_file_dynamicequil = 1;

stop_overwrite = 0;
if(save_file_dynamicequil==1)
%check if new subdirectory exists
if(wrad_const == 1)
    subdir_dyn = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
else
    subdir_dyn = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
end

if(exist(subdir_dyn)==7)

    %check to make sure you're not overwriting original TC_stats output
    if(exist(sprintf('%s/NOTANORIGINALTCSTATS.mat',subdir_dyn))==2)   %TRUE if directory exists but does NOT contain this file
    else
        sprintf('You are trying to overwrite original TC_stats output!  Bad.')
        stop_overwrite = 1;
    end

    
else
    
    if(wrad_const == 1)
        mkdir(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simplots_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst/PLOTS',T_mean,dt_final_dynamic))

        cd(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic))
        save NOTANORIGINALTCSTATS.mat T_mean tf dt_final dt_final_dynamic
        cd(dir_home)
    else
        
        mkdir(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simplots_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
        mkdir(sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic/PLOTS',T_mean,dt_final_dynamic))

        cd(sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic))
        save NOTANORIGINALTCSTATS.mat T_mean tf dt_final dt_final_dynamic
        cd(dir_home)
    end
    
end
end

if(stop_overwrite == 0)

    if(wrad_const == 1)
        file_in_dyn = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst/ax%s.mat',T_mean,tf-dt_final,tf,subdir);
    else
        file_in_dyn = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir);
    end

    
    if(exist(file_in_dyn)==2)
        sprintf('Updating equilibrium values for ax%s',subdir)

        %%Load data for given simulation to fix things within if needed
        load(file_in_dyn);

        %%Find most stable [dt_final_dynamic]-day period after day 60 in the timeseries for Vmax
        t_start = 60;   %[days]
        i_start = t_start/(dt/60/60/24) + 1;    %index of day [t_start]
        t_end = 150;    %[days]
        i_end = t_end/(dt/60/60/24) + 1;    %index of day [t_end]
        ni_final_dynamic = dt_final_dynamic/(dt/60/60/24)+1;    %number of points needed to cover period
        
        %%Extract data for [dt_final_dynamic]-day period at end of simulation and calculate its mean
        for ii = 1:((i_end-ni_final_dynamic)-i_start+2)
            %Note: ii corresponds to index for FIRST timestep in period (relative to i_start)
            
            indices = i_start+ii-1:i_start+ii-1+(ni_final_dynamic-1);
            
            t_movave_temp = t_day(indices);
            Vmax_movave_g_temp = Vmax_movave_g(indices);
            clear indices
            Vmax_equil_g_temp(ii) = nanmean(Vmax_movave_g_temp);
            Vmax_variance_g_temp(ii) = nanvar(Vmax_movave_g_temp);

        end
        
        %%Determine timeperiod of minimum variance and set equilibrium value accordingly
        i_start_equil = find(Vmax_variance_g_temp==min(Vmax_variance_g_temp));
        indices_equil = i_start+i_start_equil-1:i_start+i_start_equil-1+(ni_final_dynamic-1);
        Vmax_equil_g_dynamic = Vmax_equil_g_temp(i_start_equil)
        t0_equil = t_day(indices_equil(1))
        tf_equil = t_day(indices_equil(end))
        
        numvars = 10;     %Vmax, rmax, rmid, r0, r0Lil, Vmax_g, rmax_g, rmid_g, r0_g, r0Lil_g
        for l=6:numvars     %SKIP THE FULL WIND ONES FOR NOW

            switch l
                case 1
                    var_movave = Vmax_movave;
                case 2
                    var_movave = rmax_movave;
                case 3
                    var_movave = rmid_movave;
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
                    var_movave = rmid_movave_g;
                    %sprintf('rmid_g')
                case 9
                    var_movave = r0_movave_g;
                    %sprintf('r0_g')
                case 10
                    var_movave = r0Lil_movave_g;
                    %sprintf('r0Lil_g')
            end

            %%Extract data for (dt_final)-day period at end of simulation and calculate its mean
            var_dat_end = var_movave(indices_equil);
            var_mean_end = nanmean(var_dat_end);
            var_equil(l) = var_mean_end;

        end
%}

        %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
        %variable values
        Vmax_equil_g_sim = var_equil(6);
        rmax_equil_g_sim = var_equil(7);
        rmid_equil_g_sim = var_equil(8);
        r0_equil_g_sim = var_equil(9);
        r0Lil_equil_g_sim = var_equil(10);

        %% Save data to file
        if(save_file_dynamicequil == 1)
            save tempstuff.mat l ii
            clear l i ii
            save temp.mat
            load tempstuff.mat
            
            if(wrad_const == 1)
                movefile('temp.mat',sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst/ax%s.mat',T_mean,dt_final_dynamic,subdir))
            else
                movefile('temp.mat',sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic/ax%s.mat',T_mean,dt_final_dynamic,subdir))
            end

            
            delete('tempstuff.mat')
        end
        
    end


end

end
