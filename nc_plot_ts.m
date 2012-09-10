%nc_plot_ts.m

%Created: 30 Jan 2011, Dan Chavas

%Purpose: to plot timeseries of data from multiple netcdf files

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!


clear
clc
figure(1)
clf(1)     %clear figures without closing them

%% USER INPUT %%%%%%%%%%%%%%%%%%

subdir_pre='CTRL_icRCE/';
subdir_pre='TRANSFER/';
%subdir_pre='RCE/';    %general subdir that includes multiple runs within
%subdir_pre='';
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types= [3 1 1 1]; %1=axisym; 3=3D
%    toax = 1.688/1.326;
%    div_dat_by=[toax toax 1 4 4*toax 16 16*toax 64 64*toax 4*64*toax 1 1 1 1]; %for any domain-integrated quantities, can normalize to equivalent domain size
%    div_dat_by=[1 4 4^2 4^3 4^4 4^5 4^6]
%    div_dat_by=[1 1 4];
    div_dat_by=ones(30,1);
subdirs = {

'RCE_test'

}; %name of sub-directory with nc files

t0 = 0;    %[day], starting time for timeseries plot
tf = 150;   %[day], ending time for timeseries plot

numruns=length(subdirs);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)

plot_type = 1;  %0=no plot
                %1=single plot of time series of [var] for above runs
                %4=4-panel:
                    %axisym: vinterp, w, qv, uinterp
                    %3D: vinterp, mid-level w, qv, XZ-xsec qv
%IF plot_type=1: input variable and domain
%    pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
%    pl_clrs={'b--' 'k' 'r--' 'b' 'g--' 'r' 'c--' 'g' 'y--' 'c' 'y--' 'y'};
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
%    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    %pl_clrs={'b' 'k' 'r' 'g' 'c' 'y' 'm' 'b--' 'k--' 'r--' 'g--' 'c--' 'y--' 'm--'};
    use_stats = 1;  %0=extract from original data; 1=use cm1out_stats.nc
    var = 'svmax';  %variable of interest    
    v_scale = 1; %scale the y-axis variable by this value (note: units will be wrong)
    use_smooth = 1; %0=raw data; 1=smoothed (N-pt moving average)
        numpt_sm = 50;  %smoother uses this number of points
        %IF use_stats=0:
        ts_type=1;  %1=global max; 2=global min; 3=global mean
        x0=0;   %first x grid point [0,end]
        xf=1000;   %first y grid point [0,end]
        y0=0;   %first y grid point [0,end]
        yf=0;   %last y grid point [0,end]
        z0=0;  %first z grid point [0,end]
        zf=20;  %last z grid point [0,end]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Append ax or 3d to subdir names for plot
for i=1:length(subdirs)
    if(run_types(i)==1) %ax
        subdirs_plot{i}=sprintf('ax%s',subdirs{i})
    else    %3d
        subdirs_plot{i}=sprintf('3d%s',subdirs{i})
    end
end

switch plot_type
    case 0  %no plot

    case 1  %single plot of [var] over specified domain
        


        %% LOOP OVER RUNS TO INCUDE IN FIGURE
        for rr=1:numruns

            clear tt
            %%EXTRACT RUN_TYPE AND SUBDIR NAME
            run_type=run_types(rr);
            subdir=subdirs{rr};

            %% OPTIONS FOR EITHER AXISYM OR 3D RUNS
            if(run_type==1)
                run_type_str='axisym';
                if(ext_hd==0)
                    dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/axisym/CM1_output/%s',subdir_pre);
                else    %external harddrive
                    dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/axisym/%s',subdir_pre);
                end
                copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
            elseif(run_type==3)
                run_type_str='3d';
                if(ext_hd==0)
                    dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3d/CM1_output/%s',subdir_pre);
                else    %external harddrive
                    dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/3d/%s',subdir_pre);
                end
                copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
            end
            
            
            %% EXTRACT DATA EITHER USING cm1out_stats.nc OR ORIGINAL nc DATA
dir_in
subdir
            if(use_stats == 1)  %USE cm1out_stats.nc
                [data_ts v_def v_units tt t_units] = nc_extract_stats(dir_in,subdir,var);
                dt = tt(2)-tt(1);   %timestep in stats file
            else    %USE ORIGINAL DATA
                %%GET NUMBER OF FILES
                subdir_full=sprintf('%s%s',dir_in,subdir)
                numfiles=length(dir(sprintf('%s/cm1out_*.nc',subdir_full)))-1;  %ignore cm1out_stats.nc    %includes 6 extra: ".","..",input_sounding,namelist.input,filecontent.info,cm1out_stats.nc

                %% LOOP OVER FILES TO OBTAIN TIME SERIES DATA
                clear data_ts
                for ii=1:numfiles

                    t_file=ii;

                    %%NC FILENAME (default)
                    if(run_type==1)
                        if(t_file<10)
                            nc_file = sprintf('cm1out_000%i.nc',t_file);
                        elseif(t_file<100)
                            nc_file = sprintf('cm1out_00%i.nc',t_file);
                        else
                            nc_file = sprintf('cm1out_0%i.nc',t_file);
                        end
                    elseif(run_type==3)
                        if(t_file<10)
                            nc_file = sprintf('cm1out_t000%i.nc',t_file);
                        elseif(t_file<100)
                            nc_file = sprintf('cm1out_t00%i.nc',t_file);
                        else
                            nc_file = sprintf('cm1out_t0%i.nc',t_file);
                        end    
                    end

                    %%EXTRACT DATA
                    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
                    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

                    %%SAVE DESIRED DATA
                    tt(ii)=time;    %time

                    switch ts_type
                        case 1  %max value in domain
                            data_ts(ii)=max(max(data));
                            data_ts_str=sprintf('Global max %s',var);
                        case 2  %min value in domain
                            data_ts(ii)=min(min(data));
                            data_ts_str=sprintf('Global min %s',var);
                        case 3  %domain average
                            data_ts(ii)=mean(mean(data));
                            data_ts_str=sprintf('Global mean %s',var);
                    end

                end

            end

            %%TIME IN DAYS
            tt_day = tt ./ 86400;   %for all data

            %%PLOT TIMESERIES
            %PLOT DATA
            figure(1)
            hold on
            
            %%DETERMINE INDICES FOR DESIRED INITIAL AND FINAL TIMES
            dt = tt(2)-tt(1);   %timestep in stats file
            i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for plot
            i_tf = min(length(tt),round(tf*24*60*60/dt)+1); %timestep corresponding to end time for plot

            mean(data_ts(i_t0:i_tf))
            if(use_smooth==1)
                plot(tt_day(i_t0:i_tf)',smooth(v_scale.*data_ts(i_t0:i_tf),numpt_sm)'/div_dat_by(rr),pl_clrs{rr},'LineWidth',1)
%                plot(tt_day',smooth(v_scale.*data_ts,numpt_sm)',pl_clrs{rr},'LineWidth',2)
                smooth_str=sprintf('%i-pt smooth',numpt_sm);
            else
                plot(tt_day(i_t0:i_tf)',v_scale.*data_ts(i_t0:i_tf)'/div_dat_by(rr),pl_clrs{rr},'LineWidth',1)
                smooth_str='raw';
            end
%            dat_eq=mean(data_ts(end-30:end))
%    dan(rr)=data_ts(1)/div_dat_by(rr);
            
        end

        if(exist('data_ts_str')~=1) %needed if using cm1out_stats.nc
            data_ts_str=var;
        end
        
        %STUFF FOR PLOT
        set(gca,'fontweight','bold','fontsize',11)
        input_xlabel=sprintf('time [day]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('%s [%s]',v_def,v_units);
        ylabel(input_ylabel);
        input_title=sprintf('Time series: %s [%s] (%s)',data_ts_str,v_units,smooth_str);
        title(input_title)
        input_legend=strrep(subdirs_plot(1:numruns),'_','\_');
        legend(input_legend)
        



    otherwise
        sprintf('Invalid plot number. Fail.')
end

%}