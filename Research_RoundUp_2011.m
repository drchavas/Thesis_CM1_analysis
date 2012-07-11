%nc_plot_size.m

%Created: 05-05-11, Dan Chavas

%Purpose: to plot quasi-steady state radial wind profiles from many runs

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc
figure(1)
figure(2)
clf(1)     %clear figures without closing them
clf(2)

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRE87/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

t0 = 100;
tf = 150;

%mult_by = [8 4 2 1 1/2 1/4 1/8];

plot_type = 1;  %0=no plot
                %1=single plot of time series of [var] for above runs
                %4=4-panel:
                    %axisym: vinterp, w, qv, uinterp
                    %3D: vinterp, mid-level w, qv, XZ-xsec qv
%IF plot_type=1: input variable and domain
%    pl_clrs={'b--' 'k' 'r--' 'b' 'g--' 'r' 'c--' 'g' 'y--' 'c' 'y--' 'y'};
    %pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
    pl_clrs={'b' 'r' 'g' 'c' 'y' 'm' 'k--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
%    pl_clrs={'b--' 'k--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'b' 'k' 'r' 'g' 'c' 'y' 'm'};
    %pl_clrs={'g' 'c' 'k' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--'};
%    pl_clrs={'r' 'g' 'c' 'k' 'y' 'm' 'r--' 'g--' 'c--' 'k--' 'y--'};
    var = 'vinterp';  %vinterp; variable of interest    
    x0=0;   %first x grid point [0,end]
    xf=1000;   %first y grid point [0,end]
    y0=0;   %first y grid point [0,end]
    yf=0;   %last y grid point [0,end]
    z0=2;  %first z grid point [0,end]
    zf=2;  %last z grid point [0,end]
    vmin_plot = 0;  %lowest value plotted
    vmax_plot = 50; %highest value plotted
    plot_normalized = 0;    %1=plot V/Vmax, r/rmax; 0=plot regular V(r)
    r_ts = 0; %1=plot time series of r(v_usr) 0=don't
        v_usr = 12; %find radius of this wind speed; 'MAX' for V_max
        numpt_sm = 50;  %smoother uses this number of points
        rmin_plot = 0;  %lowest value plotted
        rmax_plot = 2000;    %highest value plotted

subdirs = {
  
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_lh12000_200days'
'CTRLv0qrhSATqdz5000_div8_dx.25_nx768_short'



}; %name of sub-directory with nc files

    
numruns=length(subdirs);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Append ax or 3d to subdir names for plot
for i=1:length(subdirs)
    if(run_types(i)==1) %ax
        subdirs_plot{i}=sprintf('ax%s',subdirs{i})
    else    %3d
        subdirs_plot{i}=sprintf('3d%s',subdirs{i})
    end
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;
for rr=1:numruns

    if(rr==2)
        t0=65;
        tf=75;
    end
%{
    if(rr==1)
        t0 = 100;    %[day], starting time for averaging
        tf = 150;   %[day], ending time for averaging
    else
        t0 = 50;
        tf = 65;
    end
%}    
    %%EXTRACT RUN_TYPE AND SUBDIR NAME
    run_type=run_types(rr);
    subdir=subdirs{rr};

    %% OPTIONS FOR EITHER AXISYM OR 3D RUNS
    if(run_type==1)
        run_type_str='axisym';
        if(ext_hd==0)
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/axisym/CM1_output/%s',subdir_pre);
        else    %external harddrive
            dir_in=sprintf('/Volumes/CHAVAS_CM1/CM1_output/axisym/%s',subdir_pre);
        end        
        copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
    elseif(run_type==3)
        run_type_str='3D';
        if(ext_hd==0)
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3D/CM1_output/%s',subdir_pre);
        else    %external harddrive
            dir_in=sprintf('/Volumes/CHAVAS_CM1/CM1_output/3D/%s',subdir_pre);
        end
        copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
    end
    
%{
%% Determine if things are reasonably steady
    [data_ts v_def v_units tt t_units] = nc_extract_stats(dir_in,subdir,'vmax');
    
    dt=(tt(end)-tt(end-1))/86400;   %[days/pt]
    numda_smooth=10; %[days]
    numpts_smooth=numda_smooth/dt;
    v_steady = smooth(data_ts,numpts_smooth);
    i_steady=find(v_steady>v_steady(end)-5 || v_steady<v_steady(end)+5);
    i_steady=i_steady(numpts_smooth+1:end);
    
    plot(tt/86400,data_ts,tt/86400,v_steady)
%}

    %%DIRECTORY WITH OUTPUT DATA
    subdir_full=sprintf('%s%s',dir_in,subdir)

    %%EXTRACT TIMESTEP SIZE
    var_dt = 'qvpert'; %doesn't matter, just need any variable
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    if(run_type==3)
        numfiles=length(dir(sprintf('%s/cm1out_t*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_t0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    elseif(run_type==1)
        numfiles=length(dir(sprintf('%s/cm1out_0*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    end
    %xres = dx; %this is done farther down now
    dt = time; %[s]; time of cm1out_t0001.nc is defined as zero
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
    i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

    %%TIME IN DAYS
%    t_min_day = dt*i_t0/86400;
%    t_max_day = dt*i_tf/86400;

    %%initialize the output vectors
    clear v_df_qv v_units_qv rprof_tmean
    rprof_tmean=zeros(xf-x0+1,1);
    
    clear t_day
    clear r_usr
    for ii=1:i_tf-i_t0+1
        fclose('all')
        t_file=i_t0+ii-1;
        
        t_day(ii) = (t_file-1)*dt/(24*60*60);
        
        
        %%NC FILENAME (default)
        if(run_type==1 || run_type==10 || run_type==30)
            if(t_file<10)
                nc_file = sprintf('cm1out_000%i.nc',t_file);
            elseif(t_file<100)
                nc_file = sprintf('cm1out_00%i.nc',t_file);
            elseif(t_file<1000)
                nc_file = sprintf('cm1out_0%i.nc',t_file);
            else
                nc_file = sprintf('cm1out_%i.nc',t_file);
            end
        elseif(run_type==3)
            if(t_file<10)
                nc_file = sprintf('cm1out_t000%i.nc',t_file);
            elseif(t_file<100)
                nc_file = sprintf('cm1out_t00%i.nc',t_file);
            elseif(t_file<1000)
                nc_file = sprintf('cm1out_t0%i.nc',t_file);
            else
                nc_file = sprintf('cm1out_t%i.nc',t_file);
            end    
        end

        %% OPEN FILE AND EXTRACT nx, ny, nz -- only need to do once!
        
        if(ii==1)
            dir_start=pwd;
            dir_tmp = strcat(dir_in,subdir);
            cd(dir_tmp)
            nc_file
            ncid = netcdf.open(nc_file,'WRITE');

            %NX
            if(run_type==3)
                dimid = netcdf.inqDimID(ncid,'nx');
            else
                dimid = netcdf.inqDimID(ncid,'ni');
            end
            [junk, nx] = netcdf.inqDim(ncid,dimid);

            %NY
            if(run_type==3)
                dimid = netcdf.inqDimID(ncid,'ny');
            else
                dimid = netcdf.inqDimID(ncid,'nj');
            end
            [junk, ny] = netcdf.inqDim(ncid,dimid);

            %NZ
            if(run_type==3)
                dimid = netcdf.inqDimID(ncid,'nz');
            else
                dimid = netcdf.inqDimID(ncid,'nk');
            end
            [junk, nz] = netcdf.inqDim(ncid,dimid);

            cd(dir_start)
            
        end

        %%Create mean radial profile for desired variable
        if(strcmp(var,'satdef'))
            %% Calculate various thermodynamic quantities
            snd_file = 'input_sounding';

            %%load initial sounding variables
            dz_frac=.5;
            dir_full = strcat(dir_in,subdir);
            [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc dz] = snd_extract(dir_full,snd_file,dz_frac);
            xres = dx;
            
            %%load prspert
            var_temp = 'prspert';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            pp_pert = squeeze(data);

            %%load thpert
            var_temp = 'thpert';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            th_pert = squeeze(data);

            %%load qvpert
            var_temp = 'qvpert';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qv_pert = squeeze(data);

            %%load upert
            var_temp = 'uinterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            u_pert = squeeze(data);

            %%calculate pressure
            p = repmat(pp00(z0+1:z0+nz_sub),nx_sub,1) + pp_pert;

            %%calculate potential temperature
            th = repmat(th00(z0+1:z0+nz_sub),nx_sub,1) + th_pert;

            %%calculate water vapor mixing ratio
            qv = repmat(qv00(z0+1:z0+nz_sub),nx_sub,1) + qv_pert;    %[kg/kg]

            %%calculate T, Tv, rho, rh
            [T rho rh the] = thermo(p,th,qv);

            %%calculate saturation deficit = Lv*(q*-q)/T
            Lv = 2.5e6;    %[J/kg]
            qvs=qv./rh; %[kg/kg]
            satdef_all=Lv*(qvs-qv)./T;  %satdef everywhere
            satdef = sum(p.*satdef_all,2)./sum(p,2); %pressure-weighted column satdef

            %%some other stuff
            data = satdef;
            var = 'satdef';
            v_def = 'Saturation deficit';
            v_units = 'J kg^{-1} K^{-1}';  

        else
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
            xres = dx;

        end
        
        %%STORE THE WIND SPEED AT THE DESIRED RADIUS
        if(r_ts==1)
            i_max = find(data==max(max(data)));  %index of max wind speed
            i_max = i_max(1);
            if(strmatch(v_usr,'MAX'))
                if(~isempty(i_max))
                    r_usr(ii)=xres*i_max-.5*xres;
                else
                    r_usr(ii)=NaN;
                end
            else
                if(strmatch(var,'vinterp'))
                    temp1=find(data<=v_usr);
                    temp2=find(temp1>i_max,1);
                    i_vusr_out=temp1(temp2);
                    clear temp1 temp2
                    if(~isempty(i_vusr_out))
                        v_out = data(i_vusr_out);   %v at first point where v<=v_usr
                        v_in = data(i_vusr_out-1);  %v at last point where v>=v_usr
                        r_usr(ii) = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                    else
                        v_out = NaN;
                        v_in = NaN;
                        r_usr(ii) = NaN;
                    end
                elseif(strmatch(var,'satdef'))
                    i_min = find(data==min(data(1:25)));  %index of min satdef
                    temp1=find(data>v_usr); %find pts with satdef > v_usr
                    i_vusr_out=temp1(find(temp1>i_min,1));    %first of these outside of eyewall satdef min
                    clear temp1
                    if(~isempty(i_vusr_out))
                        v_out = data(i_vusr_out);   %v at first point where v<=v_usr
                        v_in = data(i_vusr_out-1);  %v at last point where v>=v_usr
                        r_usr(ii) = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                    else
                        v_out = NaN;
                        v_in = NaN;
                        r_usr(ii) = NaN;
                    end
                end
            end
        end
        
        %%CALCULATE TIME-AVERAGED RADIAL PROFILE
        v_def_qv = v_def;
        v_units_qv = v_units;
        rprof_tmean=rprof_tmean(1:nx_sub);
        rprof_tmean=rprof_tmean+data/(i_tf-i_t0+1);
                

        
    end

    
    if(strcmp(var,'vinterp'))
        %%INTERPOLATE TO FIND TRUE Vmax AND rmax (important for low-res runs)
        xvals_all = xres*(.5:1:.5+length(data));
        [Vmaxs(rr) rmaxs(rr)] = VmaxRmax(xvals_all,rprof_tmean);

        %%use pure dataset's Vmax and Rmax
    %    xvals = xmin_sub:dx:xmax_sub;
    %    i_max = find(rprof_tmean==max(rprof_tmean));  %index of max wind speed
    %    Vmaxs(rr) = rprof_tmean(i_max);
    %    rmaxs(rr) = xvals(i_max);

    end
    
    %% PLOT
    %%PLOT THINGS
    %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
    xvals = xmin_sub:dx:xmax_sub;
    yvals = ymin_sub:dy:ymax_sub;
    zvals = zmin_sub:dz:zmax_sub;

    %%BASIC STATISTICS
    rprof_min = min(min(rprof_tmean));
    if(strcmp(var,'vinterp'))
        rprof_max = Vmaxs(rr);
    else
        rprof_max = max(max(rprof_tmean));
    end
    
    %1-DIM RADIAL PLOT
    if(plot_normalized==1)
        figure(2)
        plot(xvals/rmaxs(rr),rprof_tmean'/Vmaxs(rr),pl_clrs{rr},'LineWidth',2)
    %    axis([xmin_sub xmax_sub 0 Vmaxs(rr)])
        %axis([0/rmaxs(rr) rmax_plot/rmaxs(rr) 0 1])
        axis([0/rmaxs(rr) 20 0 1])
        hold on

        %%TIME SERIES OF r_usr
        if(r_ts==1)
            t_day_min = min(t_day);
            t_day_max = max(t_day_max,max(t_day));

            figure(1)
            plot(t_day',smooth(r_usr/rmaxs(rr),numpt_sm)',pl_clrs{rr},'LineWidth',2)
            %axis([t_day_min t_day_max rmin_plot/rmaxs(rr) rmax_plot/rmaxs(rr)])
            axis([t_day_min t_day_max rmin_plot/rmaxs(rr) 20])
            hold on
            smooth_str=sprintf('%i-pt smooth',numpt_sm);
            r_usr_means(rr) = mean(smooth(r_usr/rmaxs(rr),numpt_sm));
        end
        
    else        
        figure(2)
        plot(xvals,rprof_tmean',pl_clrs{rr},'LineWidth',4)
%        axis([xmin_sub xmax_sub 0 max(Vmaxs(rr))])
%        axis([xmin_sub xmax_sub vmin_plot vmax_plot])
        axis([0 rmax_plot vmin_plot vmax_plot])
        hold on
        
        %%TIME SERIES OF r_usr
        if(r_ts==1)
            t_day_min = min(t_day);
            t_day_max = max(t_day_max,max(t_day));

            figure(1)
            plot(t_day',smooth(r_usr,numpt_sm)',pl_clrs{rr},'LineWidth',2)
            axis([t_day_min t_day_max rmin_plot rmax_plot])
            hold on
            smooth_str=sprintf('%i-pt smooth',numpt_sm);
            r_usr_means(rr) = mean(smooth(r_usr,numpt_sm));
        end

    end
            
end
figure(2)
set(gca,'fontweight','bold','fontsize',24)
if(plot_normalized==1)
    input_xlabel=sprintf('r / rmax');
    xlabel(input_xlabel);
    input_ylabel=sprintf('V / Vmax',v_def);
    ylabel(input_ylabel);
    input_title1=sprintf('Normalized azimuthal velocity');
    input_title2=sprintf('z=%5.2f %s',zmin_sub,zunits);
    title({input_title1,input_title2})
else
    input_xlabel=sprintf('Radius [km]',xunits);
    xlabel(input_xlabel);
    input_ylabel=sprintf('Azimuthal velocity [m/s]');
    ylabel(input_ylabel);
    input_title1=sprintf('%s [%s]',v_def,v_units);
    if(strcmp(var,'vinterp'))
        input_title2=sprintf('z=%5.2f %s; from %5.3f to %5.3f',zmin_sub,zunits,rprof_min,max(Vmaxs));
    else
        input_title2=sprintf('z=%5.2f %s; from %5.3f to %5.3f',zmin_sub,zunits,rprof_min,rprof_max);
    end
    title('Modeled radial profiles of azimuthal wind speed ')
end
%input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
legend('Largest storm','Smallest storm')

input_title=strrep(sprintf('time-mean days %i-%i',t0,tf),'_','\_');
%annotation('textbox','String',input_title,'Position',[0 .5 .2 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')

set(findall(gca, 'Type','text'), 'FontSize', 24,'fontweight','bold')

if(r_ts==1)

    figure(1)
    set(gca,'fontweight','bold','fontsize',11)
    if(plot_normalized==1)
        input_xlabel=sprintf('time [day]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('r%5.1f / rmax',v_usr);
        ylabel(input_ylabel);
        input_title1=sprintf('Time series: Normalized radius of%5.1f m/s wind (%s)',v_usr,smooth_str);
        input_title2=sprintf('z=%5.2f %s',zmin_sub,zunits);
        title({input_title1,input_title2})
        
    else
        input_xlabel=sprintf('time [day]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('radius of%5.1f m/s wind [km]',v_usr);
        ylabel(input_ylabel);
        input_title1=sprintf('Time series: Radius of%5.1f m/s wind [km] (%s)',v_usr,smooth_str);
        input_title2=sprintf('z=%5.2f %s; from %5.1f to %5.1f',zmin_sub,zunits,min(r_usr),max(r_usr));
        title({input_title1,input_title2})
    end
    
    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
    legend(input_legend,'Location','SouthEast')

end

%%output some useful data to screen
if(strmatch(var,'vinterp'))
    Vmaxs
    rmaxs
end
v_usr
r_usr_means

%percent_ch_from_1st = (r_usr_means-r_usr_means(1))/r_usr_means(1)*100

    
%}