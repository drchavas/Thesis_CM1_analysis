%nc_plot.m

%Adapted to axisym/3d: 25 Jan 2011, Dan Chavas

%Purpose: to extract subsets of data from a netcdf file and plot them

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

%close all
clear
clc
figure(1)
clf(1)

%% USER INPUT %%%%%%%%%%%%%%%%%%
run_type=1; %1=axisym; 3=3d; 10/30=initial condition test (subdir 'code_test')

subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
%subdir_pre='TRANSFER/';
subdir_pre='STATSTEST/';
%subdir_pre='RCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive


subdir_axisym = 'DAN'; %name of sub-directory with nc files
t_file_axisym=327;  %file timestep
subdir_3d = 'RCE_nx48_SST300.00K_Tthresh200K_usfc3';  %name of sub-directory with nc files
t_file_3d=101;  %file timestep

plot_type = 1;  %0=no plotd
                %1=single plot of [var] in domain defined below
                %4=4-panel:
                    %axisym: vinterp, w, qv, uinterp
                    %3d: vinterp, mid-level w, qv, XZ-xsec qv
    fcor = .00005;  %for angular momentum plot
    
%IF plot_type=1: input variable and domain
    fig_hold=1; %1=do not clear old figure
    pl_color='r--';
    var = 'vinterp';  %variable of interest (DRC vars: 'angmom')
    usr_maxmin=0;   %0=use data max/min for plot; 1=usr input
        vmin_plot = -40;
        vmax_plot = 40;
    
 if(run_type==1 || run_type==10)   %axisym
    x0=0;   %first x grid point [0,end]
    xf=1000;   %first y grid point [0,end]
    y0=0;   %first y grid point [0,end]
    yf=0;   %last y grid point [0,end]
    z0=1;  %first z grid point [0,end]
    zf=1;  %last z grid point [0,end]
 elseif(run_type==3 || run_type==30)   %3d
    x0=-1000;   %first x grid point [0,end]
    xf=1000;   %first y grid point [0,end]
    y0=0;   %first y grid point [0,end]
    yf=0;   %last y grid point [0,end]
    z0=0;  %first z grid point [0,end]
    zf=0;  %last z grid point [0,end]
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%CLEAR OLD FIGURE?
if(fig_hold==1)
    hold on
else
    clf     %clear figures without closing them
end

quit_prog=0;
%% OPTIONS FOR EITHER AXISYM OR 3d RUNS
if(run_type==1)
    run_type_str='axisym';
    if(ext_hd==0)
        dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/axisym/CM1_output/%s',subdir_pre);
    else    %external harddrive
        dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/axisym/%s',subdir_pre);
    end
    subdir=subdir_axisym;
    t_file=t_file_axisym;
    copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract       
elseif(run_type==3)
    run_type_str='3d';
    if(ext_hd==0)
        dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3d/CM1_output/%s',subdir_pre);
    else    %external harddrive
        dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/3d/%s',subdir_pre);
    end
    subdir=subdir_3d;
    t_file=t_file_3d;
    copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
elseif(run_type==10)
    run_type_str='axisym';
    dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/code_test/CM1_output/%s',subdir_pre);
    subdir=subdir_axisym;
    t_file=t_file_axisym;
    copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
elseif(run_type==30)
    run_type_str='3d';
    dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/code_test/CM1_output/%s',subdir_pre);
    subdir=subdir_3d;
    t_file=t_file_3d;
    copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract    
else
    sprintf('ERROR: Input value for run_type must be 1 (axisym) or 3 (3d).')
    quit_prog=1;
end

if(quit_prog==0)

%%NC FILENAME (default)
if(t_file<10)
    nc_file = sprintf('cm1out_000%i.nc',t_file);
elseif(t_file<100)
    nc_file = sprintf('cm1out_00%i.nc',t_file);
else
    nc_file = sprintf('cm1out_0%i.nc',t_file);
end

%% OPEN FILE AND EXTRACT nx, ny, nz
dir_start=pwd;
dir_full = strcat(dir_in,subdir);
cd(dir_full)
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

%DEFINE y-midpoint FOR XZ cross-sections 
y_mid=round(ny/2);


%% DEFINE BOUNDARIES FOR plot_type=4 FOR EITHER AXISYM OR 3d RUNS
if(run_type==1 || run_type==10)
    %for plot=4: domain subset limits
    x01=0;   %first x grid point [0,end]
    xf1=400;   %first y grid point [0,end]
    y01=0;   %first y grid point [0,end]
    yf1=0;   %last y grid point [0,end]
    z01=0;  %first z grid point [0,end]
    zf1=35;  %last z grid point [0,end]

    x02=0;   %first x grid point [0,end]
    xf2=300;   %first y grid point [0,end]
    y02=0;   %first y grid point [0,end]
    yf2=0;   %last y grid point [0,end]
    z02=0;  %first z grid point [0,end]
    zf2=25;  %last z grid point [0,end]

    x03=0;   %first x grid point [0,end]
    xf3=400;   %first y grid point [0,end]
    y03=0;   %first y grid point [0,end]
    yf3=0;   %last y grid point [0,end]
    z03=0;  %first z grid point [0,end]
    zf3=35;  %last z grid point [0,end]

    x04=0;   %first x grid point [0,end]
    xf4=400;   %first y grid point [0,end]
    y04=0;   %first y grid point [0,end]
    yf4=0;   %last y grid point [0,end]
    z04=0;  %first z grid point [0,end]
    zf4=35;  %last z grid point [0,end]
    
elseif(run_type==3)
    %for plot=4: domain subset limits
    x01=0;   %first x grid point [0,end]
    xf1=5000;   %first y grid point [0,end]
    y01=y_mid;   %first y grid point [0,end]
    yf1=y_mid;   %last y grid point [0,end]
    z01=0;  %first z grid point [0,end]
    zf1=15;  %last z grid point [0,end]

    x02=0;   %first x grid point [0,end]
    xf2=5000;   %first y grid point [0,end]
    y02=y_mid;   %first y grid point [0,end]
    yf2=y_mid;   %last y grid point [0,end]
    z02=0;  %first z grid point [0,end]
    zf2=15;  %last z grid point [0,end]

    x03=0;   %first x grid point [0,end]
    xf3=5000;   %first y grid point [0,end]
    y03=y_mid;   %first y grid point [0,end]
    yf3=y_mid;   %last y grid point [0,end]
    z03=0;  %first z grid point [0,end]
    zf3=15;  %last z grid point [0,end]

    x04=0;   %first x grid point [0,end]
    xf4=5000;   %first y grid point [0,end]
    y04=y_mid;   %first y grid point [0,end]
    yf4=y_mid;   %last y grid point [0,end]
    z04=0;  %first z grid point [0,end]
    zf4=8;  %last z grid point [0,end]
    
else
    sprintf('ERROR: Input value for run_type must be 1 (axisym) or 3 (3d).')
    quit_prog=1;
end

%% PLOT
switch plot_type
    case 0  %no plot

    case 1  %single plot of [var] over specified domain
        %%EXTRACT DATA
        
        if(strcmp(var,'angmom'))
            %%INPUT INFO
            var = 'vinterp';  %variable of interest

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy nx_sub ny_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
            v_def = 'absolute angular momentum';
            v_units = 'm^2/s';
            
            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA: absolute angular momentum
%            pcolor(xmat,zmat,double(data.*xmat+.5*fcor*xmat.^2))
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
%            shading flat
%            shading interp
            contour(xmat,zmat,double(data.*xmat+.5*fcor*xmat.^2),40)
            colorbar            
            set(gca,'fontweight','bold','fontsize',11)
            
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})
            
        elseif(strcmp(var,'strmfnc'))

            snd_file = 'input_sounding';

            %%load initial sounding variables
            [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

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

            %%load qcpert
            var_temp = 'qc';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qc_pert = squeeze(data);

            %%load qrpert
            var_temp = 'qr';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qr_pert = squeeze(data);

            %%load qipert
            var_temp = 'qi';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qi_pert = squeeze(data);

            %%load qspert
            var_temp = 'qs';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qs_pert = squeeze(data);

            %%load qgpert
            var_temp = 'qg';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qg_pert = squeeze(data);

            %%load upert
            var_temp = 'uinterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            u_pert = squeeze(data);

            %%load wpert
            var_temp = 'winterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            w_pert = squeeze(data);

            
            %%calculate pressure
            p = repmat(pp00(z0+1:zf+1)',nx_sub,1) + pp_pert;

            %%calculate potential temperature
            th = repmat(th00(z0+1:zf+1)',nx_sub,1) + th_pert;

            %%calculate water vapor mixing ratio
            qv = repmat(qv00(z0+1:zf+1)',nx_sub,1) + qv_pert;    %[kg/kg]
            
            %%calculate liquid water mixing ratio
            ql = ql_pert;    %[kg/kg]

            %%calculate T, Tv, rho, rh
            [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);
            
            %%calculate radial velocity
            u = repmat(u00(z0+1:zf+1)',nx_sub,1) + u_pert;

            %%calculate vertical velocity
            w = w_pert;

            %%calculate streamfunction
            r = xmin_sub:dx:xmax_sub;
            r=r';
            r=repmat(r,1,nz_sub);
            psi = flowfun(r.*rho.*u,r.*rho.*w,'-');
            figure(4)
            contour(psi)
            colorbar
            
            
            
            figure(2)
            data = rh;

            
            
            
        else
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
        end
        
        %%TIME IN DAYS
        t_day = time / 86400;
        
        %%BASIC STATISTICS
        data_min = min(min(data));
        data_max = max(max(data));

        %%PLOT THINGS
        %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
        xvals = xmin_sub:dx:xmax_sub;
        yvals = ymin_sub:dy:ymax_sub;
        zvals = zmin_sub:dz:zmax_sub;

        %ORIENT DATA CORRECTLY FOR PLOTTING
        data=squeeze(data)'; %now x-data is along row and y-data is along column
        data=double(data);  %pcolor requires double format
        
        %%PLOT DATA
        %1-DIM PLOTS
        if(usr_maxmin==0)
            vmin_plot = min(min(data));
            vmax_plot = max(max(data));
        end
        
        if(nx_sub>1 && ny_sub==1 && nz_sub==1)  %X-only plot
            %%DRC 01/31/12
            %data = (data(3:end)-data(1:end-2))./(2*1000*dx);
            %data = ((1000*xvals(3:end).*data(3:end)-1000*xvals(1:end-2).*data(1:end-2))./(2*1000*dx))./(1000*xvals(2:end-1)); %for uinterp
            %plot(xvals(2:end-1),data',pl_color,'LineWidth',2)
            %end DRC
            
            plot(xvals,data',pl_color,'LineWidth',1)
            axis([xmin_sub xmax_sub vmin_plot vmax_plot])
            set(gca,'fontweight','bold','fontsize',11)
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('%s [%s]',v_def,v_units);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s, z=%5.2f %s; from %5.3f to %5.3f',t_day,xmin_sub,xunits,zmin_sub,zunits,data_min,data_max);
            title({input_title1,input_title2})
            
        elseif(ny_sub>1 && nx_sub==1 && nz_sub==1)  %Y-only plot
            plot(yvals,data',pl_color,'LineWidth',2)
            axis([ymin_sub ymax_sub vmin_plot vmax_plot])
            set(gca,'fontweight','bold','fontsize',11)
            input_xlabel=sprintf('y-distance [%s]',yunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('%s [%s]',v_def,v_units);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, x=%5.2f %s, z=%5.2f %s; from %5.3f to %5.3f',t_day,xmin_sub,xunits,zmin_sub,zunits,data_min,data_max);
            title({input_title1,input_title2})
            
        elseif(nz_sub>1 && nx_sub==1 && ny_sub==1)  %Z-only plot
            plot(data',zvals,pl_color,'LineWidth',2)
            axis([vmin_plot vmax_plot zmin_sub zmax_sub])
            set(gca,'fontweight','bold','fontsize',11)
            input_xlabel=sprintf('%s [%s]',v_def,v_units);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, x=%5.2f %s, y=%5.2f %s; from %5.3f to %5.3f',t_day,xmin_sub,xunits,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})
        
        %2-DIM PLOTS
        elseif(ny_sub==1 && nx_sub>1 && nz_sub>1)   %XZ cross-section
            %MAKE X AND Z MATRICES
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            %pcolor has x and y switched vs. geographic coordinates 
            data=squeeze(data); %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            pcolor(xmat,zmat,data)
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            caxis([vmin_plot vmax_plot]); 
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})

        elseif(nx_sub==1 && ny_sub>1 && nz_sub>1)   %YZ cross-section
            %MAKE Y AND Z MATRICES
            ymat = repmat(yvals,length(zvals),1);
            zmat = repmat(zvals,length(yvals),1)';

            %PLOT DATA
            figure(1)
            pcolor(ymat,zmat,data)
            axis([ymin_sub ymax_sub zmin_sub zmax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            caxis([vmin_plot, vmax_plot]); 
            input_xlabel=sprintf('y-distance [%s]',yunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, x=%5.2f %s; from %5.3f to %5.3f',t_day,xmin_sub,xunits,data_min,data_max);
            title({input_title1 input_title2})
            
        elseif(nz_sub==1 && nx_sub>1 && ny_sub>1)   %XY plane
            %MAKE X AND Y MATRICES
            xmat = repmat(xvals,length(yvals),1);
            ymat = repmat(yvals,length(xvals),1)';

            %PLOT DATA
            figure(1)
            pcolor(xmat,ymat,data)
            axis([xmin_sub xmax_sub ymin_sub ymax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            caxis([vmin_plot, vmax_plot]); 
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('y-distance [%s]',yunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, z=%5.2f %s; from %5.3f to %5.3f',t_day,zmin_sub,zunits,data_min,data_max);
            title({input_title1 input_title2})

        else    
            sprintf('Cant make plot, check subset dimensions.')
        end
      
        %% FIGURE TITLE
        fig_title = sprintf('%s, %s',subdir,run_type_str);
        annotation('textbox','String',fig_title,'Position',[0 .5 .4 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')
            
    case 4  %4-panel plot: rain, mid-level w, qv, XZ-xsec qv
%        if(run_type==1) %axisym: vinterp, w, qv, uinterp

            %% Calculate various thermodynamic quantities
            snd_file = 'input_sounding';

            %%load initial sounding variables
            [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

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

            %%load qcpert
            var_temp = 'qc';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qc_pert = squeeze(data);

            %%load qrpert
            var_temp = 'qr';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qr_pert = squeeze(data);

            %%load qipert
            var_temp = 'qi';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qi_pert = squeeze(data);

            %%load qspert
            var_temp = 'qs';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qs_pert = squeeze(data);

            %%load qgpert
            var_temp = 'qg';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qg_pert = squeeze(data);
            
            %%load upert
            var_temp = 'uinterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            u_pert = squeeze(data);

            %%calculate pressure
            p = repmat(pp00(z0+1:z0+nz_sub-1),nx_sub,1) + pp_pert(:,1:end-1);

            %%calculate potential temperature
            th = repmat(th00(z0+1:z0+nz_sub-1),nx_sub,1) + th_pert(:,1:end-1);

            %%calculate water vapor mixing ratio
            qv = repmat(qv00(z0+1:z0+nz_sub-1),nx_sub,1) + qv_pert(:,1:end-1);    %[kg/kg]

            %%calculate liquid water mixing ratio
            ql = ql_pert;    %[kg/kg]

            %%calculate T, Tv, rho, rh
            [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);
           

            %% PANEL 1: vinterp XZ x-sec
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'vinterp';  %variable of interest
            x0=x01;   %first x grid point [0,end]
            xf=xf1;   %first y grid point [0,end]
            y0=y01;   %first y grid point [0,end]
            yf=yf1;   %last y grid point [0,end]
            z0=z01;  %first z grid point [0,end]
            zf=zf1;  %last z grid point [0,end]

            %%EXTRACT DATA
            pwd
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
                        
            %[xx yy]=find(data==max(max(data)))
            %%TIME IN DAYS
            t_day = time / 86400;
            
            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            hold off
            subplot(2,2,1)
            pcolor(xmat,zmat,data)
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
            %caxis([0 70])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})


            %% PANEL 2: the XZ x-sec
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
%            var = 'thpert';  %variable of interest
            x0=x02;   %first x grid point [0,end]
            xf=xf2;   %first y grid point [0,end]
            y0=y02;   %first y grid point [0,end]
            yf=yf2;   %last y grid point [0,end]
            z0=z02;  %first z grid point [0,end]
            zf=zf2;  %last z grid point [0,end]

%% Calculate various thermodynamic quantities
            snd_file = 'input_sounding';

            %%load initial sounding variables
            [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

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

            %%load qcpert
            var_temp = 'qc';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qc_pert = squeeze(data);

            %%load qrpert
            var_temp = 'qr';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qr_pert = squeeze(data);

            %%load qipert
            var_temp = 'qi';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qi_pert = squeeze(data);

            %%load qspert
            var_temp = 'qs';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qs_pert = squeeze(data);

            %%load qgpert
            var_temp = 'qg';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qg_pert = squeeze(data);

            %%ql_pert
            ql_pert = qc_pert + qr_pert + qi_pert + qs_pert + qg_pert;
            
            %%load upert
            var_temp = 'uinterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            u_pert = squeeze(data);

            %%calculate pressure
            p = repmat(pp00(z0+1:z0+nz_sub-1),nx_sub,1) + pp_pert(:,1:end-1);

            %%calculate potential temperature
            th = repmat(th00(z0+1:z0+nz_sub-1),nx_sub,1) + th_pert(:,1:end-1);

            %%calculate water vapor mixing ratio
            qv = repmat(qv00(z0+1:z0+nz_sub-1),nx_sub,1) + qv_pert(:,1:end-1);    %[kg/kg]

            %%calculate liquid water mixing ratio
            ql = ql_pert;    %[kg/kg]
            
            %%calculate T, Tv, rho, rh
            [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);
            
            %%TIME IN DAYS
            t_day = time / 86400;
            
            for vv = 1:1
                if(vv==1)
                    data = the;
                    var = 'the';
                    v_def = 'Equivalent potential temperature';
                    v_units = 'K';
                elseif(vv==2)
                    data = rh;
                    var = 'rh';
                    v_def = 'Relative humidity';
                    v_units = '';  
                end
            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals(1:end-1)),1);
            zmat = repmat(zvals(1:end-1),length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            hold on
            subplot(2,2,2)
            if(vv==1)
%            pcolor(xmat,zmat,data)
%            caxis([305 360])
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
%            shading flat
%            shading interp
            contour(xmat,zmat,data,20)
            colorbar            
            set(gca,'fontweight','bold','fontsize',11)
            elseif(vv==2)
            contour(xmat,zmat,data)
            end
            
            end
            
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})

            %% PANEL 3: winterp XZ x-sec
%{
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'winterp';  %variable of interest
            x0=x03;   %first x grid point [0,end]
            xf=xf3;   %first y grid point [0,end]
            y0=y03;   %first y grid point [0,end]
            yf=yf3;   %last y grid point [0,end]
            z0=z03;  %first z grid point [0,end]
            zf=zf3;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            hold on
            subplot(2,2,3)
            pcolor(xmat,zmat,data)
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})
%}
            
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'vinterp';  %variable of interest
            x0=x03;   %first x grid point [0,end]
            xf=xf3;   %first y grid point [0,end]
            y0=y03;   %first y grid point [0,end]
            yf=yf3;   %last y grid point [0,end]
            z0=z03;  %first z grid point [0,end]
            zf=zf3;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
            v_def = 'absolute angular momentum';
            v_units = 'm^2/s';
            
            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA: absolute angular momentum
            hold on
            subplot(2,2,3)
%            pcolor(xmat,zmat,double(data.*xmat+.5*fcor*xmat.^2))
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
%            shading flat
%            shading interp
            contour(xmat,zmat,double(data.*xmat+.5*fcor*xmat.^2),40)
            colorbar            
            set(gca,'fontweight','bold','fontsize',11)
            
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})

            %% PANEL 4: rh XZ x-sec
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            x0=x04;   %first x grid point [0,end]
            xf=xf4;   %first y grid point [0,end]
            y0=y04;   %first y grid point [0,end]
            yf=yf4;   %last y grid point [0,end]
            z0=z04;  %first z grid point [0,end]
            zf=zf4;  %last z grid point [0,end]
            
            %% Calculate various thermodynamic quantities
            snd_file = 'input_sounding';

            %%load initial sounding variables
            [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

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

            %%load qcpert
            var_temp = 'qc';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qc_pert = squeeze(data);

            %%load qrpert
            var_temp = 'qr';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qr_pert = squeeze(data);

            %%load qipert
            var_temp = 'qi';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qi_pert = squeeze(data);

            %%load qspert
            var_temp = 'qs';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qs_pert = squeeze(data);

            %%load qgpert
            var_temp = 'qg';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            qg_pert = squeeze(data);
            
            %%load upert
            var_temp = 'uinterp';
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
            u_pert = squeeze(data);

            %%calculate pressure
            p = repmat(pp00(z0+1:z0+nz_sub-1),nx_sub,1) + pp_pert(:,1:end-1);

            %%calculate potential temperature
            th = repmat(th00(z0+1:z0+nz_sub-1),nx_sub,1) + th_pert(:,1:end-1);

            %%calculate water vapor mixing ratio
            qv = repmat(qv00(z0+1:z0+nz_sub-1),nx_sub,1) + qv_pert(:,1:end-1);    %[kg/kg]

            %%calculate liquid water mixing ratio
            ql = ql_pert;    %[kg/kg]
            
            %%calculate T, Tv, rho, rh
            [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);
            
            %{
            %%calculate radial velocity
            u = repmat(u00(z0+1:zf+1)',nx_sub,1) + u_pert;

            %%calculate vertical velocity
            w = w_pert;

            %%calculate streamfunction
            r = xmin_sub:dx:xmax_sub;
            r=r';
            r=repmat(r,1,nz_sub);
            psi = flowfun(r.*rho.*u,r.*rho.*w,'-');
            contour(psi)
            colorbar
            %}
            
            %%TIME IN DAYS
            t_day = time / 86400;

            for vv = 2:2
                if(vv==1)
                    data = the;
                    var = 'the';
                    v_def = 'Equivalent potential temperature';
                    v_units = 'K';
                elseif(vv==2)
                    data = rh;
                    var = 'rh';
                    v_def = 'Relative humidity';
                    v_units = '';  
                end
            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals(1:end-1)),1);
            zmat = repmat(zvals(1:end-1),length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            hold on
            subplot(2,2,4)
%            if(vv==1)
            pcolor(xmat,zmat,data)
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
%            elseif(vv==2)
%            contour(xmat,zmat,data)
%            end
            
            end
            
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s]',v_def,v_units);
            input_title2=sprintf('t=%5.3f days, y=%5.2f %s; from %5.3f to %5.3f',t_day,ymin_sub,yunits,data_min,data_max);
            title({input_title1,input_title2})
            
            


            %% FIGURE TITLE
            fig_title = sprintf('%s, %s',subdir,run_type_str);
            annotation('textbox','String',fig_title,'Position',[.1 .5 .8 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')

    case 41 %3d ONLY: VARIOUS PLOTS
            
        if(run_type==1 || run_type==10)
            sprintf('This plot_type works for 3d data only.')
        elseif(run_type==3)
        %{
            %% PANEL 1: rain
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'rain';  %variable of interest
            x0=0;   %first x grid point [0,end]
            xf=10000;   %first y grid point [0,end]
            y0=0;   %first y grid point [0,end]
            yf=10000;   %last y grid point [0,end]
            z0=0;  %first z grid point [0,end]
            zf=0;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals yvals xmat ymat
            xvals = xmin_sub:dx:xmax_sub;
            yvals = ymin_sub:dy:ymax_sub;
            xmat = repmat(xvals,length(yvals),1);
            ymat = repmat(yvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            title('test')
            hold off
            subplot(2,2,1)
            pcolor(xmat,ymat,data)
            axis([xmin_sub xmax_sub ymin_sub ymax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('y-distance [%s]',yunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s] @ z = %5.2f %s',var,v_units,zmin_sub,zunits);
            input_title2=sprintf('t=%5.3f days, z=%5.2f %s; from %5.3f to %5.3f',t_day,zmin_sub,zunits,data_min,data_max);
            title({input_title1 input_title2})
            %}
            
            %% PANEL 1: vinterp
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'vinterp';  %variable of interest
            x0=0;   %first x grid point [0,end]
            xf=200;   %first y grid point [0,end]
            y0=0;   %first y grid point [0,end]
            yf=200;   %last y grid point [0,end]
            z0=0;  %first z grid point [0,end]
            zf=0;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals yvals xmat ymat
            xvals = xmin_sub:dx:xmax_sub;
            yvals = ymin_sub:dy:ymax_sub;
            xmat = repmat(xvals,length(yvals),1);
            ymat = repmat(yvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            title('test')
            hold off
            subplot(2,2,1)
            pcolor(xmat,ymat,data)
            axis([xmin_sub xmax_sub ymin_sub ymax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('y-distance [%s]',yunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s] @ z = %5.2f %s',var,v_units,zmin_sub,zunits);
            input_title2=sprintf('t=%5.3f days, z=%5.2f %s; from %5.3f to %5.3f',t_day,zmin_sub,zunits,data_min,data_max);
            title({input_title1 input_title2})

            %% PANEL 2: winterp
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'winterp';  %variable of interest
            x0=70;   %first x grid point [0,end]
            xf=130;   %first y grid point [0,end]
            y0=70;   %first y grid point [0,end]
            yf=130;   %last y grid point [0,end]
            z0=4;  %first z grid point [0,end]
            zf=4;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals yvals xmat ymat
            xvals = xmin_sub:dx:xmax_sub;
            yvals = ymin_sub:dy:ymax_sub;
            xmat = repmat(xvals,length(yvals),1);
            ymat = repmat(yvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            hold on
            subplot(2,2,2)
            pcolor(xmat,ymat,data)
            axis([xmin_sub xmax_sub ymin_sub ymax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('y-distance [%s]',yunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s] @ z = %5.2f %s',var,v_units,zmin_sub,zunits);
            input_title2=sprintf('t=%5.3f days, z=%5.2f %s; from %5.3f to %5.3f',t_day,zmin_sub,zunits,data_min,data_max);
            title({input_title1 input_title2})

            %% PANEL 3: qv
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'qv';  %variable of interest
            x0=70;   %first x grid point [0,end]
            xf=130;   %first y grid point [0,end]
            y0=70;   %first y grid point [0,end]
            yf=130;   %last y grid point [0,end]
            z0=4;  %first z grid point [0,end]
            zf=4;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals yvals xmat ymat
            xvals = xmin_sub:dx:xmax_sub;
            yvals = ymin_sub:dy:ymax_sub;
            xmat = repmat(xvals,length(yvals),1);
            ymat = repmat(yvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            subplot(2,2,3)
            pcolor(xmat,ymat,data)
            axis([xmin_sub xmax_sub ymin_sub ymax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('y-distance [%s]',yunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s] @ z = %5.2f %s',var,v_units,zmin_sub,zunits);
            input_title2=sprintf('t=%5.3f days, z=%5.2f %s; from %5.3f to %5.3f',t_day,zmin_sub,zunits,data_min,data_max);
            title({input_title1 input_title2})

            %% PANEL 4: qv XZ cross-section
            %%INPUT INFO
            clear var x0 xf y0 yf z0 zf
            var = 'qv';  %variable of interest
            x0=70;   %first x grid point [0,end]
            xf=130;   %first y grid point [0,end]
            y0=100;   %first y grid point [0,end]
            yf=100;   %last y grid point [0,end]
            z0=0;  %first z grid point [0,end]
            zf=5;  %last z grid point [0,end]

            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            %%TIME IN DAYS
            t_day = time / 86400;

            %%BASIC STATISTICS
            clear data_min data_max
            data_min = min(min(data));
            data_max = max(max(data));

            %%PLOT THINGS
            %MAKE MATRIX OF x AND y VALUES FOR PCOLOR
            clear xvals zvals xmat zmat
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1);
            zmat = repmat(zvals,length(xvals),1)';

            %ORIENT DATA CORRECTLY FOR PLOTTING
            data=squeeze(data)'; %now x-data is along row and y-data is along column
            data=double(data);  %pcolor requires double format

            %PLOT DATA
            figure(1)
            subplot(2,2,4)
            pcolor(xmat,zmat,data)
            axis([xmin_sub xmax_sub zmin_sub zmax_sub])
            shading flat
            shading interp
            set(gca,'fontweight','bold','fontsize',11)
            colorbar
            input_xlabel=sprintf('x-distance [%s]',xunits);
            xlabel(input_xlabel);
            input_ylabel=sprintf('height [%s]',zunits);
            ylabel(input_ylabel);
            input_title1=sprintf('%s [%s] @ y = %5.2f %s',var,v_units,ymin_sub,yunits);
            input_title2=sprintf('from %5.3f to %5.3f',data_min,data_max);
            title({input_title1 input_title2})

            %% FIGURE TITLE
            fig_title = sprintf('%s, %s',subdir,run_type_str);
            annotation('textbox','String',fig_title,'Position',[.1 .5 .8 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')
        end
    otherwise
        sprintf('Invalid plot number. Fail.')
end

end

%}