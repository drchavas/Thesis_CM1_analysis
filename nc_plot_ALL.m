%nc_plot_steady_xz.m

%Created: 10-05-11, Dan Chavas

%Purpose: This file plots any subset of data (X,Y,Z,XY,XZ,YZ) averaged over
%any time interval.  Very nice.

clear
clc
figure(1)
hold on
figure(2)
clf(1)     %clear figures without closing them
clf(2)     %clear figures without closing them


%% USER INPUT %%%%%%%%%%%%%%%%%%
%subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
subdir_pre='RCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_type=3; %1=axisym; 3=3d

t0 = 99;
tf = 100;

subdir = 'RCE_nx48_SST300.00K_Tthresh200K_usfc3_drag'; %name of sub-directory with nc files

plot_type = 1;  %0=no plot
                %1=single plot of time series of [var] for above runs
                %4=4-panel:
                    %axisym: vinterp, w, qv, uinterp
                    %3d: vinterp, mid-level w, qv, XZ-xsec qv
%IF plot_type=1: input variable and domain
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    vars = {'p'} %'vinterp' 'Ri' 'vg' 'T' 'rho' 'rh' 'the' 's' 's_sat' 'thv' 'p' 'th' 'qv' 'ql' 'qc' 'qr' 'qi' 'qs' 'qg'};  %pcolor; contour on top

%%Define subset of points to plot
xmin_plot = -1000;  %-1000 [km]; lowest value plotted
xmax_plot = 1000;    %1000 [km]; highest value plotted
ymin_plot = -1000;  %-1000 [km]; lowest value plotted
ymax_plot = 1000;    %1000 [km]; highest value plotted
zmin_plot = 0;  %0 [km]; lowest value plotted
zmax_plot = 20; %1000 [km]; highest value plotted
datamin_plot = -500000;   %-5000 minimum data value plotted
datamax_plot = 500000;   %500000 minimum data value plotted

%%Spatial averaging?
ave_x = 1;  %average data in x
ave_y = 1;  %average data in y
ave_z = 0;  %average data in z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=100000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=100000;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=1000;  %last z grid point [0,end]


numvars = length(vars);

%%Append ax or 3d to subdir names for plot
if(run_type==1) %ax
    subdir_plot=sprintf('ax%s',subdir)
else    %3d
    subdir_plot=sprintf('3d%s',subdir)
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;

%% OPTIONS FOR EITHER AXISYM OR 3d RUNS
if(run_type==1)
    run_type_str='axisym';
    copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
elseif(run_type==3)
    run_type_str='3d';
    copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
end

switch ext_hd
    case 0,
        dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/%s/CM1_output/%s',run_type_str,subdir_pre);
    case 1,
        dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/%s/%s',run_type_str,subdir_pre);
    case 2,
        dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL_nodrag/CM1_output/%s/%s',run_type_str,subdir_pre);
    otherwise
        assert(2==3,'Invalid number for ext_hd!')
end 

%%DIRECTORY WITH OUTPUT DATA
subdir_full=sprintf('%s%s',dir_in,subdir)

%%EXTRACT TIMESTEP SIZE
var_dt = 'qvpert'; %doesn't matter, just need any variable
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
numfiles=length(dir(sprintf('%s/cm1out_0*.nc',subdir_full)));
[data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);

dt = time; %[s]; time of cm1out_0001.nc is defined as zero
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy xunits yunits zunits v_def v_units time t_units
i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

%%TIME IN DAYS
%    t_min_day = dt*i_t0/86400;
%    t_max_day = dt*i_tf/86400;

%%DIRECTORY WITH OUTPUT DATA
subdir_full=sprintf('%s%s',dir_in,subdir)

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding';

%%load initial sounding variables
dir_full = strcat(dir_in,subdir);
[zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

%%initialize the output vectors
clear v_df_qv v_units_qv data_tmean
for vv=1:numvars
    data_tmean{vv}=zeros(nx_sub,ny_sub,nz_sub);
end
rho_tmean=zeros(nx_sub,ny_sub,nz_sub);


clear t_day
clear r_usr
for ii=1:i_tf-i_t0+1
    fclose('all')
    t_file=i_t0+ii-1;

    t_day(ii) = (t_file-1)*dt/(24*60*60);


    %%NC FILENAME (default)
    if(t_file<10)
        nc_file = sprintf('cm1out_000%i.nc',t_file);
    elseif(t_file<100)
        nc_file = sprintf('cm1out_00%i.nc',t_file);
    elseif(t_file<1000)
        nc_file = sprintf('cm1out_0%i.nc',t_file);
    else
        nc_file = sprintf('cm1out_%i.nc',t_file);
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
 
    %% Extract initial sounding data needed to calculate full (rather than just perturbation) values for analysis
    %%load prspert
    var_temp = 'prspert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    pp_pert = data;

    %%load thpert
    var_temp = 'thpert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    th_pert = data;

    %%load qvpert
    var_temp = 'qvpert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qv_pert = data;

    %%load qcpert
    var_temp = 'qc';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qc_pert = data;

    %%load qrpert
    var_temp = 'qr';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qr_pert = data;

    %%load qipert
    var_temp = 'qi';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qi_pert = data;

    %%load qspert
    var_temp = 'qs';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qs_pert = data;

    %%load qgpert
    var_temp = 'qg';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qg_pert = data;
    
    %%ql_pert
    ql = qc_pert + qr_pert + qi_pert + qs_pert + qg_pert;
    
    %%load upert
    var_temp = 'uinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    u_pert = data;

    %%load vpert
    var_temp = 'vinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    v_pert = data;

    %%calculate pressure
    pp00=pp00(z0+1:z0+nz);
    [prof_3d] = sounding_to_3d(pp00,nx_sub,ny_sub);    %expand vertical profile into 3d matrix of same size as data
    p = prof_3d + pp_pert;
    clear prof_3d

    %%calculate potential temperature
    th00=th00(z0+1:z0+nz);  %keep portion of sounding within domain requested by user
    [prof_3d] = sounding_to_3d(th00,nx_sub,ny_sub);    %expand vertical profile into 3d matrix of same size as data
    th = prof_3d + th_pert;
    clear prof_3d

    %%calculate water vapor mixing ratio
    qv00=qv00(z0+1:z0+nz);  %keep portion of sounding within domain requested by user
    [prof_3d] = sounding_to_3d(qv00,nx_sub,ny_sub);    %expand vertical profile into 3d matrix of same size as data
    qv = prof_3d + qv_pert;
    clear prof_3d

    %%calculate various thermodynamic quantities
    [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);

    
    %%Calculate final variable data
    for vv = 1:numvars
        var = vars{vv};
        if(strcmp(var,'Ri'))    %composite variable: Richardson Number
            %%NOTE: size(s)=[384 40]; col 1 = lowest model level
            ds_dz = NaN*s;   %initialize as NaN
            ds_dz(:,2:end-1) = (s(:,3:end)-s(:,1:end-2))/(2000*dz); %[J/K/m]

            %%calculate wind speeds
            u = repmat(u00(z0+1:z0+nz_sub),nx_sub,1) + u_pert;
            v = repmat(v00(z0+1:z0+nz_sub),nx_sub,1) + v_pert;

            du_dz = (u(:,3:end)-u(:,1:end-2))/(2*1000*dz); %[s-1]
            du_dz(:,2:end+1)=du_dz;  %add extra level at bottom equal to nearest level value
            du_dz(:,end+1)=du_dz(:,end);  %add extra level at top equal to nearest level value
            dv_dz = (v(:,3:end)-v(:,1:end-2))/(2*1000*dz); %[s-1]
            dv_dz(:,2:end+1)=dv_dz;  %add extra level at bottom equal to nearest level value
            dv_dz(:,end+1)=dv_dz(:,end);  %add extra level at top equal to nearest level value

            data = gam_m.*ds_dz./(du_dz.^2 + dv_dz.^2);    %Richardson number, eqn (23) Emanuel11
            data(data>10) = NaN;
            data(data<0) = NaN;
            v_defs{vv} = 'Richardson Number';
            v_unitss{vv} = '-';

        elseif(strcmp(var,'angmom'))    %composite variable: angular momentum
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);
            
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1)';
            data = squeeze(data).*xmat*1000 + .5*fcor*(xmat*1000).^2;   %M=rV + .5fr^2
            %data = s;
            %v_defs{vv} = 'entropy';
            %v_unitss{vv} = '-';
            v_defs{vv} = 'Angular momentum';
            v_unitss{vv} = 'm^2/s^2';

        elseif(strcmp(var,'vg'))    %composite variables: gradient wind
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);
            
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1)';
            dp_dr = (p(3:end,:)-p(1:end-2,:))/(2*1000*dx); %[m/s^2]
            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value
            dp_dr(end+1,:)=dp_dr(end,:);  %add extra level at outer edge equal to nearest level value

            %calculate at mid-points between grid points
%            dp_dr = (p(2:end,:)-p(1:end-1,:))/(1000*dx); %[m/s^2]
%            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value

            
            vinterp = data;
            temp2=(xmat*1000*fcor).^2+4*xmat*1000.*dp_dr;
            data =  .5*(-xmat*1000*fcor+sqrt((xmat*1000*fcor).^2+4*xmat*1000.*dp_dr));   %Vg^2 + r*f*Vg-r*(dp_dr); solve for Vg
            data = real(data);  %remove imaginary numbers, which are always very small
            %data(imag(data)~=0)=NaN;    %ignore imaginary numbers
%            data = data - squeeze(vinterp);
            v_defs{vv} = 'Gradient wind';
            v_unitss{vv} = 'm/s';
            
            
            
        elseif(strcmp(var,'rho'))    %composite variables: full density
            data = rho;
            v_defs{vv} = 'Density';
            v_unitss{vv} = 'kg/m3';

        elseif(strcmp(var,'T'))    %composite variables: Temperature
            data = T;
            v_defs{vv} = 'Temperature';
            v_unitss{vv} = 'K';

        elseif(strcmp(var,'rh'))    %composite variables: relative humidity
            data = rh;
            v_defs{vv} = 'Relative humidity';
            v_unitss{vv} = '%';

        elseif(strcmp(var,'the'))    %composite variables: theta-e
            data = the;
            v_defs{vv} = 'Theta-e';
            v_unitss{vv} = 'K';
            
        elseif(strcmp(var,'s'))    %composite variables: entropy
            data = s;
            v_defs{vv} = 'Entropy';
            v_unitss{vv} = 'J/(kgK)';
            
        elseif(strcmp(var,'s_sat'))    %composite variables: saturation entropy
            data = s_sat;
            v_defs{vv} = 'Saturation entropy';
            v_unitss{vv} = 'J/(kgK)';
            
        elseif(strcmp(var,'thv'))    %composite variables: virtual potential temperature
            data = thv;
            v_defs{vv} = 'Theta-v';
            v_unitss{vv} = 'K';

        elseif(strcmp(var,'p'))    %composite variables: full pressure
            data = p;
            v_defs{vv} = 'Pressure';
            v_unitss{vv} = 'kg/(ms2)';

        elseif(strcmp(var,'th'))    %composite variables: potential temperature
            data = th;
            v_defs{vv} = 'Theta';
            v_unitss{vv} = 'K';
            
        elseif(strcmp(var,'qv'))    %composite variables: water vapor mixing ratio
            data = qv;
            v_defs{vv} = 'Water vapor mixing ratio';
            v_unitss{vv} = 'kg/kg';

        elseif(strcmp(var,'ql'))    %composite variables: liquid water mixing ratio
            data = ql;
            v_defs{vv} = 'Liquid water mixing ratio';
            v_unitss{vv} = 'kg/kg';
            
        elseif(strcmp(var,'qc'))    %composite variables: cloud water mixing ratio
            data = qc_pert;
            v_defs{vv} = 'Cloud water mixing ratio';
            v_unitss{vv} = 'kg/kg';
            
        elseif(strcmp(var,'qi'))    %composite variables: ice water mixing ratio
            data = qi_pert;
            v_defs{vv} = 'Ice water mixing ratio';
            v_unitss{vv} = 'kg/kg';
            
        elseif(strcmp(var,'qg'))    %composite variables: graupel water mixing ratio
            data = qg_pert;
            v_defs{vv} = 'Graupel water mixing ratio';
            v_unitss{vv} = 'kg/kg';

            
        else %just extract data for the input variable itself
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            v_defs{vv} = v_def;
            v_unitss{vv} = v_units;
            
            %data = s;
        end
%{
        %%STORE THE WIND SPEED AT THE DESIRED RADIUS
        if(r_ts==1)
            xres = dx;
            i_max = find(data==max(max(data)));  %index of max wind speed
            i_max = i_max(1);
            if(strmatch(v_usr,'MAX'))
                if(~isempty(i_max))
                    r_usr(ii)=xres*i_max-.5*xres;
                else
                    r_usr(ii)=NaN;
                end
            else
                if(strcmp(var,'vinterp') || strcmp(var,'vg'))
                    temp1=find(data<=v_usr);
                    temp2=find(temp1>i_max,1);
                    i_vusr_out=temp1(temp2);
                    clear temp1 temp2
                    if(~isempty(i_vusr_out))
                        v_out = data(i_vusr_out);   %v at first point where v<=v_usr
                        v_in = data(i_vusr_out-1);  %v at last point where v>=v_usr
                        r_usr(ii) = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                        if(r_usr(ii)<0 || r_usr(ii)>1000000 || max(max(data))<20)
                            r_usr(ii)=NaN;
                        end
                    else
                        v_out = NaN;
                        v_in = NaN;
                        r_usr(ii) = NaN;
                    end
                end
            end
        end
%}        
        
        %%CALCULATE TIME-AVERAGED xz-cross-section
        data=squeeze(data);
        data_tmean{vv}=data_tmean{vv}(1:nx_sub,1:ny_sub,1:nz_sub);
        data_tmean{vv}=data_tmean{vv}+data/(i_tf-i_t0+1);
        
        %%other time-averaged quantities (e.g. T rho rh the s s_sat gam_m thv)
%        rho_tmean=rho_tmean(1:nx_sub,1:ny_sub,1:nz_sub);
%        rho_tmean=rho_tmean+rho/(i_tf-i_t0+1);
        
    end

        
end



%% Extract subset of data for plotting
%ALL DATA
xvals = xmin_sub:dx:xmax_sub;
yvals = ymin_sub:dy:ymax_sub;
zvals = zmin_sub:dz:zmax_sub;
i_xvals = max(1,x0+1):max(1,x0+1)+nx_sub;
i_yvals = max(1,y0+1):max(1,y0+1)+ny_sub;
i_zvals = max(1,z0+1):max(1,z0+1)+nz_sub;

%SUBSET DATA
indices = find(xvals>=xmin_plot & xvals<=xmax_plot);
if(isempty(indices))
    indices = find(xvals>=xmin_plot,1);
end
xvals = xvals(indices);
i_xvals = i_xvals(indices);
clear indices
indices = find(yvals>=ymin_plot & yvals<=ymax_plot);
if(isempty(indices))
    indices = find(yvals>=ymin_plot,1);
end
yvals = yvals(indices);
i_yvals = i_yvals(indices);
clear indices
indices = find(zvals>=zmin_plot & zvals<=zmax_plot);
if(isempty(indices))
    indices = find(zvals>=zmin_plot,1);
end
zvals = zvals(indices);
i_zvals = i_zvals(indices);
clear indices

for vv=1:numvars
    data_tmean_sub{vv} = data_tmean{vv}(i_xvals',i_yvals',i_zvals');
    data_tmean_sub{vv}=double(data_tmean_sub{vv});  %pcolor requires double format
end

%% Spatial averaging?
if(ave_x==1)
    data_tmean_sub{vv} = nanmean(data_tmean_sub{vv},1);
    ave_x_str = sprintf('x-mean [%3.1f, %3.1f]',min(xvals),max(xvals));
    xvals = nanmean(xvals);
end
if(ave_y==1)
    data_tmean_sub{vv} = nanmean(data_tmean_sub{vv},2);
    ave_y_str = sprintf('y-mean [%3.1f, %3.1f]',min(yvals),max(yvals));
    yvals = nanmean(yvals);
end
if(ave_z==1)
    data_tmean_sub{vv} = nanmean(data_tmean_sub{vv},3);
    ave_z_str = sprintf('z-mean [%3.1f, %3.1f]',min(zvals),max(zvals));
    zvals = nanmean(zvals);
end

data_tmean_sub{vv}=squeeze(data_tmean_sub{vv});

%%BASIC STATISTICS
for vv=1:numvars
    data_min(vv) = min(min(data_tmean_sub{vv}));
    data_max(vv) = max(max(data_tmean_sub{vv}));
end

%%BLOCK OUT UNDESIRED DATA VALUES
if(max(max(max(data_tmean_sub{1})))>datamax_plot)
    sprintf('WARNING: MAX VALUE ABOVE THRESHOLD!')
end
if(min(min(min(data_tmean_sub{1})))<datamin_plot)
    sprintf('WARNING: MIN VALUE ABOVE THRESHOLD!')
end
data_tmean_sub{1}(data_tmean_sub{1}>datamax_plot)=datamax_plot;
data_tmean_sub{1}(data_tmean_sub{1}<datamin_plot)=datamin_plot;


%% PLOT THINGS
set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
h=figure(1);

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 15 15];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.15    0.15    0.80    0.70]);

plot_dims = [length(xvals) length(yvals) length(zvals)];
if(plot_dims(1)>1 && plot_dims(2)==1 && plot_dims(3)>1)   %xz (rz) color plot

    %turn vector data into matrices for plotting
    xmat = repmat(xvals,length(zvals),1);
    zmat = repmat(zvals,length(xvals),1)';

    %make contour color plot
    contourf(xmat,zmat,data_tmean_sub{1}',20)
    
    %colorbar
    temp=caxis;
    caxis(temp);    %are you serious?
    clear temp
    colorbar
    
    %%for plot labeling
    input_xlabel=sprintf('x-distance [%s]',xunits);
    input_ylabel=sprintf('z-distance [%s]',zunits);
    if(ave_y==1)
        input_title3=sprintf('%s %s',ave_y_str,yunits);
    else
        input_title3=sprintf('y=%5.2f %s',yvals,yunits);
    end

elseif(plot_dims(1)==1 && plot_dims(2)>1 && plot_dims(3)>1)   %yz color plot

    %turn vector data into matrices for plotting
    ymat = repmat(yvals,length(zvals),1);
    zmat = repmat(zvals,length(yvals),1)';

    %make contour color plot
    contourf(ymat,zmat,data_tmean_sub{1}',20)
    
    %colorbar
    temp=caxis;
    caxis(temp);    %are you serious?
    clear temp
    colorbar
    
    %%for plot labeling
    input_xlabel=sprintf('y-distance [%s]',yunits);
    input_ylabel=sprintf('z-distance [%s]',zunits);
    if(ave_x==1)
        input_title3=sprintf('%s %s',ave_x_str,xunits);
    else
        input_title3=sprintf('x=%5.2f %s',xvals,xunits);
    end

elseif(plot_dims(1)>1 && plot_dims(2)>1 && plot_dims(3)==1)   %xy color plot

    %turn vector data into matrices for plotting
    xmat = repmat(xvals,length(yvals),1);
    ymat = repmat(yvals,length(xvals),1)';

    %make contour color plot
    contourf(xmat,ymat,data_tmean_sub{1}',20)
    
    %colorbar
    temp=caxis;
    caxis(temp);    %are you serious?
    clear temp
    colorbar
    
    %%for plot labeling
    input_xlabel=sprintf('x-distance [%s]',xunits);
    input_ylabel=sprintf('y-distance [%s]',yunits);
    input_title3=sprintf('%s=%5.2f %s',ave_z_str,zvals,zunits);
    if(ave_z==1)
        input_title3=sprintf('%s %s',ave_z_str,zunits);
    else
        input_title3=sprintf('z=%5.2f %s',zvals,zunits);
    end

elseif(plot_dims(1)>1 && plot_dims(2)==1 && plot_dims(3)==1)   %x (r) only plot
    
    plot(xvals,data_tmean_sub{1})
    hold on
    x_zeroline = min(0,1.1*min(xvals)):(max(xvals)-min(xvals))/100:1.1*max(xvals);
    plot(x_zeroline,0*x_zeroline,'k')
    input_xlabel=sprintf('x-distance [%s]',xunits);
    input_ylabel=sprintf('%s [%s]',v_defs{1},v_unitss{1});
    axis([min(xvals) max(xvals) min(data_tmean_sub{1}) max(data_tmean_sub{1})])
    
    if(ave_y==1 && ave_z==1)
        input_title3=sprintf('%s %s; %s %s',ave_y_str,yunits,ave_z_str,zunits);
    elseif(ave_y==1)
        input_title3=sprintf('%s %s; z=%5.2f %s',ave_y_str,yunits,zvals,zunits);
    elseif(ave_z==1)
        input_title3=sprintf('y=%5.2f %s; %s %s',yvals,yunits,ave_z_str,zunits);
    else
        input_title3=sprintf('y=%5.2f %s; z=%5.2f %s',yvals,yunits,zvals,zunits);
    end

elseif(plot_dims(1)==1 && plot_dims(2)>1 && plot_dims(3)==1)   %y (r) only plot
    
    plot(yvals,data_tmean_sub{1})
    hold on
    y_zeroline = min(0,1.1*min(yvals)):(max(yvals)-min(yvals))/100:1.1*max(yvals);
    plot(y_zeroline,0*y_zeroline,'k')
    input_xlabel=sprintf('y-distance [%s]',xunits);
    input_ylabel=sprintf('%s [%s]',v_defs{1},v_unitss{1});
    axis([min(yvals) max(yvals) min(data_tmean_sub{1}) max(data_tmean_sub{1})])
    
    if(ave_x==1 && ave_z==1)
        input_title3=sprintf('%s %s; %s %s',ave_x_str,xunits,ave_z_str,zunits);
    elseif(ave_x==1)
        input_title3=sprintf('%s %s; z=%5.2f %s',ave_x_str,xunits,zvals,zunits);
    elseif(ave_z==1)
        input_title3=sprintf('x=%5.2f %s; %s %s',xvals,xunits,ave_z_str,zunits);
    else
        input_title3=sprintf('x=%5.2f %s; z=%5.2f %s',xvals,xunits,zvals,zunits);
    end

elseif(plot_dims(1)==1 && plot_dims(2)==1 && plot_dims(3)>1)   %z only plot
    
    plot(data_tmean_sub{1},zvals)
    hold on
    z_zeroline = min(0,1.1*min(zvals)):(max(zvals)-min(zvals))/100:1.1*max(zvals);
    plot(0*z_zeroline,z_zeroline,'k')
    input_xlabel=sprintf('x-distance [%s]',xunits); %%FIX ME
    input_ylabel=sprintf('%s [%s]',v_defs{1},v_unitss{1});
    axis([min(data_tmean_sub{1}) max(data_tmean_sub{1}) min(zvals) max(zvals)])
    
    if(ave_x==1 && ave_y==1)
        input_title3=sprintf('%s %s; %s %s',ave_x_str,xunits,ave_y_str,yunits);
    elseif(ave_x==1)
        input_title3=sprintf('%s %s; y=%5.2f %s',ave_x_str,xunits,yvals,yunits);
    elseif(ave_y==1)
        input_title3=sprintf('x=%5.2f %s; %s %s',xvals,xunits,ave_y_str,yunits);
    else
        input_title3=sprintf('x=%5.2f %s; y=%5.2f %s',xvals,xunits,yvals,yunits);
    end
    
else
    assert(2==3,'Your requested sub-domain doesnt make any sense!')
end

figure(1)
set(gca,'fontweight','bold','fontsize',11)
xlabel(input_xlabel);
ylabel(input_ylabel);
input_title1=sprintf('%s [%s]',v_defs{1},v_unitss{1});
input_title2=sprintf('min=%5.2f max=%5.2f',data_min(1),data_max(1));
title({input_title1,input_title2,input_title3})

%input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%legend(input_legend)

%% FIGURE TITLE
input_title=strrep(sprintf('time-mean days %i-%i',t0,tf),'_','\_');
annotation('textbox','String',input_title,'Position',[0 .5 .2 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')
fig_title = sprintf('%s',subdir);
annotation('textbox','String',fig_title,'Position',[0 -.45 .5 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')
