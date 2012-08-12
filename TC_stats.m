%TC_stats.m

%Created: 10-05-11, Dan Chavas
%Updated: 12-02-11, Dan Chavas

%Purpose: This file performs all the calculations to characterize the
%structure and evolution of the storms of interest.  The output is then
%saved in individual files (one per simulation) in directory simdata_Tmean#_#_#/

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

%%TESTING
%{
subdir_pre = 'CTRL_icRCE/';
ext_hd = 1;
run_type = 1;
t0 = 0;
tf = 150;
tmean0_usr = 100;
tmeanf_usr = 150;
v_usr_fracVp = 0.1000;
T_mean = 5;
dt_equil = 30;
dt_final = 50;
save_file = 0;
subdir = 'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh3000';
x0 = 0;
xf = 100000;
y0 = 0;
yf = 0;
z0 = 0;
zf = 1000;
rmin_sub = 0;
rmax_sub = 100000;
zmin_subsub = 0.5000;
zmax_subsub = 1;
%}

function [junk] = TC_stats(subdir_pre,ext_hd,run_type,t0,tf,tmean0_usr,tmeanf_usr,v_usr_fracVp,T_mean,dt_equil,dt_final,save_file,subdir,x0,xf,y0,yf,z0,zf,rmin_sub,rmax_sub,zmin_subsub,zmax_subsub,dir_home,moist);

%for calculating the outer radius using control values of the constant parameters
wrad_ctrl = .0027;   %control run value
Cd_in_ctrl = .0015;
fcor_ctrl = 5e-5;   %%THIS IS NO LONGER USED IN CALCULATION OF r0Lil_Lilctrl* !!

%%Write out to screen whats going on
sprintf('TC_stats for: %s',subdir)

junk='junk';
%clear
%clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
%subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
%ext_hd = 1; %0=local hard drive; 1=external hard drive

%run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

%tmean0_usr = 100;    %[day]
%tmeanf_usr = 150;    %[day]

%v_usr_fracVp = .1;  %wind speed as fraction of Vp; beyond this radius, radiative subsidence radial wind profile should apply.
%T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
%dt_equil = 30;  %[days]; how long must be quasi-steady to define equilibrium
%dt_final = 50;  %[days]; length of period from end of simulation over which equilibrium is calculated

%save_file = 1;

%%CONSTANTS
Cpd = 1004; %[J/K/kg]

file_in = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir);

if(exist(file_in)==2 && save_file == 1)
    sprintf('Data already saved for ax%s',subdir)
    
    %%Load data for given simulation to fix things within if needed
%    load(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdirs_save{rr}));

%{    
%    clear Vpots
%    subdirs_save = subdirs_save_drc;
%    clear subdirs_save_drc
%    save(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdirs_save{rr}));
    
    %Extract and save MPI in existing data file
%    mpi = mpi_retrieve(subdir) %[ms-1]
%    save(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir),'mpi','-append')
    
    %Extract and save lh from simulation namelist.input file
    %%EXTRACT RUN_TYPE AND SUBDIR NAME
    run_type=run_types(rr);

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
        run_type_str='3D';
        if(ext_hd==0)
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3D/CM1_output/%s',subdir_pre);
        else    %external harddrive
            dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/3D/%s',subdir_pre);
        end
        copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
    end

    %%DIRECTORY WITH OUTPUT DATA
    subdir_full=sprintf('%s%s',dir_in,subdir)
    
    [lh] = lh_retrieve(subdir_full)
    save(sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir),'lh','-append')
%}
    
    
else

%% EXTRACT MPI
mpi = mpi_retrieve(subdir); %[ms-1]
    
%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;

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
    run_type_str='3D';
    if(ext_hd==0)
        dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3D/CM1_output/%s',subdir_pre);
    else    %external harddrive
        dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/3D/%s',subdir_pre);
    end
    copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
end

%%DIRECTORY WITH OUTPUT DATA
subdir_full=sprintf('%s%s',dir_in,subdir);

%%EXTRACT lh
[lh] = lh_retrieve(subdir_full);

%%EXTRACT Cd
[Cd_in] = Cd_retrieve(subdir_full);

%%EXTRACT TIMESTEP SIZE
var_dt = 'vinterp'; %doesn't matter, just need any variable
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
if(run_type==3)
    numfiles=length(dir(sprintf('%s/cm1out_t*.nc',subdir_full)));
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_t0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
elseif(run_type==1)
    numfiles=length(dir(sprintf('%s/cm1out_*.nc',subdir_full)));
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
end

dt = time; %[s]; time of cm1out_t0001.nc is defined as zero
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub dx dy nx_sub ny_sub xunits yunits zunits v_def v_units time t_units
i_t0 = max(0,round(t0*24*60*60/dt)+1); %start time for analysis
i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %end time for analysis
i_t0_usr = max(0,round(tmean0_usr*24*60*60/dt)+1); %timestep corresponding to start time for averaging
i_tf_usr = min(numfiles,round(tmeanf_usr*24*60*60/dt)+1); %timestep corresponding to end time for averaging

%%TIME IN DAYS
%    t_min_day = dt*i_t0/86400;
%    t_max_day = dt*i_tf/86400;

%%DIRECTORY WITH OUTPUT DATA
subdir_full=sprintf('%s%s',dir_in,subdir);

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding';

%%load initial sounding variables
dir_full = strcat(dir_in,subdir);

[zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

%%Estimate w_rad (clear-sky subsidence rate) from initial RCE sounding
th_lowtrop = th00(zz00<=5000 & zz00>=1500);
zz_lowtrop = zz00(zz00<=5000 & zz00>=1500);
dthdz_wrad = (th_lowtrop(end)-th_lowtrop(1))/(zz_lowtrop(end)-zz_lowtrop(1));

%Extract Qrad [K/day]
temp1 = subdir(strfind(subdir,'rad')+3:end);
if(isempty(temp1))
    Qrad = 1;
else
    itemp = strfind(temp1,'K');
    itemp = itemp(1);
    Qrad = temp1(1:itemp-1); %K/halfday
    Qrad = 2*str2num(Qrad); %K/day
end
wrad = (Qrad/86400) / dthdz_wrad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FULL WIND %%%%%%%
%{
var = 'vinterp';    %'vinterp'=full wind; 'vg'=gradient wind

%%initialize the output vectors
clear v_df_qv v_units_qv data_tmean_usr data_tmean
data_tmean_usr=zeros(xf-x0+1,1);
data_tmean=zeros(xf-x0+1,tf-t0);    %1 profile per day
Vmax_movave=nan(i_tf-i_t0+1,1);
rmax_movave=nan(i_tf-i_t0+1,1);
rrad_movave=nan(i_tf-i_t0+1,1);
r0_movave=nan(i_tf-i_t0+1,1);
r0Lil_movave=nan(i_tf-i_t0+1,1);
num_prof = 1;   %initialize iteration of time average profiles
numerr(rr)=0;

clear t_day
clear rrad r0
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

        cd(dir_home)

    end

    %% Calculate various thermodynamic quantities: NOT NEEDED FOR WIND STATS :)
%{        
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

    %%load vpert
    var_temp = 'vinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    v_pert = squeeze(data);

    %%calculate pressure
    p = repmat(pp00(z0+1:z0+nz_sub),nx_sub,1) + pp_pert;

    %%calculate potential temperature
    th = repmat(th00(z0+1:z0+nz_sub),nx_sub,1) + th_pert;

    %%calculate water vapor mixing ratio
    qv = repmat(qv00(z0+1:z0+nz_sub),nx_sub,1) + qv_pert;    %[kg/kg]

    %%calculate T, Tv, rho, rh
    [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);
%}        

    if(strcmp(var,'vg'))    %gradient wind

        %%load prspert
        var_temp = 'prspert';
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
        pp_pert = squeeze(data);

        %%calculate pressure
        p = repmat(pp00(z0+1:z0+nz_sub),nx_sub,1) + pp_pert;

        %%EXTRACT DATA: just need the model parameters here actually
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);

        xvals = xmin_sub:dx:xmax_sub;
        zvals = zmin_sub:dz:zmax_sub;
        xmat = repmat(xvals,length(zvals),1)';

        %%Calculate radial pressure gradients at grid-points
        dp_dr = (p(3:end,:)-p(1:end-2,:))/(2*1000*dx); %[m/s^2]
        dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value
        dp_dr(end+1,:)=dp_dr(end,:);  %add extra level at outer edge equal to nearest level value

        %%alternative: calculate at mid-points between grid points
%            dp_dr = (p(2:end,:)-p(1:end-1,:))/(1000*dx); %[m/s^2]
%            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value

        %%Calculate gradient wind from pressure field
        vinterp = data;
        temp2=(xmat*1000*fcor).^2+4*xmat*1000.*dp_dr;
        data =  .5*(-xmat*1000*fcor+sqrt((xmat*1000*fcor).^2+4*xmat*1000.*dp_dr));   %Vg^2 + r*f*Vg-r*(dp_dr); solve for Vg
        data = real(data);  %remove imaginary numbers, which are always very small
        %data(imag(data)~=0)=NaN;    %ignore imaginary numbers
%            data = data - squeeze(vinterp);

        v_def = 'Gradient wind';
        v_units = 'm/s';

    elseif(strcmp(var,'vinterp'))   %full azimuthal wind

        %%EXTRACT DATA
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

    end

    %% Find quasi-steady states

    %%Extract subset of data corresponding to desired wind level
    %ALL DATA
    xvals = xmin_sub:dx:xmax_sub;
    zvals = zmin_sub:dz:zmax_sub;
    i_xvals = max(1,x0+1):max(1,x0+1)+nx_sub;
    i_zvals = max(1,z0+1):max(1,z0+1)+nz_sub;

    %SUBSET DATA
    indices = find(xvals>=rmin_sub & xvals<=rmax_sub);
    i_xvals = i_xvals(indices);
    clear indices
    indices = find(zvals>=zmin_subsub & zvals<=zmax_subsub,1);
    indices = indices(end); %pick the highest level below z=1km
    i_zvals = i_zvals(indices);
    clear indices

    data_temp = squeeze(data);
    data_sub = data_temp(i_xvals',i_zvals');
    xvals_sub = xvals(i_xvals);
    data_tmean=data_tmean(i_xvals,:);

    %%Calculate and save instantaneous Vmax time-series
    i_temp = round(length(data_sub)/8); %should be at small r relative to domain size
    Vmax(ii) = max(data_sub(1:i_temp));  %index of max wind speed

    %%Calculate and save instantaneous rmax time-series
    xres = dx;
    i_max = find(data_sub(1:i_temp)==max(data_sub(1:i_temp)));  %index of max wind speed
    i_max = i_max(1);
    if(~isempty(i_max))
        rmax(ii)=xres*i_max-.5*xres;
    else
        rmax(ii)=NaN;
    end

    %%Calculate and save instantaneous rrad time-series
    v_usr = v_usr_fracVp*mpi;
    data_sub_temp = smooth(data_sub,30);
    temp1=find(data_sub_temp<=v_usr);
    temp2=find(temp1>i_max,1);
    i_vusr_out=temp1(temp2);
    clear temp1 temp2
    if(~isempty(i_vusr_out))
        v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
        v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
        r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
        if(r_usr<0 || r_usr>1000000 || max(data_sub)<15)
            r_usr=NaN;
        end
    else
        r_usr = NaN;
    end
    rrad(ii)=r_usr;
    clear i_vusr_out v_out v_in

    %%Calculate and save instantaneous r0 time-series
    v_usr = 0;
    data_sub_temp = smooth(data_sub,30);
    temp1=find(data_sub_temp<=v_usr);
    temp2=find(temp1>i_max,1);
    i_vusr_out=temp1(temp2);
    clear temp1 temp2
    if(~isempty(i_vusr_out))
        v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
        v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
        r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
        if(r_usr<0 || r_usr>1000000 || max(data_sub)<15)
            r_usr=NaN;
        end
    else
        r_usr = NaN;
    end
    r0(ii)=r_usr;
    clear i_vusr_out v_out v_in i_max

    %%Calculate T-day running mean time-series of rmax and r0
    nfile_mean = T_mean*(24*60*60/dt)+1;  %corresponding number of files

    if(ii==1)
        v_r_store = data_sub;
        v_r_mean = 0*data_sub; %initialize moving average radial profile
    else
        v_r_store = [v_r_store data_sub];
    end

    if(ii>=nfile_mean)   %if have enough data to cover T_mean period, then begin to save
        v_r_mean = mean(v_r_store,2);   %represents the middle of the period
        v_r_store = v_r_store(:,2:end); %drop the earliest radial wind profile to make room for the next iteration

        %%Calculate and save T-day running mean Vmax
        Vmax_movave_temp = max(v_r_mean);  %index of max wind speed
        Vmax_movave(ii-floor(nfile_mean/2))=Vmax_movave_temp;

        %%Calculate and save T-day running mean rmax
        i_max = find(v_r_mean(1:i_temp)==max(v_r_mean(1:i_temp)));  %index of max wind speed
        i_max = i_max(1);
        if(~isempty(i_max))
            rmax_movave_temp=xres*i_max-.5*xres;
        else
            rmax_movave_temp=NaN;
        end
        rmax_movave(ii-floor(nfile_mean/2))=rmax_movave_temp;

        %%Calculate and save T-day running mean rrad
        v_usr = v_usr_fracVp*mpi;
        v_r_mean_temp = smooth(v_r_mean,30);
        temp1=find(v_r_mean_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<15)
                r_usr_movave=NaN;
            end
        else
            r_usr_movave = NaN;
        end
        rrad_movave(ii-floor(nfile_mean/2)) = r_usr_movave;
        if(~isnan(r_usr_movave))
            V_user = v_usr_fracVp*mpi;
            r_user = rrad_movave(ii-floor(nfile_mean/2));
            [r_0_full,res_error,i_error] = r0_calc(r_user,V_user,fcor,Cd_in,wrad);
            r_0_full/1000
            numerr = numerr+i_error;
            clear junk1 junk2

            r0Lil_movave(ii-floor(nfile_mean/2)) = r_0_full/1000;
        else
            r0Lil_movave(ii-floor(nfile_mean/2)) = NaN;
        end
        clear i_vusr_out v_out v_in

        %%Calculate and save T-day running mean r0
        v_usr = 0;
        v_r_mean_temp = smooth(v_r_mean,30);
        temp1=find(v_r_mean_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<15)
                r_usr_movave=NaN;
            end
        else
            r_usr_movave = NaN;
        end
        r0_movave(ii-floor(nfile_mean/2)) = r_usr_movave;
        clear i_vusr_out v_out v_in i_max
    end

        
    %%CALCULATE days tmean0-tmeanf TIME-AVERAGED xz-cross-section
    if(t_day(ii)>=tmean0_usr & t_day(ii)<=tmeanf_usr)
        data_tmean_usr=data_tmean_usr(i_xvals,1);
        data_tmean_usr=data_tmean_usr+data_sub/(i_tf_usr-i_t0_usr+1);
    end

    %%CALCULATE TIME-AVERAGED xz-cross-section for set time periods
    data_tmean(:,num_prof)=data_tmean(:,num_prof)+data_sub/(T_mean*(60*60*24/dt));
    if(ii>1)
        if(mod(t_day(ii),1)<mod(t_day(ii-1),1))
            num_prof=num_prof+1;
        end
    end

end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GRADIENT WIND %%%%%%%
var = 'vg';    %'vinterp'=full wind; 'vg'=gradient wind

%%initialize the output vectors
clear v_df_qv v_units_qv data_tmean_usr_g data_tmean_g
Vmax_movave_g=nan(i_tf-i_t0+1,1);
rmax_movave_g=nan(i_tf-i_t0+1,1);
rrad_movave_g=nan(i_tf-i_t0+1,1);
r0_movave_g=nan(i_tf-i_t0+1,1);
r0Lil_movave_g=nan(i_tf-i_t0+1,1);
r0Lil_Lilctrl_movave_g=nan(i_tf-i_t0+1,1);
r0ER11_movave_g=nan(i_tf-i_t0+1,1);
num_prof = 1;   %initialize iteration of time average profiles
numerr=0;

clear t_day
clear rrad r0
for ii=1:i_tf-i_t0+1
    ii
    fclose('all');
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
        dir_tmp = strcat(dir_in,subdir);
        cd(dir_tmp);
        nc_file;
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

        cd(dir_home);
        
        data_tmean_usr_g=zeros(nx,1);
        data_tmean_g=zeros(nx,tf-t0);    %1 profile per day
        v_r_all = NaN*zeros(nx,i_tf-i_t0+1);

    end

    %% Calculate various thermodynamic quantities: NOT NEEDED FOR WIND STATS :)      
    %%load prspert
    var_temp = 'prspert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    pp_pert = squeeze(data);

    %%load thpert
    var_temp = 'thpert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    th_pert = squeeze(data);

    if(moist == 1)
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
    else    %DRY
        qv_pert = 0;
        qc_pert = 0;
        qr_pert = 0;
        qi_pert = 0;
        qs_pert = 0;
        qg_pert = 0;
    end
    
    %%ql_pert
    ql_pert = qc_pert + qr_pert + qi_pert + qs_pert + qg_pert;
    
    %%load upert
    var_temp = 'uinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    u_pert = squeeze(data);

    %%load vpert
    var_temp = 'vinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    v_pert = squeeze(data);

    %%calculate pressure
    p = repmat(pp00(z0+1:z0+nz_sub),nx_sub,1) + pp_pert;

    %%calculate potential temperature
    th = repmat(th00(z0+1:z0+nz_sub),nx_sub,1) + th_pert;

    %%calculate water vapor mixing ratio
    if(moist == 1)
        qv = repmat(qv00(z0+1:z0+nz_sub),nx_sub,1) + qv_pert;    %[kg/kg]

        %%calculate liquid water mixing ratio
        ql = ql_pert;    %[kg/kg]

    else    %DRY
        qv = 0;
        ql = 0;
    end

    %%calculate T, Tv, rho, rh
    [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql);

    if(strcmp(var,'vg'))    %gradient wind

        var_temp = 'pipert';
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
        pi_pert = squeeze(data);

        %%calculate pi
        pi = repmat(pi00(z0+1:z0+nz_sub),nx_sub,1) + pi_pert;

        %%EXTRACT DATA: just need the model parameters here actually
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);

        xvals = xmin_sub:dx:xmax_sub;
        zvals = zmin_sub:dz:zmax_sub;
        xmat = repmat(xvals,length(zvals),1)';

        %%Calculate radial pi gradients at grid-points
        dpi_dr = (pi(3:end,:)-pi(1:end-2,:))/(2*1000*dx); %[m/s^2]
        dpi_dr(2:end+1,:)=dpi_dr;  %add extra level at center equal to nearest level value
        dpi_dr(end+1,:)=dpi_dr(end,:);  %add extra level at outer edge equal to nearest level value

        %%alternative: calculate at mid-points between grid points
%            dp_dr = (p(2:end,:)-p(1:end-1,:))/(1000*dx); %[m/s^2]
%            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value

        %%Calculate gradient wind from pi field (BR09b p.3)
        data =  .5*(-xmat*1000*fcor+sqrt((xmat*1000*fcor).^2+4*xmat*1000*Cpd.*thv.*dpi_dr));   %Vg^2 + r*f*Vg-r*(dp_dr); solve for Vg
        data = real(data);  %remove imaginary numbers, which are always very small
        %data(imag(data)~=0)=NaN;    %ignore imaginary numbers
%            data = data - squeeze(vinterp);

        v_def = 'Gradient wind';
        v_units = 'm/s';

    elseif(strcmp(var,'vinterp'))   %full azimuthal wind

        %%EXTRACT DATA
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

    end

    %% Find quasi-steady states

    %%Extract subset of data corresponding to desired wind level
    %ALL DATA
    xvals = xmin_sub:dx:xmax_sub;
    zvals = zmin_sub:dz:zmax_sub;
    i_xvals = max(1,x0+1):max(1,x0+1)+nx_sub;
    i_zvals = max(1,z0+1):max(1,z0+1)+nz_sub;

    %SUBSET DATA
    indices = find(xvals>=rmin_sub & xvals<=rmax_sub);
    i_xvals = i_xvals(indices);
    clear indices
    indices = find(zvals>=zmin_subsub & zvals<=zmax_subsub);
    indices = indices(end); %pick the highest level below z=1km
    i_zvals = i_zvals(indices);
    clear indices

    data_temp = squeeze(data);
    data_sub = data_temp(i_xvals',i_zvals');
    xvals_sub = xvals(i_xvals);
    data_tmean_g=data_tmean_g(i_xvals,:);

    %% Calculate and save instantaneous Vmax time-series
    i_temp = round(length(data_sub)/8); %CHECKS ONLY FIRST 1/8 OF DOMAIN
    Vmax_g(ii) = max(data_sub(1:i_temp));  %index of max wind speed

    %% Calculate and save instantaneous rmax time-series
    xres = dx;
    i_max = find(data_sub(1:i_temp)==max(data_sub(1:i_temp)));  %index of max wind speed
    i_max = i_max(1);
    if(~isempty(i_max))
        rmax_g(ii)=xres*i_max-.5*xres;
    else
        rmax_g(ii)=NaN;
    end

    %% Calculate and save instantaneous r[v_usr] time-series
    v_usr = v_usr_fracVp*mpi;
    Vrad = v_usr;
    data_sub_temp=smooth(data_sub,10);  %use 10-point smoother on profile to calculate
    temp1=find(data_sub_temp<=v_usr);
    temp2=find(temp1>i_max,1);
    i_vusr_out=temp1(temp2);
    clear temp1 temp2
    if(~isempty(i_vusr_out))
        v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
        v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
        r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
        if(r_usr<0 || r_usr>1000000 || max(data_sub)<Vrad)
            r_usr=NaN;
        end
    else
        r_usr = NaN;
    end
    rrad_g(ii)=r_usr;
    clear i_vusr_out v_out v_in
 
    %% Calculate and save instantaneous r0 time-series
    v_usr = 0;
    data_sub_temp=smooth(data_sub,10);  %use 10-point smoother on profile to calculate
    temp1=find(data_sub_temp<=v_usr);
    temp2=find(temp1>i_max,1);
    i_vusr_out=temp1(temp2);
    clear temp1 temp2
    if(~isempty(i_vusr_out))
        v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
        v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
        r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
        if(r_usr<0 || r_usr>1000000 || max(data_sub)<Vrad)
            r_usr=NaN;
        end
    else
        r_usr = NaN;
    end
    r0_g(ii)=r_usr;
    clear i_vusr_out v_out v_in i_max

    %% Calculate T-day running mean time-series of same variables
    nfile_mean = T_mean*(24*60*60/dt)+1;  %corresponding number of files

    if(ii==1)
        v_r_store = data_sub;
        v_r_mean = 0*data_sub; %initialize moving average radial profile
    else
        v_r_store = [v_r_store data_sub];
    end

    if(ii>=nfile_mean)   %if have enough data to cover T_mean period, then begin to save
        v_r_mean = mean(v_r_store,2);   %represents the middle of the period
        v_r_store = v_r_store(:,2:end); %drop the earliest radial wind profile to make room for the next iteration

        %%Calculate and save T-day running mean Vmax
        Vmax_movave_g_temp = max(v_r_mean);  %index of max wind speed
        Vmax_movave_g(ii-floor(nfile_mean/2))=Vmax_movave_g_temp;

        %%Calculate and save T-day running mean rmax
        i_max = find(v_r_mean(1:i_temp)==max(v_r_mean(1:i_temp)));  %index of max wind speed
        i_max = i_max(1);
        if(~isempty(i_max))
            rmax_movave_g_temp=xres*i_max-.5*xres;
        else
            rmax_movave_g_temp=NaN;
        end
        rmax_movave_g(ii-floor(nfile_mean/2))=rmax_movave_g_temp;

        %%Calculate and save T-day running mean rrad, and calculate r0Lil from rrad
        v_usr = v_usr_fracVp*mpi;
        v_r_mean_temp=smooth(v_r_mean,10);  %use 10-point smoother on profile to calculate
        temp1=find(v_r_mean_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<Vrad)
                r_usr_movave=NaN;
            end
        else
            r_usr_movave = NaN;
        end

        rrad_movave_g(ii-floor(nfile_mean/2)) = r_usr_movave;
        if(~isnan(r_usr_movave))
            V_user = v_usr_fracVp*mpi;
            r_user = rrad_movave_g(ii-floor(nfile_mean/2));
            [r_0_full,res_error,i_error] = r0_calc(r_user,V_user,fcor,Cd_in,wrad);
            numerr = numerr+i_error;
            [r_0_full_Lilctrl,res_error,i_error] = r0_calc(r_user,V_user,fcor,Cd_in_ctrl,wrad_ctrl);
            [r0ER11] = r0ER11_calc(r_user,v_usr_fracVp,mpi,fcor,Cd_in);    %calculate outer radius using ER11 Eq 36
            clear junk1 junk2

            r0Lil_movave_g(ii-floor(nfile_mean/2)) = r_0_full/1000;
            r0Lil_Lilctrl_movave_g(ii-floor(nfile_mean/2)) = r_0_full_Lilctrl/1000;
            r0ER11_movave_g(ii-floor(nfile_mean/2)) = r0ER11/1000;
        else
            r0Lil_movave_g(ii-floor(nfile_mean/2)) = NaN;
            r0Lil_Lilctrl_movave_g(ii-floor(nfile_mean/2)) = NaN;
            r0ER11_movave_g(ii-floor(nfile_mean/2)) = NaN;
        end
        clear i_vusr_out v_out v_in

        %%Calculate and save T-day running mean r0
        v_usr = 0;
        v_r_mean_temp=smooth(v_r_mean,10);  %use 10-point smoother on profile to calculate
        temp1=find(v_r_mean_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<Vrad)
                r_usr_movave=NaN;
            end
        else
            r_usr_movave = NaN;
        end
        r0_movave_g(ii-floor(nfile_mean/2)) = r_usr_movave;
        clear i_vusr_out v_out v_in i_max
    end

%}  
    %%Save all radial profile data (i.e. with no averaging) -- need this to
    %%calculate mean radial profile for dynamic equilibrium
    v_r_all(1:length(data_sub),ii) = data_sub;
    
    %%CALCULATE days tmean0-tmeanf TIME-AVERAGED xz-cross-section
    if(t_day(ii)>=tmean0_usr & t_day(ii)<=tmeanf_usr)
        data_tmean_usr_g=data_tmean_usr_g(i_xvals,1);
        data_tmean_usr_g=data_tmean_usr_g+data_sub/(i_tf_usr-i_t0_usr+1);
    end

    %%CALCULATE TIME-AVERAGED xz-cross-section for set time periods
    data_tmean_g(:,num_prof)=data_tmean_g(:,num_prof)+data_sub/(T_mean*(60*60*24/dt));
    if(ii>1)
        if(mod(t_day(ii),1)<mod(t_day(ii-1),1))
            num_prof=num_prof+1;
        end
    end

end

%% Transient statistics %%%%%%%%%%
%%Find time to genesis, tau_gen, defined as the first time vmax_g  > vmax_equil_g

%Calculate mean Vmax at end of simulation
tdat_end = t_day(end-dt_final/(dt/60/60/24):end);
Vmax_dat_end = Vmax_movave_g(end-dt_final/(dt/60/60/24):end);
indices = ~isnan(Vmax_dat_end);
tdat_end = tdat_end(indices);
Vmax_dat_end = Vmax_dat_end(indices);
clear indices
Vmax_mean_end = mean(Vmax_dat_end);

%time to genesis and stats at genesis
indices = find(Vmax_movave_g>Vmax_mean_end);
i_gen = indices(1);
tau_gen_sim = t_day(i_gen); %defined using the GRADIENT wind
%Vmax_gen_sim = Vmax_movave(i_gen);
%rmax_gen_sim = rmax_movave(i_gen);
%rrad_gen_sim = rrad_movave(i_gen);
%r0_gen_sim = r0_movave(i_gen);
%r0Lil_gen_sim = r0Lil_movave(i_gen);
Vmax_gen_g_sim = Vmax_movave_g(i_gen);
rmax_gen_g_sim = rmax_movave_g(i_gen);
rrad_gen_g_sim = rrad_movave_g(i_gen);
r0_gen_g_sim = r0_movave_g(i_gen);
r0Lil_gen_g_sim = r0Lil_movave_g(i_gen);
r0Lil_Lilctrl_gen_g_sim = r0Lil_Lilctrl_movave_g(i_gen);
r0ER11_gen_g_sim = r0ER11_movave_g(i_gen);

%time to peak value and value
indices = find(Vmax_movave_g(i_gen:end)==max(Vmax_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
Vmax_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
Vmax_max_g_sim = Vmax_movave_g(i_peak);

indices = find(rmax_movave_g(i_gen:end)==max(rmax_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
rmax_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
rmax_max_g_sim = rmax_movave_g(i_peak);

indices = find(rrad_movave_g(i_gen:end)==max(rrad_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
rrad_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
rrad_max_g_sim = rrad_movave_g(i_peak);

indices = find(r0_movave_g(i_gen:end)==max(r0_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
r0_max_g_sim = r0_movave_g(i_peak);

indices = find(r0Lil_movave_g(i_gen:end)==max(r0Lil_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0Lil_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
r0Lil_max_g_sim = r0Lil_movave_g(i_peak);

indices = find(r0Lil_Lilctrl_movave_g(i_gen:end)==max(r0Lil_Lilctrl_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0Lil_Lilctrl_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
r0Lil_Lilctrl_max_g_sim = r0Lil_Lilctrl_movave_g(i_peak);

indices = find(r0ER11_movave_g(i_gen:end)==max(r0ER11_movave_g(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0ER11_tau_max_g_sim = t_day(i_peak); %defined using the GRADIENT wind
r0ER11_max_g_sim = r0ER11_movave_g(i_peak);

%{
indices = find(Vmax_movave(i_gen:end)==max(Vmax_movave(i_gen:end)));
i_peak = i_gen-1+indices(1);
Vmax_tau_max_sim = t_day(i_peak); %defined using the GRADIENT wind
Vmax_max_sim = Vmax_movave(i_peak);

indices = find(rmax_movave(i_gen:end)==max(rmax_movave(i_gen:end)));
i_peak = i_gen-1+indices(1);
rmax_tau_max_sim = t_day(i_peak); %defined using the GRADIENT wind
rmax_max_sim = rmax_movave(i_peak);

indices = find(rrad_movave(i_gen:end)==max(rrad_movave(i_gen:end)));
i_peak = i_gen-1+indices(1);
rrad_tau_max_sim = t_day(i_peak); %defined using the GRADIENT wind
rrad_max_sim = rrad_movave(i_peak);

indices = find(r0_movave(i_gen:end)==max(r0_movave(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0_tau_max_sim = t_day(i_peak); %defined using the GRADIENT wind
r0_max_sim = r0_movave(i_peak);

indices = find(r0Lil_movave(i_gen:end)==max(r0Lil_movave(i_gen:end)));
i_peak = i_gen-1+indices(1);
r0Lil_tau_max_sim = t_day(i_peak); %defined using the GRADIENT wind
r0Lil_max_sim = r0Lil_movave(i_peak);
%}

numvars = 12;     %Vmax, rmax, rrad, r0, r0Lil, Vmax_g, rmax_g, rrad_g, r0_g, r0Lil_g
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
        case 7
            var_movave = rmax_movave_g;
        case 8
            var_movave = rrad_movave_g;
        case 9
            var_movave = r0_movave_g;
        case 10
            var_movave = r0Lil_movave_g;
        case 11
            var_movave = r0Lil_Lilctrl_movave_g;
        case 12
            var_movave = r0ER11_movave_g;
    end

    %%Extract data for (dt_final)-day period at end of simulation and calculate its mean
    tdat_end = t_day(end-dt_final/(dt/60/60/24):end);
    var_dat_end = var_movave(end-dt_final/(dt/60/60/24):end);
    indices = ~isnan(var_dat_end);
    tdat_end = tdat_end(indices);
    var_dat_end = var_dat_end(indices);
    clear indices
    var_mean_end = mean(var_dat_end);

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

%        else
%            var_tau_equil(l) = NaN;

%        end
end

%% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
%variable values
Vmax_equil_sim = var_equil(1);
rmax_equil_sim = var_equil(2);
rrad_equil_sim = var_equil(3);
r0_equil_sim = var_equil(4);
r0Lil_equil_sim = var_equil(5);
Vmax_equil_g_sim = var_equil(6);
rmax_equil_g_sim = var_equil(7);
rrad_equil_g_sim = var_equil(8);
r0_equil_g_sim = var_equil(9);
r0Lil_equil_g_sim = var_equil(10);
r0Lil_Lilctrl_equil_g_sim = var_equil(11);
r0ER11_equil_g_sim = var_equil(12);

%timescales to those values
Vmax_tau_equil_sim = var_tau_equil(1);
rmax_tau_equil_sim = var_tau_equil(2);
rrad_tau_equil_sim = var_tau_equil(3);
r0_tau_equil_sim = var_tau_equil(4);
r0Lil_tau_equil_sim = var_tau_equil(5);
Vmax_tau_equil_g_sim = var_tau_equil(6);
rmax_tau_equil_g_sim = var_tau_equil(7);
rrad_tau_equil_g_sim = var_tau_equil(8);
r0_tau_equil_g_sim = var_tau_equil(9);
r0Lil_tau_equil_g_sim = var_tau_equil(10);
r0Lil_Lilctrl_tau_equil_g_sim = var_tau_equil(11);
r0ER11_tau_equil_g_sim = var_tau_equil(12);

%% Save data for all simulations
%Vmax_movave_sim=Vmax_movave;
%rmax_movave_sim=rmax_movave;
%rrad_movave_sim=rrad_movave;
%r0_movave_sim=r0_movave;
%r0Lil_movave_sim=r0Lil_movave;

Vmax_movave_g_sim=Vmax_movave_g;
rmax_movave_g_sim=rmax_movave_g;
rrad_movave_g_sim=rrad_movave_g;
r0_movave_g_sim=r0_movave_g;
r0Lil_movave_g_sim=r0Lil_movave_g;
r0Lil_Lilctrl_movave_g_sim=r0Lil_Lilctrl_movave_g;
r0ER11_movave_g_sim=r0ER11_movave_g;

%% User profile %%%%%%%%%%%%%%%%%%%%%%%%
xvals_sub_sim = xvals_sub;
%data_tmean_usr_sim = data_tmean_usr;
data_tmean_usr_g_sim = data_tmean_usr_g;

%% Save data to file
if(save_file == 1)
    
    save tempstuff.mat l ii
    clear l i ii
    save temp.mat
    load tempstuff.mat
    movefile('temp.mat',sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir))
    
    delete('tempstuff.mat')

end

end

end
