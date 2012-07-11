%TC_stats.m

%Created: 10-05-11, Dan Chavas

%Purpose: Plot sensitivity of storm structure to individual scale changes
%1) 

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc



%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D
%dz_frac=1;  %ratio of dz_model to dz_sounding

tmean0_usr = 70;    %[day]
tmeanf_usr = 100;    %[day]
pwd

T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
dt_equil = 30;  %[days]; how long must be quasi-steady to define equilibrium
Cd_in = 1.5e-3; %only used to calculate r0_Lilly

%% Plots on/off %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_gradient = 1;  %1=make plots with gradient wind
plot_full = 1;  %1=make plots with full wind (note: can do both)

sim_set = 'mpi';  %if file exists, will load automatically: domain; cor; qro; also the name of the output directory for plots

%SINGLE SIMULATION
plot_ts = 0;    %0=no; plot time series of vmax, rmax, r12, r0, r0Lil
plot_proftevol = 0; %0=no; plot radial wind profile at multiple times

%SINGLE/MULTI SIMULATION
plot_usrprof = 1;   %0=no; plot radial wind profile for user-defined time-means

%MULTI SIMULATION
plot_ts_multi = 1;    %0=no; plot time series of vmax, rmax, r12, r0
plot_stats = 1; %0=no; plot comparisons of basic statistics: Vmax, Rmax, R0, tau_equil for each; can be done at "equil time", genesis time, time of first V~Vmax, final (day 70-100 ave)
    save_file = 1;  %0=don't save; 1=save file with file name [out_file]
    CTRL_val = 93;   %CTRL value of quantity varied across simulations
    units = 'ms-1';
    multipliers = log2([40 49 80 93 120 125]/CTRL_val);  %log2(1.5) = x1.5
    %multipliers = [0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
   
%Plotting domain
rmin_plot = 0;  %[km]
rmax_plot = 2500;    %[km]
%zmin_plot = .5;
%zmax_plot = 1;
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

%%Simulations
subdirs = {


%'CTRLv0qro100000qrhSATqdz5000_nx3072'
%'CTRLv0qrhSATqdz5000_nx3072'
%'CTRLv0qro400000qrhSATqdz5000_nx3072'
%'CTRLv12.5qrh0qdz5000_nx3072'

%{
%%NON-DIM NUMBER SCALING Vpot/(f*lh)
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'
'CTRLv0qrhSATqdz5000_nx3072_fx2'
'CTRLv0qrhSATqdz5000_nx3072_lh3000'

'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2'

'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_lh750'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh750'
%}


%{
%%DOMAIN SIZE
'CTRLv0qrhSATqdz5000_nx192'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_nx768'
'CTRLv0qrhSATqdz5000_nx1536'
'CTRLv0qrhSATqdz5000_nx3072'
%}

%{
%%HORIZONTAL RESOLUTION
%'CTRLv0qrhSATqdz5000_nx12288_dx1'
'CTRLv0qrhSATqdz5000_nx6144_dx2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx1536_dx8'
'CTRLv0qrhSATqdz5000_nx768_dx16'
'CTRLv0qrhSATqdz5000_nx384_dx32'
%%'CTRLv0qrhSATqdz5000_nx6144_dx1'
%%'CTRLv0qrhSATqdz5000_nx1536_dx4'
%%'CTRLv0qrhSATqdz5000_nx3072_dx2'
%}

%{
%%VERTICAL RESOLUTION
%'CTRLv0qrhSATqdz5000_nx3072_dz156.25'
'CTRLv0qrhSATqdz5000_nx3072_dz312.5'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_dz1250'
%%'CTRLv0qrhSATqdz5000_nx1536_dx4_dz156.25'
%%'CTRLv0qrhSATqdz5000_nx1536_dx4_dz312.5'
%}

%{
%%HORIZONTAL MIXING LENGTH
'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
'CTRLv0qrhSATqdz5000_nx3072_lh375'
'CTRLv0qrhSATqdz5000_nx3072_lh750'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_lh6000'
'CTRLv0qrhSATqdz5000_nx3072_lh12000'
%}

%{
%%HORIZ SCALE OF MOISTURE PERTURBATION
'CTRLv0qro25000qrhSATqdz5000_nx3072'
'CTRLv0qro50000qrhSATqdz5000_nx3072'
'CTRLv0qro100000qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qro400000qrhSATqdz5000_nx3072'
%'CTRLv0qro500000qrhSATqdz5000_nx3072'
%'CTRLv0qro600000qrhSATqdz5000_nx3072'
%'CTRLv0qro700000qrhSATqdz5000_nx3072'
'CTRLv0qro800000qrhSATqdz5000_nx3072'
'CTRLv0qro1600000qrhSATqdz5000_nx3072'
%}

%{
%%HORIZ SCALE OF WIND PERTURBATION
'CTRLro50000v12.5qrh0qdz5000_nx3072'
'CTRLro100000v12.5qrh0qdz5000_nx3072'
'CTRLro200000v12.5qrh0qdz5000_nx3072'
'CTRLv12.5qrh0qdz5000_nx3072'
'CTRLro800000v12.5qrh0qdz5000_nx3072'
'CTRLro1600000v12.5qrh0qdz5000_nx3072'
'CTRLro3200000v12.5qrh0qdz5000_nx3072'
%}

%{
%%RATIO OF INITIAL RMAX TO R0
%'CTRLrodrm2.5v12.5qrh0qdz5000_nx3072'
'CTRLv12.5qrh0qdz5000_nx3072'
'CTRLrodrm10v12.5qrh0qdz5000_nx3072'
%}

%{
%%CORIOLIS
%'CTRLv0qrhSATqdz5000_nx3072_fdiv8'
'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_fx2'
'CTRLv0qrhSATqdz5000_nx3072_fx4'
%'CTRLv0qrhSATqdz5000_nx3072_fx8'
%}


%%POTENTIAL INTENSITY
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5' %40 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'    %49 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc5' %80 m/s
'CTRLv0qrhSATqdz5000_nx3072'   %93 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc1' %120 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'    %125 m/s
%'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1' %160 m/s

%}


%{
%%AMPLITUDE OF INITIAL MOISTURE ENVELOPE WITHIN WHICH VORTEX IS EMBEDDED
'CTRLv12.5qro800000qrh.1qdz5000_nx3072'
'CTRLv12.5qro800000qrh.2qdz5000_nx3072'
'CTRLv12.5qro800000qrh.3qdz5000_nx3072'
%'CTRLv12.5qro800000qrhSATqdz5000_nx3072'
%}


}; %name of sub-directory with nc files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_file = sprintf('simsets/%s.mat',sim_set);

if(exist(out_file)==2)  %data already exists for this simulation set, don't need to run it again; just go straight to plotting
    load(out_file);
    
else

numruns=length(subdirs);  %total number of runs you want to plot (taken from below)

%%want data for entire simulation
t0 = 0;
tf = 100;

%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=100000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=100;  %last z grid point [0,end]

%%Define subset of points for analysis
rmin_sub = 0;  %[km]; lowest value plotted
rmax_sub = 10000;    %[km]; highest value plotted
zmin_subsub = .5;  %[km]; lowest value plotted
zmax_subsub = 1; %[km]; highest value plotted

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

    %%DIRECTORY WITH OUTPUT DATA
    subdir_full=sprintf('%s%s',dir_in,subdir)

    %%EXTRACT TIMESTEP SIZE
    var_dt = 'qvpert'; %doesn't matter, just need any variable
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    if(run_type==3)
        numfiles=length(dir(sprintf('%s/cm1out_t*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_t0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    elseif(run_type==1)
        numfiles=length(dir(sprintf('%s/cm1out_0*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    end

    dt = time; %[s]; time of cm1out_t0001.nc is defined as zero
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
    i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging
    i_t0_usr = max(0,round(tmean0_usr*24*60*60/dt)+1); %timestep corresponding to start time for averaging
    i_tf_usr = min(numfiles,round(tmeanf_usr*24*60*60/dt)+1); %timestep corresponding to end time for averaging

    %%TIME IN DAYS
    %    t_min_day = dt*i_t0/86400;
    %    t_max_day = dt*i_tf/86400;

    %%DIRECTORY WITH OUTPUT DATA
    subdir_full=sprintf('%s%s',dir_in,subdir)

    %% EXTRACT DATA FROM input_sounding
    snd_file = 'input_sounding';

    %%load initial sounding variables
    dir_full = strcat(dir_in,subdir);
    dz_frac = dz / .625;
    clear dz
    [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc dz] = snd_extract(dir_full,snd_file,dz_frac);


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% FULL WIND %%%%%%%
    var = 'vinterp';    %'vinterp'=full wind; 'vg'=gradient wind

    %%initialize the output vectors
    clear v_df_qv v_units_qv data_tmean_usr data_tmean
    data_tmean_usr=zeros(xf-x0+1,1);
    data_tmean=zeros(xf-x0+1,tf-t0);    %1 profile per day
    Vmax_movave=nan(i_tf-i_t0+1,1);
    rmax_movave=nan(i_tf-i_t0+1,1);
    r12_movave=nan(i_tf-i_t0+1,1);
    r0_movave=nan(i_tf-i_t0+1,1);
    r0Lil_movave=nan(i_tf-i_t0+1,1);
    num_prof = 1;   %initialize iteration of time average profiles


    clear t_day
    clear r12 r0
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
        [T rho rh the s s_sat gam_m] = thermo(p,th,qv);
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
        indices = find(zvals>=zmin_subsub & zvals<=zmax_subsub);
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

        %%Calculate and save instantaneous r12 time-series
        v_usr = 12;
        data_sub_temp = smooth(data_sub,30);
        temp1=find(data_sub_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr<0 || r_usr>1000000 || max(data_sub)<20)
                r_usr=NaN;
            end
        else
            r_usr = NaN;
        end
        r12(ii)=r_usr;
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
            if(r_usr<0 || r_usr>1000000 || max(data_sub)<20)
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
            
            %%Calculate and save T-day running mean r12
            v_usr = 12;
            v_r_mean_temp = smooth(v_r_mean,30);
            temp1=find(v_r_mean_temp<=v_usr);
            temp2=find(temp1>i_max,1);
            i_vusr_out=temp1(temp2);
            clear temp1 temp2
            if(~isempty(i_vusr_out))
                v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
                v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
                r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<20)
                    r_usr_movave=NaN;
                end
            else
                r_usr_movave = NaN;
            end
            r12_movave(ii-floor(nfile_mean/2)) = r_usr_movave;
            if(~isnan(r_usr_movave))
                V_user = 12;
                r_user = r12_movave(ii-floor(nfile_mean/2));
                [junk1,junk2,r_0_full,res_error] = r0_calc(r_user,V_user,fcor,Cd_in);
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
                if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<20)
                    r_usr_movave=NaN;
                end
            else
                r_usr_movave = NaN;
            end
            r0_movave(ii-floor(nfile_mean/2)) = r_usr_movave;
            clear i_vusr_out v_out v_in i_max
        end

%}        
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
  
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% GRADIENT WIND %%%%%%%
    var = 'vg';    %'vinterp'=full wind; 'vg'=gradient wind

    %%initialize the output vectors
    clear v_df_qv v_units_qv data_tmean_usr_g data_tmean_g
    data_tmean_usr_g=zeros(xf-x0+1,1);
    data_tmean_g=zeros(xf-x0+1,tf-t0);    %1 profile per day
    Vmax_movave_g=nan(i_tf-i_t0+1,1);
    rmax_movave_g=nan(i_tf-i_t0+1,1);
    r12_movave_g=nan(i_tf-i_t0+1,1);
    r0_movave_g=nan(i_tf-i_t0+1,1);
    r0Lil_movave_g=nan(i_tf-i_t0+1,1);
    num_prof = 1;   %initialize iteration of time average profiles


    clear t_day
    clear r12 r0
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
        [T rho rh the s s_sat gam_m] = thermo(p,th,qv);
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
        indices = find(zvals>=zmin_subsub & zvals<=zmax_subsub);
        i_zvals = i_zvals(indices);
        clear indices
        
zvals_all{rr} = zvals(i_zvals);

        data_temp = squeeze(data);
        data_sub = data_temp(i_xvals',i_zvals');
        xvals_sub = xvals(i_xvals);
        data_tmean_g=data_tmean_g(i_xvals,:);
        
        %%Calculate and save instantaneous Vmax time-series
        i_temp = round(length(data_sub)/8);
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

        %%Calculate and save instantaneous r12 time-series
        v_usr = 12;
        data_sub_temp=smooth(data_sub,30);  %use 30-point smoother on profile to calculate
        temp1=find(data_sub_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr<0 || r_usr>1000000 || max(data_sub)<20)
                r_usr=NaN;
            end
        else
            r_usr = NaN;
        end
        r12(ii)=r_usr;
        clear i_vusr_out v_out v_in
        
        %%Calculate and save instantaneous r0 time-series
        v_usr = 0;
        data_sub_temp=smooth(data_sub,30);  %use 30-point smoother on profile to calculate
        temp1=find(data_sub_temp<=v_usr);
        temp2=find(temp1>i_max,1);
        i_vusr_out=temp1(temp2);
        clear temp1 temp2
        if(~isempty(i_vusr_out))
            v_out = data_sub_temp(i_vusr_out);   %v at first point where v<=v_usr
            v_in = data_sub_temp(i_vusr_out-1);  %v at last point where v>=v_usr
            r_usr = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
            if(r_usr<0 || r_usr>1000000 || max(data_sub)<20)
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
            
            %%Calculate and save T-day running mean r12
            v_usr = 12;
            v_r_mean_temp=smooth(v_r_mean,30);  %use 30-point smoother on profile to calculate
            temp1=find(v_r_mean_temp<=v_usr);
            temp2=find(temp1>i_max,1);
            i_vusr_out=temp1(temp2);
            clear temp1 temp2
            if(~isempty(i_vusr_out))
                v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
                v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
                r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<20)
                    r_usr_movave=NaN;
                end
            else
                r_usr_movave = NaN;
            end
            
            r12_movave_g(ii-floor(nfile_mean/2)) = r_usr_movave;
            if(~isnan(r_usr_movave))
                V_user = 12;
                r_user = r12_movave_g(ii-floor(nfile_mean/2));
                [junk1,junk2,r_0_full,res_error] = r0_calc(r_user,V_user,fcor,Cd_in);
                clear junk1 junk2
                
                r0Lil_movave_g(ii-floor(nfile_mean/2)) = r_0_full/1000;
            else
                r0Lil_movave_g(ii-floor(nfile_mean/2)) = NaN;
            end
            clear i_vusr_out v_out v_in
            
            %%Calculate and save T-day running mean r0
            v_usr = 0;
            v_r_mean_temp=smooth(v_r_mean,30);  %use 30-point smoother on profile to calculate
            temp1=find(v_r_mean_temp<=v_usr);
            temp2=find(temp1>i_max,1);
            i_vusr_out=temp1(temp2);
            clear temp1 temp2
            if(~isempty(i_vusr_out))
                v_out = v_r_mean_temp(i_vusr_out);   %v at first point where v<=v_usr
                v_in = v_r_mean_temp(i_vusr_out-1);  %v at last point where v>=v_usr
                r_usr_movave = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                if(r_usr_movave<0 || r_usr_movave>1000000 || max(v_r_mean)<20)
                    r_usr_movave=NaN;
                end
            else
                r_usr_movave = NaN;
            end
            r0_movave_g(ii-floor(nfile_mean/2)) = r_usr_movave;
            clear i_vusr_out v_out v_in i_max
        end

%}        
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
    tdat_end = t_day(end-dt_equil/(dt/60/60/24):end);
    Vmax_dat_end = Vmax_movave_g(end-dt_equil/(dt/60/60/24):end);
    indices = ~isnan(Vmax_dat_end);
    tdat_end = tdat_end(indices);
    Vmax_dat_end = Vmax_dat_end(indices);
    clear indices
    Vmax_mean_end = mean(Vmax_dat_end);

    %time to genesis and stats at genesis
    indices = find(Vmax_movave_g>Vmax_mean_end);
    i_gen = indices(1);
    tau_gen(rr) = t_day(i_gen); %defined using the GRADIENT wind
    Vmax_gen(rr) = Vmax_movave(i_gen);
    rmax_gen(rr) = rmax_movave(i_gen);
    r12_gen(rr) = r12_movave(i_gen);
    r0_gen(rr) = r0_movave(i_gen);
    r0Lil_gen(rr) = r0Lil_movave(i_gen);
    Vmax_gen_g(rr) = Vmax_movave_g(i_gen);
    rmax_gen_g(rr) = rmax_movave_g(i_gen);
    r12_gen_g(rr) = r12_movave_g(i_gen);
    r0_gen_g(rr) = r0_movave_g(i_gen);
    r0Lil_gen_g(rr) = r0Lil_movave_g(i_gen);
    
    %time to peak value and value
    indices = find(Vmax_movave_g(i_gen:end)==max(Vmax_movave_g(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    Vmax_tau_max_g(rr) = t_day(i_peak); %defined using the GRADIENT wind
    Vmax_max_g(rr) = Vmax_movave_g(i_peak);

    indices = find(rmax_movave_g(i_gen:end)==max(rmax_movave_g(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    rmax_tau_max_g(rr) = t_day(i_peak); %defined using the GRADIENT wind
    rmax_max_g(rr) = rmax_movave_g(i_peak);
    
    indices = find(r12_movave_g(i_gen:end)==max(r12_movave_g(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r12_tau_max_g(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r12_max_g(rr) = r12_movave_g(i_peak);
    
    indices = find(r0_movave_g(i_gen:end)==max(r0_movave_g(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r0_tau_max_g(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r0_max_g(rr) = r0_movave_g(i_peak);
    
    indices = find(r0Lil_movave_g(i_gen:end)==max(r0Lil_movave_g(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r0Lil_tau_max_g(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r0Lil_max_g(rr) = r0Lil_movave_g(i_peak);
    
    indices = find(Vmax_movave(i_gen:end)==max(Vmax_movave(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    Vmax_tau_max(rr) = t_day(i_peak); %defined using the GRADIENT wind
    Vmax_max(rr) = Vmax_movave(i_peak);

    indices = find(rmax_movave(i_gen:end)==max(rmax_movave(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    rmax_tau_max(rr) = t_day(i_peak); %defined using the GRADIENT wind
    rmax_max(rr) = rmax_movave(i_peak);
    
    indices = find(r12_movave(i_gen:end)==max(r12_movave(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r12_tau_max(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r12_max(rr) = r12_movave(i_peak);
    
    indices = find(r0_movave(i_gen:end)==max(r0_movave(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r0_tau_max(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r0_max(rr) = r0_movave(i_peak);
    
    indices = find(r0Lil_movave(i_gen:end)==max(r0Lil_movave(i_gen:end)));
    i_peak = i_gen-1+indices(1);
    r0Lil_tau_max(rr) = t_day(i_peak); %defined using the GRADIENT wind
    r0Lil_max(rr) = r0Lil_movave(i_peak);
    
    numvars = 10;     %Vmax, rmax, r12, r0, r0Lil, Vmax_g, rmax_g, r12_g, r0_g, r0Lil_g
    for l=1:numvars
        
        switch l
            case 1
                var_movave = Vmax_movave;
            case 2
                var_movave = rmax_movave;
            case 3
                var_movave = r12_movave;
            case 4
                var_movave = r0_movave;
            case 5
                var_movave = r0Lil_movave;
            case 6
                var_movave = Vmax_movave_g;
            case 7
                var_movave = rmax_movave_g;
            case 8
                var_movave = r12_movave_g;
            case 9
                var_movave = r0_movave_g;
            case 10
                var_movave = r0Lil_movave_g;
        end
        
        %%Extract data for (dt_equil)-day period at end of simulation and calculate its mean
        tdat_end = t_day(end-dt_equil/(dt/60/60/24):end);
        var_dat_end = var_movave(end-dt_equil/(dt/60/60/24):end);
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
        disequil_frac = .1;    %define disequilibrium as having mean that deviates from equilibrium value by greater than 5%
        istep = 1;  %steps backwards in time to find disequilibrium point
        t_day_temp = t_day;
        var_movave_temp = var_movave;
        equil = 0;
        dt_equil_temp = dt_equil; %30 days

        while(dt_equil_temp > T_mean)
            tdat_sub = t_day_temp(istep:ceil(dt_equil_temp/(dt/60/60/24))+istep);
            var_dat_sub = var_movave_temp(istep:ceil(dt_equil_temp/(dt/60/60/24))+istep);
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

            istep = istep + 1;    %shift up 1 timestep at a time until reach equilibrium
        end

        var_tau_equil(l) = t_day_equilf;

%        else
%            var_tau_equil(l) = NaN;

%        end
    end
    
    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil(rr) = var_equil(1);
    rmax_equil(rr) = var_equil(2);
    r12_equil(rr) = var_equil(3);
    r0_equil(rr) = var_equil(4);
    r0Lil_equil(rr) = var_equil(5);
    Vmax_equil_g(rr) = var_equil(6);
    rmax_equil_g(rr) = var_equil(7);
    r12_equil_g(rr) = var_equil(8);
    r0_equil_g(rr) = var_equil(9);
    r0Lil_equil_g(rr) = var_equil(10);

    %timescales to those values
    Vmax_tau_equil(rr) = var_tau_equil(1);
    rmax_tau_equil(rr) = var_tau_equil(2);
    r12_tau_equil(rr) = var_tau_equil(3);
    r0_tau_equil(rr) = var_tau_equil(4);
    r0Lil_tau_equil(rr) = var_tau_equil(5);
    Vmax_tau_equil_g(rr) = var_tau_equil(6);
    rmax_tau_equil_g(rr) = var_tau_equil(7);
    r12_tau_equil_g(rr) = var_tau_equil(8);
    r0_tau_equil_g(rr) = var_tau_equil(9);
    r0Lil_tau_equil_g(rr) = var_tau_equil(10);
    
    %% Save data for all simulations
    Vmax_movave_all(:,rr)=Vmax_movave;
    rmax_movave_all(:,rr)=rmax_movave;
    r12_movave_all(:,rr)=r12_movave;
    r0_movave_all(:,rr)=r0_movave;
    r0Lil_movave_all(:,rr)=r0Lil_movave;
    
    Vmax_movave_g_all(:,rr)=Vmax_movave_g;
    rmax_movave_g_all(:,rr)=rmax_movave_g;
    r12_movave_g_all(:,rr)=r12_movave_g;
    r0_movave_g_all(:,rr)=r0_movave_g;
    r0Lil_movave_g_all(:,rr)=r0Lil_movave_g;
    
    %% User profile %%%%%%%%%%%%%%%%%%%%%%%%
    xvals_sub_all{rr} = xvals_sub;
    data_tmean_usr_all{rr} = data_tmean_usr;
    data_tmean_usr_g_all{rr} = data_tmean_usr_g;


end

%% Save data to file
if(save_file == 1)
    save tempstuff.mat plot_usrprof plot_ts_multi plot_ts plot_proftevol pl_clrs rmin_plot rmax_plot datamin_plot datamax_plot plot_stats plot_gradient plot_full
    clear plot_usrprof plot_ts_multi plot_ts plot_proftevol pl_clrs rmin_plot rmax_plot datamin_plot datamax_plot plot_stats plot_gradient plot_full
    save temp.mat
    movefile('temp.mat',out_file)
    load tempstuff.mat
    delete('tempstuff.mat')
    
    cd simsets/PLOTS
    mkdir(sim_set)
    cd ../..
end

end


%% PLOTTING %%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',8,'defaultaxesfontweight','bold')
%Single simulation: Plot time-series of Vmax, rmax, r12, r0 for both V and Vg
if(plot_ts==1)  
    %Plot time series
    figure(1)
    clf(1)
    
    for i=1:1%numruns
        
        subplot(2,2,1)
        plot(t_day,Vmax_movave_all(:,i),'Color',pl_clrs{i})
        hold on
    end

    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
    legend(input_legend)


    for i=1:1%numruns
        figure(1)
        
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(Vmax_tau_equil(i)))
            h=plot(t_day(Vmax_tau_equil(i)/(dt/60/60/24)),Vmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        plot(t_day,Vmax_movave_g_all(:,i),'--','Color',pl_clrs{i})
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        hold on
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(rmax_tau_equil(i)))
            h=plot(t_day(rmax_tau_equil(i)/(dt/60/60/24)),rmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        plot(t_day,rmax_movave_g_all(:,i),'--','Color',pl_clrs{i})
        if(~isnan(rmax_tau_equil_g(i)))
            h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        hold on
%{        
        subplot(2,2,3)
        plot(t_day,r12_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r12_tau_equil(i)))
            h=plot(t_day(r12_tau_equil(i)/(dt/60/60/24)),r12_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        plot(t_day,r12_movave_g_all(:,i),'--','Color',pl_clrs{i})
        if(~isnan(r12_tau_equil_g(i)))
            h=plot(t_day(r12_tau_equil_g(i)/(dt/60/60/24)),r12_equil_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        hold on
%}        
        subplot(3,1,3)
%{
        plot(t_day,r0_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil(i)))
            h=plot(t_day(r0_tau_equil(i)/(dt/60/60/24)),r0_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
        plot(t_day,r0_movave_g_all(:,i),'--','Color',pl_clrs{i})
        if(~isnan(r0_tau_equil_g(i)))
            h=plot(t_day(r0_tau_equil_g(i)/(dt/60/60/24)),r0_equil_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',[.49 1 .63]);
        end
%}        
        %Lilly model r0
        plot(t_day,r0Lil_movave_all(:,i),pl_clrs{i})
        hold on
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil(i)/(dt/60/60/24)),r0Lil_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',.5*[.49 1 .63]);
        end
        plot(t_day,r0Lil_movave_g_all(:,i),'r--')
        if(~isnan(r0Lil_tau_equil_g(i)))
            h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',.5*[.49 1 .63]);
        end

        hold on
    end

    subplot(3,1,1)
    hold on
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    hold on
    grid on
    subplot(3,1,2)
    axis([t0 tf 0 1.1*max(max(rmax_movave_all(round(size(rmax_movave_all,1)/10):end,1)),max(rmax_movave_g_all(round(size(rmax_movave,1)/10):end,1)))])
    input_title=sprintf('rmax [km]');
    title(input_title)
    hold on
    grid on
%    subplot(2,2,3)
%    axis([t0 tf 0 1.1*max(max(r12_movave_all(round(size(r12_movave_all,1)/10):end,1)),max(r12_movave_g_all(round(size(r12_movave_all,1)/10):end,1)))])
%    input_title=sprintf('r12 [km]');
%    title(input_title)
%    hold on
    subplot(3,1,3)
    axis([t0 tf 0 1.1*max(max(r0Lil_movave_all(round(size(r0Lil_movave_all,1)/10):end,1)),max(r0Lil_movave_g_all(round(size(r0Lil_movave_all,1)/10):end,1)))])
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    input_title1=sprintf('Time evolution for FULL (solid) and GRADIENT (dashed)');
    input_title2=sprintf('%i -day moving averages',T_mean);
    suptitle({input_title1,input_title2})
    grid on
end

%Single simulation: Comparing radial wind profile with time
if(plot_proftevol==1)
    num_profs = 5;
    num_days = size(data_tmean,2);
    t_days_plot = round(size(data_tmean,2)/num_profs):round(size(data_tmean,2)/num_profs):size(data_tmean,2);

    i_xvals_pl = find(xvals_sub_all{1}>=rmin_plot & xvals_sub_all{1}<=rmax_plot);
    xvals_pl = xvals_sub_all{1}(i_xvals_pl);
    xvals_mat = repmat(xvals_pl',1,num_profs);

    
    rprofs=data_tmean(i_xvals_pl,t_days_plot);
    figure(2)
    clf(2)
    
    hold off
    %plot(data_tmean_usr(1:600))
    plot(xvals_mat,rprofs)
    axis([rmin_plot rmax_plot max(min(min(rprofs)),datamin_plot) min(max(max(rprofs)),datamax_plot)])
    hold on
    input_legend = {};
    for i=1:length(t_days_plot)
        input_legend{end+1}=num2str(t_days_plot(i));
    end
    legend(input_legend)

    %plot(data_tmean_usr_g(1:600),'r')
    %plot(smooth(data_tmean_usr_g(1:600),50),'m')
    plot(0*(1:xvals(end)),'k')
    input_title = sprintf('FULL Radial wind profiles at various times [day]; z=%5.2f %s',zvals(i_zvals),zunits);
    title(input_title)
    xlabel('Radius [km]')
    ylabel('Azimuthal wind speed [m/s]')
    grid on
    
    rprofs=data_tmean_g(i_xvals_pl,t_days_plot);
    figure(3)
    clf(3)
    
    hold off
    %plot(data_tmean_usr(1:600))
    plot(xvals_mat,rprofs)
    axis([rmin_plot rmax_plot max(min(min(rprofs)),datamin_plot) min(max(max(rprofs)),datamax_plot)])
    hold on
    input_legend = {};
    for i=1:length(t_days_plot)
        input_legend{end+1}=num2str(t_days_plot(i));
    end
    legend(input_legend)
    
    %plot(data_tmean_usr_g(1:600),'r')
    %plot(smooth(data_tmean_usr_g(1:600),50),'m')
    plot(0*(1:xvals(end)),'k')
    input_title = sprintf('GRADIENT Radial wind profiles at various times [day]; z=%5.2f %s',zvals(i_zvals),zunits);
    title(input_title)
    xlabel('Radius [km]')
    ylabel('Azimuthal wind speed [m/s]')
    grid on
    
end


%Single/Multi simulation: Plot user-defined radial wind profiles
if(plot_usrprof==1)

    if(plot_full == 1)
    figure(4)
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
    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
    legend(input_legend)
    grid on
    end

    if(plot_gradient == 1)
    figure(5)
    clf(5)
    
    hold off
    for i=1:numruns
        i_xvals_pl = find(xvals_sub_all{i}>=rmin_plot & xvals_sub_all{i}<=rmax_plot);
        xvals_pl = xvals_sub_all{i}(i_xvals_pl);

        plot(xvals_pl,data_tmean_usr_g_all{i}(i_xvals_pl),pl_clrs{i},'LineWidth',1)
        hold on
    end
    plot(0*(1:xvals_pl(end)),'k')
    input_title = sprintf('Time-mean radial GRADIENT wind profiles: days %i - %i; z=%5.2f %s',tmean0_usr,tmeanf_usr,zvals(i_zvals),zunits);
    title(input_title)
    xlabel('Radius [km]')
    ylabel('Azimuthal wind speed [m/s]')
    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
    legend(input_legend)
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'radprof','jpeg')
    cd ../../..
    end
    
end

%Multi simulation: Plot comparison of various statistical quantities
if(plot_stats==1)

    %xvals_pl = CTRL_val*multipliers;    %values defined by user at top
    xvals_pl = multipliers;    %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    
    %%Equilibrium storm characteristics
    if(plot_full == 1)
    figure(6)
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

%    data_temp = r12_equil;
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

    legend({'Vmax','rmax','r0Lil','tauvmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_equil,tf,sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_equil(i_ctrl),rmax_equil(i_ctrl),r0Lil_equil(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    figure(61)
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

%    data_temp = r12_tau_equil;
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

    legend({'Vmax','rmax','r0Lil','tauvmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_equil,tf,sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_equil(i_ctrl),rmax_tau_equil(i_ctrl),r0Lil_tau_equil(i_ctrl));
    title({input_title1,input_title3})
    grid on
    end

    if(plot_gradient == 1)
    figure(7)
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

%    data_temp = r12_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = r0Lil_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on

    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
    legend({'Vmax','rmax','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_equil,tf,sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_equil_g(i_ctrl),rmax_equil_g(i_ctrl),r0Lil_equil_g(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'equil','jpeg')
    cd ../../..
    
    %%Timescales to equilibrium
    figure(71)
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

%    data_temp = r12_tau_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-.')
%    hold on

%    data_temp = r0_tau_equil_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-.')
%    hold on
    
    data_temp = r0Lil_tau_equil_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-.')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
    legend({'tauvmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Equilibrium storm %i - %i days: log_2(Y/Y*) vs. %s [%s]',tf-dt_equil,tf,sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_equil_g(i_ctrl),rmax_tau_equil_g(i_ctrl),r0Lil_tau_equil_g(i_ctrl));
    title({input_title1,input_title3})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'tau_equil','jpeg')
    cd ../../..

    end
    
    %%Transient timescales
    if(plot_full == 1)
    figure(8)
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

%    data_temp = r12_tau_equil;
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

%    data_temp = r12_gen;
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
    
%    legend({'Vmax','rmax','r12','r0','r0Lil'},'Location','SouthEast')
    legend({'Vmax','rmax','r0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL storm at Genesis: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_gen(i_ctrl),rmax_gen(i_ctrl),r0Lil_gen(i_ctrl));
    input_title3=sprintf('Y*: taugen* = %5.2f',tau_gen(i_ctrl));
    title({input_title1,input_title2,input_title3})
    grid on
    end
    
    if(plot_gradient == 1)
    figure(11)
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

%    data_temp = r12_gen_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'kx-')
%    hold on

%    data_temp = r0_gen_g;
%    data_pl = log2(data_temp./data_temp(i_ctrl));
%    dat_max = max(dat_max,max(data_pl));
%    plot(xvals_pl,data_pl,'gx-')
%    hold on

    data_temp = r0Lil_gen_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
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
    legend({'Vmax','rmax','r0Lil','taugen'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT storm at Genesis: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_gen_g(i_ctrl),rmax_gen_g(i_ctrl),r0Lil_gen_g(i_ctrl));
    input_title3=sprintf('Y*: taugen* = %5.2f',tau_gen(i_ctrl));
    title({input_title1,input_title2,input_title3})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'gen','jpeg')
    cd ../../..

    end
    
    %%Peak values of Vmax, rmax, r0; time scales to each
    if(plot_full == 1)
    figure(12)
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
    legend({'Vmaxmax','rmaxmax','r0Lilmax','tauVmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_max(i_ctrl),rmax_max(i_ctrl),r0Lil_max(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    figure(121)
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
    legend({'Vmaxmax','rmaxmax','r0Lilmax','tauVmax','taurmax','taur0Lil'},'Location','SouthEast')
    input_title1=sprintf('FULL Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_max(i_ctrl),rmax_tau_max(i_ctrl),r0Lil_tau_max(i_ctrl));
    title({input_title1,input_title3})
    grid on
    end
    
    if(plot_gradient == 1)
    figure(13)
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

    data_temp = r0Lil_max_g;
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    plot(xvals_pl,data_pl,'cx-')
    hold on
    
    axis([-4 4 -4 4])
    ylabel('log_2(Y/Y*)')
    xlabel({sprintf('log_2(X/X*)'),sprintf('%s: X* = %5.2f [%s]',sim_set,CTRL_val,units)})
    
%    legend({'Vmax_g','rmax_g','r0_g','r0Lil_g'},'Location','SouthEast')
    legend({'Vmaxmax','rmaxmax','r0Lilmax'},'Location','SouthEast')
    input_title1=sprintf('GRADIENT Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title2=sprintf('Y*: Vmax* = %5.2f m/s; rmax* = %5.2f km; r0Lil* = %5.2f km',Vmax_max_g(i_ctrl),rmax_max_g(i_ctrl),r0Lil_max_g(i_ctrl));
    title({input_title1,input_title2})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'max','jpeg')
    cd ../../..

    figure(131)
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
    
    data_temp = r0Lil_tau_max_g;
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
    input_title1=sprintf('GRADIENT Maximum values: log_2(Y/Y*) vs. %s [%s]',sim_set,units);
    input_title3=sprintf('Y*: tauV* = %5.2f d; taurmax* = %5.2f d; taur0* = %5.2f d',Vmax_tau_max_g(i_ctrl),rmax_tau_max_g(i_ctrl),r0Lil_tau_max_g(i_ctrl));
    title({input_title1,input_title3})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
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
    clf(14)
    
    
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_all(:,i),'Color',pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(Vmax_movave_all))])

    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),Vmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(Vmax_tau_max(i)))
            h=plot(t_day(Vmax_tau_max(i)/(dt/60/60/24)),Vmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(Vmax_tau_equil(i)))
            h=plot(t_day(Vmax_tau_equil(i)/(dt/60/60/24)),Vmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_all(:,i),'Color',pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(rmax_movave_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),rmax_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(rmax_tau_max(i)))
            h=plot(t_day(rmax_tau_max(i)/(dt/60/60/24)),rmax_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(rmax_tau_equil(i)))
            h=plot(t_day(rmax_tau_equil(i)/(dt/60/60/24)),rmax_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
%{
    %r12
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day,r12_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r12_tau_equil(i)))
            h=plot(t_day(r12_tau_equil(i)/(dt/60/60/24)),r12_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    input_title=sprintf('r12 [km]');
    title(input_title)
    axis([0 100 0 1.1*max(max(r12_movave_all(100:end,:)))])
%}
%{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day,r0_movave_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil(i)))
            h=plot(t_day(r0_tau_equil(i)/(dt/60/60/24)),r0_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    input_title=sprintf('r0 [km]');
    title(input_title)
    axis([0 100 0 1.1*max(max(r0_movave_all(100:end,:)))])
%}
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day,r0Lil_movave_all(:,i),'Color',pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(r0Lil_movave_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),r0Lil_gen(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(r0Lil_tau_max(i)))
            h=plot(t_day(r0Lil_tau_max(i)/(dt/60/60/24)),r0Lil_max(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil(i)/(dt/60/60/24)),r0Lil_equil(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
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
%    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day(1:i_gen(i))/tau_gen(i),Vmax_movave_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
        plot(t_day(1:i_gen(i))/tau_gen(i),rmax_movave_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
        plot(t_day(1:i_gen(i))/tau_gen(i),r0Lil_movave_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
    %Plot time series
    figure(15)
    clf(15)
    
    %legend info
    runs_pl = CTRL_val*(2.^multipliers);        %values defined by user at top
    i_ctrl = find(multipliers==0,1);
    input_legend = {};
    for i=1:length(runs_pl)
        input_legend{end+1}=num2str(runs_pl(i));
    end
%    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day,Vmax_movave_g_all(:,i),'Color',pl_clrs{i})
        hold on
        
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('Vmax [m/s]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(Vmax_movave_g_all))])

    for i=1:numruns
        
        subplot(3,1,1)
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),Vmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(Vmax_tau_max_g(i)))
            h=plot(t_day(Vmax_tau_max_g(i)/(dt/60/60/24)),Vmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(Vmax_tau_equil_g(i)))
            h=plot(t_day(Vmax_tau_equil_g(i)/(dt/60/60/24)),Vmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    %rmax
    for i=1:numruns
        
        subplot(3,1,2)
        plot(t_day,rmax_movave_g_all(:,i),'Color',pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('rmax [km]');
    title(input_title)
    ylabel(input_title);
%    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(rmax_movave_g_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),rmax_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(rmax_tau_max(i)))
            h=plot(t_day(rmax_tau_max_g(i)/(dt/60/60/24)),rmax_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(rmax_tau_equil(i)))
            h=plot(t_day(rmax_tau_equil_g(i)/(dt/60/60/24)),rmax_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
%{
    %r12
    for i=1:numruns
        
        subplot(2,2,3)
        plot(t_day,r12_movave_g_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r12_tau_equil_g(i)))
            h=plot(t_day(r12_tau_equil_g(i)/(dt/60/60/24)),r12_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    input_title=sprintf('r12 [km]');
    title(input_title)
    axis([0 100 0 1.1*max(max(r12_movave_g_all(100:end,:)))])
%}
%{
    %r0
    for i=1:numruns
        
        subplot(2,2,4)
        plot(t_day,r0_movave_g_all(:,i),'Color',pl_clrs{i})
        hold on
        if(~isnan(r0_tau_equil_g(i)))
            h=plot(t_day(r0_tau_equil_g(i)/(dt/60/60/24)),r0_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
    
    input_title=sprintf('r0 [km]');
    title(input_title)
    axis([0 100 0 1.1*max(max(r0_movave_g_all(100:end,:)))])
%}
    %r0Lil
    for i=1:numruns
        
        subplot(3,1,3)
        plot(t_day,r0Lil_movave_g_all(:,i),'Color',pl_clrs{i})
        hold on
    end
    grid on
    legend(input_legend,'Location','EastOutside')
    input_title=sprintf('r0Lil [km]');
    title(input_title)
    ylabel(input_title);
    xlabel('time [days]')
    axis([0 100 0 1.1*max(max(r0Lil_movave_g_all(20:end,:)))])

    for i=1:numruns
        
        hold on
        if(~isnan(tau_gen(i)))
            h=plot(t_day(tau_gen(i)/(dt/60/60/24)),r0Lil_gen_g(i),'d');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(r0Lil_tau_max_g(i)))
            h=plot(t_day(r0Lil_tau_max_g(i)/(dt/60/60/24)),r0Lil_max_g(i),'^');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
        if(~isnan(r0Lil_tau_equil(i)))
            h=plot(t_day(r0Lil_tau_equil_g(i)/(dt/60/60/24)),r0Lil_equil_g(i),'s');
            set(h,'markersize',10,'MarkerFaceColor',pl_clrs{i});
        end
    end
        
    hold on
    input_title1=sprintf('GRADIENT storm evolutions for variable %s [%s]',sim_set,units);
    input_title2=sprintf('%i -day moving averages',T_mean);
    input_title3='';
    suptitle({input_title1,input_title2,input_title3})
    grid on
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'ts','jpeg')
    cd ../../..
    
%% Zoom in on genesis period; plot from 0:tau_gen
    figure(17)
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
%    input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%    legend(input_legend)

    
    %Vmax
    dat_max = 0;
    for i=1:numruns
        
        subplot(3,1,1)
        plot(t_day(1:i_gen(i))/tau_gen(i),Vmax_movave_g_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
        plot(t_day(1:i_gen(i))/tau_gen(i),rmax_movave_g_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
        plot(t_day(1:i_gen(i))/tau_gen(i),r0Lil_movave_g_all(1:i_gen(i),i),'Color',pl_clrs{i})
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
    
    cd(sprintf('simsets/PLOTS/%s',sim_set))
    saveas(gcf,'pregen','jpeg')
    cd ../../..

    end
end

