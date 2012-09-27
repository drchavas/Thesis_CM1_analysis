%nc_extract_axisym.m

%Adapted to axisym output: 25 Jan 2011, Dan Chavas

%Purpose: to extract subsets of data from a netcdf file and plot them

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!


function [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf)

%% USER INPUT %%%%%%%%%%%%%%%%%%
%direc = 'CTRL1'; %name of directory with nc files
%nc_file = 'cm1out_0200.nc';    %nc file of interest
%var = 'winterp';  %variable of interest
%x0=0;   %first x grid point [0,end]
%xf=100;   %first y grid point [0,end]
%y0=0;   %first y grid point [0,end]
%yf=0;   %last z grid point [0,end]
%z0=0;  %first z grid point [0,end]
%zf=20;  %last z grid point [0,end]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHANGE TO PROPER DIRECTORY
dir_start=pwd;
dir_tmp = strcat(dir_in,subdir);
cd(dir_tmp);
dir_tmp;

%% OPEN FILE AND EXTRACT FILE ID#
%ncid = netcdf.open(nc_file,'NOWRITE');
ncid = netcdf.open(nc_file,'WRITE');

%%TESTING TO FIGURE OUT WHATS IN THE FILE
%[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)
%dimIDs = netcdf.inqDimIDs(ncid)


%% GET # GRIDPTS IN EACH DIMENSION (CM1:nx(0),ny(1),nz(2),time(6))
%NX
dimid = netcdf.inqDimID(ncid,'ni');
[junk, nx] = netcdf.inqDim(ncid,dimid);

%NY
dimid = netcdf.inqDimID(ncid,'nj');
[junk, ny] = netcdf.inqDim(ncid,dimid);

%NZ
dimid = netcdf.inqDimID(ncid,'nk');
[junk, nz] = netcdf.inqDim(ncid,dimid);

%MODIFY xf, yf, zf IF SMALLER THAN DIMENSION BOUNDS
if(x0 < 0)
    x0=0;
    sprintf('X0 input cannot be negative; reset to 0');
end
if(y0 < 0)
    y0=0;
    sprintf('Y0 input cannot be negative; reset to 0');
end
if(z0 < 0)
    z0=0;
    sprintf('Z0 input cannot be negative; reset to 0');
end
if(xf < 0)
    xf=0;
    sprintf('XF input cannot be negative; reset to 0');
end
if(yf < 0)
    yf=0;
    sprintf('YF input cannot be negative; reset to 0');
end
if(zf < 0)
    zf=0;
    sprintf('ZF input cannot be negative; reset to 0');
end

%MODIFY xf, yf, zf IF LARGER THAN DIMENSION BOUNDS
if(xf > nx-1)
    xf=nx-1;
    sprintf('XF input too large; reset to max = %i',xf);
end
if(yf > ny-1)
    yf=ny-1;
    sprintf('YF input too large; reset to max = %i',yf);
end
if(zf > nz-1)
    zf=nz-1;
    sprintf('ZF input too large; reset to max = %i',zf);
end
if(x0 > nx-1)
    x0=nx-1;
    sprintf('XF input too large; reset to max = %i',x0);
end
if(y0 > ny-1)
    y0=ny-1;
    sprintf('YF input too large; reset to max = %i',y0);
end
if(z0 > nz-1)
    z0=nz-1;
    sprintf('ZF input too large; reset to max = %i',z0);
end

%% GET GLOBAL ATTRIBUTES FOR DATASET
%XMIN
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),0);
xmin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'x_min');

%XMAX
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),1);
xmax = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'x_max');

%YMIN
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),6);
ymin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'y_min');

%YMAX
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),7);
ymax = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'y_max');

%XUNITS
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),3);
xunits = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'x_units');

%YUNITS
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),11);
yunits = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'y_units');

%ZMIN
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),12);
zmin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'z_min');

%ZMAX
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),13);
zmax = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'z_max');

%ZUNITS
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),14);
zunits = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'z_units');

%DX
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),2);
dx = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'x_delta');

%DY
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),8);
dy = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'y_delta');

%DZ
%name_tmp = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),14);
dz = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'z_delta');

%%  IF SUBSET, CALCULATE nx_sub, ny_sub, nz_sub (# pts in subset domain)
nx_sub=xf-x0+1;
ny_sub=yf-y0+1;
nz_sub=zf-z0+1;

%% EXTRACT TIME SINCE START OF SIMULATION
t_id = netcdf.inqVarID(ncid,'time');

%VARIABLE DEFINITION
name_tmp = netcdf.inqAttName(ncid,t_id,0);
t_def = netcdf.getAtt(ncid,t_id,name_tmp);

%UNITS
name_tmp = netcdf.inqAttName(ncid,t_id,1);
t_units = netcdf.getAtt(ncid,t_id,name_tmp);

%GET THE TIME
time = netcdf.getVar(ncid,t_id,0,1);
%t_day = time / 86400;

%% EXTRACT CORIOLIS PARAMETER
varid = netcdf.inqVarID(ncid,'f_cor');
fcor = netcdf.getVar(ncid,varid);

%%%%%%%%%%%%%%% END GLOBAL STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXTRACT VARIABLE ID#
varid = netcdf.inqVarID(ncid,var);
[varname,xtype,dimids,natts]=netcdf.inqVar(ncid,varid);
%length(dimids)

%% EXTRACT VARIABLE INFO
%NUMBER OF DIMENSIONS OF DATA (E.G. RAIN IS xy ONLY)
if(length(dimids) == 2)
    num_dim = 2;
elseif(length(dimids) == 3)
    num_dim = 3;    %3d time-INdependent data only (e.g. initial base state)
else
    num_dim = 4;    %note: 3d time-dependent data has length(dimids)=4 even
                        %if datafile is only one snapshop in time!
end

%VARIABLE DEFINITION
v_def = netcdf.getAtt(ncid,varid,'def');

%UNITS
v_units = netcdf.getAtt(ncid,varid,'units');

%XMIN OF VARIABLE
v_xmin = netcdf.getAtt(ncid,varid,'x_min');

%YMIN OF VARIABLE
v_ymin = netcdf.getAtt(ncid,varid,'y_min');

%ZMIN OF VARIABLE
if(num_dim>3)
    v_zmin = netcdf.getAtt(ncid,varid,'z_min');
else
    v_zmin=NaN;
end

%%TESTING
%[varname,xtype,dimids,natts]=netcdf.inqVar(ncid,varid)
%% EXTRACT DATA
%EXTRACT SUBSET OF DATA FOR DESIRED VARIABLE: netcdf.getVar(ncid,varid, [start], [count])
%NOTE: if error, be sure dimensions are correct!
if(num_dim == 4)
    %4D data
    data = netcdf.getVar(ncid,varid,[x0 y0 z0 0],[nx_sub ny_sub nz_sub 1]);
        %start = [x0 y0 z0 t0]
        %count = [nx ny nz nt] (including first point, i.e. 1=single point)
        %NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!
elseif(num_dim == 3)
    %3D data
    if(sum(dimids==[0 1 6])==length(dimids))    %XY; e.g. rain
        data = netcdf.getVar(ncid,varid,[x0 y0 0],[nx_sub ny_sub 1]);
        %start = [x0 y0 t0]
        %count = [nx ny t1] (including first point, i.e. 1=single point)
        %NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!
%    elseif((sum(dimids==[0 1 2])==length(dimids)) || (sum(dimids==[3 1 2])==length(dimids)))    %3D constant in t; e.g. q0
    else
        data = netcdf.getVar(ncid,varid,[x0 y0 z0],[nx_sub ny_sub nz_sub]);
        %start = [x0 y0 z0]
        %count = [nx ny nz] (including first point, i.e. 1=single point)
        %NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!
    end
else
    %2D data
    data = netcdf.getVar(ncid,varid,[x0 y0],[nx_sub ny_sub]);
end

%% OUTPUT TO SCREEN FOR VERIFICATION
%sprintf('Dimensions of your dataset are:'),size(data)
%lowest_level_data_xy=data(:,:,1)
%data_min = min(min(data))
%data_max = max(max(data))

%% CLOSE NETCDF FILE
netcdf.close(ncid);

%% SOME SIMPLE POST-PROCESSING
%AXIS LABELS OF SUBSET REGION
xmin_sub=v_xmin+dx*x0;
xmax_sub=v_xmin+dx*xf;
ymin_sub=v_ymin+dy*y0;
ymax_sub=v_ymin+dy*yf;
zmin_sub=v_zmin+dz*z0;
zmax_sub=v_zmin+dz*zf;

%% RETURN TO ORIGINAL DIRECTORY
cd(dir_start);

end