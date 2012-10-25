%nc_extract.m

%Created: 02 Dec 2010, Dan Chavas

%Purpose: to extract subsets of data from a netcdf file

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!


function [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf)

%% USER INPUT %%%%%%%%%%%%%%%%%%
%direc = 'test_96_thrandpert'; %name of directory with nc files
%nc_file = 'cm1out_t0080.nc';    %nc file of interest
%var = 'winterp';  %variable of interest
%x0=0;   %first x grid point [0,end]
%xf=100;   %first y grid point [0,end]
%y0=0;   %first y grid point [0,end]
%yf=100;   %last z grid point [0,end]
%z0=3;  %first z grid point [0,end]
%zf=3;  %last z grid point [0,end]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHANGE TO PROPER DIRECTORY
dir_start=pwd;
dir_tmp = strcat(dir_in,subdir);
cd(dir_tmp)
dir_tmp

%% OPEN FILE AND EXTRACT FILE ID#
%ncid = netcdf.open(nc_file,'NOWRITE');
ncid = netcdf.open(nc_file,'WRITE');
%[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)

%% GET # GRIDPTS IN EACH DIMENSION (CM1:nx(0),ny(1),nz(2),time(6))
%NX
dimid = netcdf.inqDimID(ncid,'nx');
[junk, nx] = netcdf.inqDim(ncid,dimid);

%NY
dimid = netcdf.inqDimID(ncid,'ny');
[junk, ny] = netcdf.inqDim(ncid,dimid);

%NZ
dimid = netcdf.inqDimID(ncid,'nz');
[junk, nz] = netcdf.inqDim(ncid,dimid);

%MODIFY xf, yf, zf IF SMALLER THAN DIMENSION BOUNDS
if(x0 < 0)
    x0=0;
    sprintf('X0 input cannot be negative; reset to 0')
end
if(y0 < 0)
    y0=0;
    sprintf('Y0 input cannot be negative; reset to 0')
end
if(z0 < 0)
    z0=0;
    sprintf('Z0 input cannot be negative; reset to 0')
end
if(xf < 0)
    xf=0;
    sprintf('XF input cannot be negative; reset to 0')
end
if(yf < 0)
    yf=0;
    sprintf('YF input cannot be negative; reset to 0')
end
if(zf < 0)
    zf=0;
    sprintf('ZF input cannot be negative; reset to 0')
end

%MODIFY xf, yf, zf IF LARGER THAN DIMENSION BOUNDS
if(xf > nx-1)
    xf=nx-1;
    sprintf('XF input too large; reset to max = %i',xf)
end
if(yf > ny-1)
    yf=ny-1;
    sprintf('YF input too large; reset to max = %i',yf)
end
if(zf > nz-1)
    zf=nz-1;
    sprintf('ZF input too large; reset to max = %i',zf)
end
if(x0 > nx-1)
    x0=nx-1;
    sprintf('XF input too large; reset to max = %i',x0)
end
if(y0 > ny-1)
    y0=ny-1;
    sprintf('YF input too large; reset to max = %i',y0)
end
if(z0 > nz-1)
    z0=nz-1;
    sprintf('ZF input too large; reset to max = %i',z0)
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
if(strmatch(netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),14),'z_delta'))
    dz = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'z_delta');
else
    dz = (zmax-zmin)/nz;
end

%% EXTRACT TIME SINCE START OF SIMULATION
%% EXTRACT TIME
varid_t = netcdf.inqVarID(ncid,'time');
%[varname_t,xtype_t,dimids_t,natts_t]=netcdf.inqVar(ncid,varid_t);
time = netcdf.getVar(ncid,varid_t)
t_units='s';
%t_day = time / 86400;

%% EXTRACT CORIOLIS PARAMETER
varid = netcdf.inqVarID(ncid,'f_cor');
fcor = netcdf.getVar(ncid,varid);

%%%%%%%%%%%%%%% END GLOBAL STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALCULATE nx_sub, ny_sub, nz_sub (# pts in subset domain)
nx_sub=xf-x0+1;
ny_sub=yf-y0+1;
nz_sub=zf-z0+1;

%% EXTRACT VARIABLE ID#
varid = netcdf.inqVarID(ncid,var);
[varname,xtype,dimids,natts]=netcdf.inqVar(ncid,varid);
varname

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
v_def_flag=0;
for m=0:natts-1     %check if variable as "def" attribute
    if ( strmatch(netcdf.inqAttName(ncid,varid,m),'def')==1 )
        v_def = netcdf.getAtt(ncid,varid,'def');
        v_def_flag = 1;
    end
end
if(v_def_flag == 0)
    v_def = varname;
end

%UNITS
v_units = netcdf.getAtt(ncid,varid,'units');

%XMIN OF VARIABLE
v_xmin = netcdf.getAtt(ncid,varid,'x_min');

%YMIN OF VARIABLE
v_ymin = netcdf.getAtt(ncid,varid,'y_min');
v_ymin = v_ymin - dz/2 + dy/2;   %DRC 25 Oct 2012; COMPENSATE FOR DUMB ERROR IN 3D STITCH FILE THAT ADDS DZ/2 INSTEAD OF DY/2

%ZMIN OF VARIABLE
if(num_dim>=3)
    v_zmin = netcdf.getAtt(ncid,varid,'z_min');
else
    v_zmin=NaN;
end

%% EXTRACT SUBSET DATA
%SYNTAX: netcdf.getVar(ncid,varid, [start], [count])
%NOTE: if error, be sure dimensions are correct!
if(num_dim == 4)
    %4D data (e.g. qv)
    data = netcdf.getVar(ncid,varid,[x0 y0 z0 0],[nx_sub ny_sub nz_sub 1]);
        %start = [x0 y0 z0 t0]
        %count = [nx ny nz nt] (including first point, i.e. 1=single point)
        %NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!
elseif(num_dim == 3)
    %3D data (e.g. qv0)
    data = netcdf.getVar(ncid,varid,[x0 y0 z0],[nx_sub ny_sub nz_sub]);
        %start = [x0 y0 z0]
        %count = [nx ny nz] (including first point, i.e. 1=single point)
        %NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!
else
    %2D data (e.g. rain)
    data = netcdf.getVar(ncid,varid,[x0 y0],[nx_sub ny_sub]);
end

%% OUTPUT TO SCREEN FOR VERIFICATION
%sprintf('Size of your dataset (X*Y*Z) is %i * %i * %i',nx_sub,ny_sub,nz_sub)
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
cd(dir_start)
end