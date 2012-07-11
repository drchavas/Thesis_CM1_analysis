%nc_extract_stats.m

%Create: 31 Jan 2011, Dan Chavas

%Purpose: to extract data from cm1out_stats.m

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!


%clear
%clc

function [data v_def v_units time t_units] = nc_extract_stats(dir_in,subdir,var)

%% USER INPUT %%%%%%%%%%%%%%%%%%
%dir_in='/Users/drchavas/Documents/Research/Thesis/CM1/v15/axisym/CM1_output/'
%subdir = 'CTRL_ptype6_ax'; %name of directory with nc files
%var = 'vmax';  %variable of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NAME OF STATS FILE
nc_file = 'cm1out_stats.nc';    %nc file of interest

%% CHANGE TO PROPER DIRECTORY
dir_start=pwd;
dir_tmp = strcat(dir_in,subdir);
cd(dir_tmp)
dir_tmp

%% OPEN FILE AND EXTRACT FILE ID#
%ncid = netcdf.open(nc_file,'NOWRITE');
ncid = netcdf.open(nc_file,'WRITE');

%%TESTING TO FIGURE OUT WHATS IN THE FILE
%[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid)
%dimIDs = netcdf.inqDimIDs(ncid)

%% EXTRACT TIME
t_id = netcdf.inqVarID(ncid,'time');

%VARIABLE DEFINITION
name_tmp = netcdf.inqAttName(ncid,t_id,0);
t_def = netcdf.getAtt(ncid,t_id,name_tmp);

%UNITS
name_tmp = netcdf.inqAttName(ncid,t_id,1);
t_units = netcdf.getAtt(ncid,t_id,name_tmp);

%GET TIME VECTOR (for entire simulation)
time = netcdf.getVar(ncid,t_id);
%t_day = time / 86400;

%%%%%%%%%%%%%%% END GLOBAL STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXTRACT VARIABLE ID#
varid = netcdf.inqVarID(ncid,var);
[varname,xtype,dimids,natts]=netcdf.inqVar(ncid,varid);
varname
%length(dimids)

%% EXTRACT VARIABLE INFO

%VARIABLE DEFINITION
name_tmp = netcdf.inqAttName(ncid,varid,0);
v_def = netcdf.getAtt(ncid,varid,name_tmp);

%UNITS
if(natts==2)
    name_tmp = netcdf.inqAttName(ncid,varid,1);
    v_units = netcdf.getAtt(ncid,varid,name_tmp);
else
    v_units = '';   %some data do not have units assigned
end
%GET DATA
data = netcdf.getVar(ncid,varid);

%% CLOSE NETCDF FILE
netcdf.close(ncid);

%% RETURN TO ORIGINAL DIRECTORY
cd(dir_start)

%}

end