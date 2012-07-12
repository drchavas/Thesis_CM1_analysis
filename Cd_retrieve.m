%Cd_retrieve.m

%Created: 12 Jul 2012, Dan Chavas

%Purpose: To look up the value of C_d from a given simulation. Called by
%TC_stats.m

%%TESTING
%{
subdir_full = '/Volumes/CHAVAS_CM1_FINAL/CM1_output/axisym/CTRL_icRCE/CTRLv0qrhSATqdz5000_nx3072'
%}

function [Cd] = Cd_retrieve(subdir_full);

%clc
%clear

%sample input filename: CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375
%sample mpi filename: input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_rad1K_mpi

%subdir_full='/Volumes/CHAVAS_CM1/CM1_output/axisym/CTRL_icRCE/CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375/';
Cd_file = sprintf('%s/namelist.input',subdir_full)

%%Open namelist.input file

FID = fopen(Cd_file, 'r');
Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
CStr = Data{1};
clear Data
fclose(FID);

temp = strfind(CStr,' cnstcd ');
nonemptyCells = ~cellfun(@isempty,temp);
i_Cd = find(nonemptyCells==1);
temp = CStr{i_Cd};
Cd = sscanf(temp, '%*s %*s %f', [1, inf]);

end
