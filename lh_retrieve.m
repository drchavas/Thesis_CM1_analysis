%lh_retrieve.m

%Created: 17 Feb 2012, Dan Chavas

%Purpose: To look up the value of l_h from a given simulation. Called by
%TC_stats.m

function [lh] = lh_retrieve(subdir_full);

%clc
%clear

%sample input filename: CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375
%sample mpi filename: input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_rad1K_mpi

%subdir_full='/Volumes/CHAVAS_CM1/CM1_output/axisym/CTRL_icRCE/CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375/';
lh_file = sprintf('%s/namelist.input',subdir_full);

%%Open namelist.input file

FID = fopen(lh_file, 'r');
Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
CStr = Data{1};
clear Data
fclose(FID);

temp = strfind(CStr,' l_h ');
nonemptyCells = ~cellfun(@isempty,temp);
i_lh = find(nonemptyCells==1);
temp = CStr{i_lh};
lh = sscanf(temp, '%*s %*s %f', [1, inf]);

end
