%mpi_retrieve.m

%Created: 07 Feb 2012, Dan Chavas

%Purpose: To look up the appropriate RCE MPI file and retrieve the mpi value (Vmax, [ms-1])

function [mpi] = mpi_retrieve(file_in);

cd input_soundings/

%clc
%clear

%sample input filename: CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375
%sample mpi filename: input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_rad1K_mpi

%file_in='CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375_rad1K'

usedrag=strfind(file_in,'drag');
if(isempty(usedrag))
    drag_pre='soundings_nodrag/'
    drag_str='';
else
    drag_pre='';
    drag_str='_drag';
end

isDRY=strfind(file_in,'DRY');
if(isempty(isDRY)) 
    dry_str='';
else    %dry simulation
    dry_str='_DRY';
end

i_sst=strfind(file_in,'SST');
if(isempty(i_sst))
    sst_str='300.00';
else
    i_sstK=strfind(file_in,'K');
    i_sstK=i_sstK(find(i_sstK>i_sst,1));
    sst_str=file_in(i_sst+3:i_sstK-1);
end

i_tpp=strfind(file_in,'Tthresh');
if(isempty(i_tpp))
    tpp_str='200';
else
    i_tppK=strfind(file_in,'K');
    i_tppK=i_tppK(find(i_tppK>i_tpp,1));
    tpp_str=file_in(i_tpp+7:i_tppK-1);
end

i_usfc=strfind(file_in,'usfc');
if(isempty(i_usfc))
    usfc_str='3';
else
    i_usfcend=strfind(file_in,'_');
    i_usfcend=i_usfcend(find(i_usfcend>i_usfc,1));
    if(isempty(i_usfcend))  %occurs if "usfc?" is the end of the filename
        usfc_str=file_in(i_usfc+4:end);
    else
        usfc_str=file_in(i_usfc+4:i_usfcend-1);
    end
end

i_rad=strfind(file_in,'rad');
if(isempty(i_rad))
    rad_str='0.5';
else
    i_radK=strfind(file_in,'K');
    i_radK=i_radK(find(i_radK>i_rad,1));
    rad_str=file_in(i_rad+3:i_radK-1);
end

if(isempty(i_rad))
    mpi_file = sprintf('%sinput_sounding_3dRCE_nx48_SST%sK_Tthresh%sK_usfc%s%s%s_mpi',drag_pre,sst_str,tpp_str,usfc_str,drag_str,dry_str);
else
    mpi_file = sprintf('%sinput_sounding_3dRCE_nx48_SST%sK_Tthresh%sK_usfc%s_rad%sK%s%s_mpi',drag_pre,sst_str,tpp_str,usfc_str,rad_str,drag_str,dry_str);
end

cd ..

[mpi] = mpi_getfromfile(mpi_file);


end
