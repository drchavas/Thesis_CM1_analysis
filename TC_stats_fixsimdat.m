%TC_stats.m

%Created: 10-05-11, Dan Chavas
%Updated: 12-02-11, Dan Chavas

%Purpose: This file performs all the calculations to characterize the
%structure and evolution of the storms of interest.  The output is then
%saved in individual files (one per simulation) in directory simdata/

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc



%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

%%Simulations
subdirs = {

%%DOMAIN SIZE
'CTRLv0qrhSATqdz5000_nx192'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_nx768'
'CTRLv0qrhSATqdz5000_nx1536'
'CTRLv0qrhSATqdz5000_nx3072'
%}

%%HORIZONTAL RESOLUTION
%'CTRLv0qrhSATqdz5000_nx12288_dx1'
'CTRLv0qrhSATqdz5000_nx6144_dx2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx1536_dx8'
'CTRLv0qrhSATqdz5000_nx768_dx16'
'CTRLv0qrhSATqdz5000_nx384_dx32'
%}

%%VERTICAL RESOLUTION
%'CTRLv0qrhSATqdz5000_nx3072_dz156.25'
'CTRLv0qrhSATqdz5000_nx3072_dz312.5'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_dz1250'
%}

%%HORIZONTAL MIXING LENGTH
'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
'CTRLv0qrhSATqdz5000_nx3072_lh375'
'CTRLv0qrhSATqdz5000_nx3072_lh750'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_lh6000'
'CTRLv0qrhSATqdz5000_nx3072_lh12000'
%}

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

%%HORIZ SCALE OF WIND PERTURBATION
'CTRLro50000v12.5qrh0qdz5000_nx3072'
'CTRLro100000v12.5qrh0qdz5000_nx3072'
'CTRLro200000v12.5qrh0qdz5000_nx3072'
'CTRLv12.5qrh0qdz5000_nx3072'
'CTRLro800000v12.5qrh0qdz5000_nx3072'
'CTRLro1600000v12.5qrh0qdz5000_nx3072'
'CTRLro3200000v12.5qrh0qdz5000_nx3072'
%}

%%RATIO OF INITIAL RMAX TO R0
%'CTRLrodrm2.5v12.5qrh0qdz5000_nx3072'
'CTRLv12.5qrh0qdz5000_nx3072'
'CTRLrodrm10v12.5qrh0qdz5000_nx3072'
%}

%%COMBINDATIONS OF VORTEX AND MOISTURE
'CTRLro800000v12.5qro100000qrhSATqdz5000_nx3072'
'CTRLv12.5qrhSATqdz5000_nx3072'
'CTRLv12.5qro50000qrhSATqdz5000_nx3072'
'CTRLv12.5qro800000qrh.1qdz5000_nx3072'
'CTRLv12.5qro800000qrh.2qdz5000_nx3072'
'CTRLv12.5qro800000qrh.3qdz5000_nx3072'
'CTRLv12.5qro800000qrhSATqdz5000_nx3072'
%}

%%NO INITIAL PERTURBATION
'CTRLv0qrh0qdz5000_nx3072'
%}

%%CORIOLIS
%'CTRLv0qrhSATqdz5000_nx3072_fdiv8'
'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_fx2'
'CTRLv0qrhSATqdz5000_nx3072_fx4'
%'CTRLv0qrhSATqdz5000_nx3072_fx8'
%}

%%POTENTIAL INTENSITY (via Qrad, usfc, SST, H (i.e. Tthresh))
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5' %40 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'    %49 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc5' %80 m/s
'CTRLv0qrhSATqdz5000_nx3072'   %93 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc1' %120 m/s
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'    %125 m/s
%'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1' %160 m/s
%}

%%NON-DIM NUMBER SCALING Vpot/(f*lh)
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'
'CTRLv0qrhSATqdz5000_nx3072_fx2'
'CTRLv0qrhSATqdz5000_nx3072_lh3000'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_lh750'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh750'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc1_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc1_fx2'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc5_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh200K_usfc5_fx2'
%}

%% HALF-DOMAIN %%%%%
'CTRLv0qrhSATqdz5000_nx1536_dx4'
'CTRLv0qrhSATqdz5000_nx1536_dx4_dz312.5'
'CTRLv0qrhSATqdz5000_nx1536_dx4'
'CTRLv0qrhSATqdz5000_nx3072_dx2'


}; %name of sub-directory with nc files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rr=1:length(subdirs)

%%Append ax or 3d to subdir names for plot
for i=1:length(subdirs)
    if(run_types(i)==1) %ax
        subdirs_save{i}=sprintf('ax%s',subdirs{i});
    else    %3d
        subdirs_save{i}=sprintf('3d%s',subdirs{i});
    end
end

if(exist(sprintf('simdata/%s.mat',subdirs_save{rr}))==2)
    sprintf('Data already saved for %s',subdirs_save{rr})

    file_in=sprintf('simdata/%s.mat',subdirs_save{rr});
    save tempstuff.mat
    clearvars -except file_in
    
    %%Load data for given simulation
    load(file_in);
    clear Vpots Vpot rr l i ii
    save(file_in);
    clear
    
    %%Restore old variables
    load tempstuff.mat
    delete('tempstuff.mat')

end

end