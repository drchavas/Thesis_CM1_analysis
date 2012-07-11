%TC_stats_newequil.m

%Created: 01 June 12, Dan Chavas

%Purpose: This file copies the file output by TC_stats, opens it,
%recalculates the equilibrium values for each structural variable given a
%new set of days, and then resaves the file in a new sub-directory

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

T_mean = 5; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
dt_final = 50;  %[days]; original length of period over which equilibrium is calculated
tf = 150;   %[days]; original end day of equilibrium calculation
dt_final_neww = 30;  %[days]; new length of period over which equilibrium is calculated
tf_neww = 60;   %[days]; new end day of equilibrium calculation

save_file_newequil = 1;

%%Simulations
%Vpots = 93*ones(1000,1);

subdirs_neweq = {

%%DOMAIN SIZE
'CTRLv0qrhSATqdz5000_nx192'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_nx768'
'CTRLv0qrhSATqdz5000_nx1536'
'CTRLv0qrhSATqdz5000_nx3072'
%}

%%HORIZONTAL RESOLUTION
'CTRLv0qrhSATqdz5000_nx6144_dx2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx1536_dx8'
'CTRLv0qrhSATqdz5000_nx768_dx16'
'CTRLv0qrhSATqdz5000_nx384_dx32'
%}

%%VERTICAL RESOLUTION
'CTRLv0qrhSATqdz5000_nx3072_dz156.25'
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

%%VERTICAL MIXING LENGTH
'CTRLv0qrhSATqdz5000_nx3072_lv12.5'
'CTRLv0qrhSATqdz5000_nx3072_lv25'
'CTRLv0qrhSATqdz5000_nx3072_lv50'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_lv200'
'CTRLv0qrhSATqdz5000_nx3072_lv400'
'CTRLv0qrhSATqdz5000_nx3072_lv800'
%}
%{
%%HORIZ SCALE OF MOISTURE PERTURBATION AS FRAC OF Vp/f
'CTRLv0qrhSATqroVpfdiv128qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv64qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv32qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv16qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv8qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv4qdz5000_nx3072'
'CTRLv0qrhSATqroVpfdiv2qdz5000_nx3072'
'CTRLv0qrhSATqroVpfqdz5000_nx3072'
%}

%%HORIZ SCALE OF MOISTURE PERTURBATION
'CTRLv0qro25000qrhSATqdz5000_nx3072'
'CTRLv0qro50000qrhSATqdz5000_nx3072'
'CTRLv0qro100000qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qro400000qrhSATqdz5000_nx3072'
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
%{
%%RATIO OF INITIAL RMAX TO R0
%'CTRLrodrm2.5v12.5qrh0qdz5000_nx3072'
'CTRLv12.5qrh0qdz5000_nx3072'
'CTRLrodrm10v12.5qrh0qdz5000_nx3072'
%}
%{
%%COMBINDATIONS OF VORTEX AND MOISTURE AS FRAC OF Vp/f
'CTRLv12.5qrhSATqroVpfdiv128qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv64qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv32qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv16qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv8qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv4qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfdiv2qdz5000_nx3072'
'CTRLv12.5qrhSATqroVpfqdz5000_nx3072'
%}
%{
'CTRLro800000v12.5qro100000qrhSATqdz5000_nx3072'
'CTRLv12.5qrhSATqdz5000_nx3072'
'CTRLv12.5qro50000qrhSATqdz5000_nx3072'
'CTRLv12.5qro800000qrh.1qdz5000_nx3072'
'CTRLv12.5qro800000qrh.2qdz5000_nx3072'
'CTRLv12.5qro800000qrh.3qdz5000_nx3072'
'CTRLv12.5qro800000qrhSATqdz5000_nx3072'
%}

%%NO INITIAL PERTURBATION
%%'CTRLv0qrh0qdz5000_nx3072'
%}

%%CORIOLIS
'CTRLv0qrhSATqdz5000_nx3072_fdiv8'
'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_fx2'
'CTRLv0qrhSATqdz5000_nx3072_fx4'
%%'CTRLv0qrhSATqdz5000_nx3072_fx8' %line1275 error Tmean5
%}

%%POTENTIAL INTENSITY (via Qrad, usfc, SST, H (i.e. Tthresh))
%%T_sst
'CTRLv0qrhSATqdz5000_nx3072_SST275.00K'
'CTRLv0qrhSATqdz5000_nx3072_SST285.00K'
'CTRLv0qrhSATqdz5000_nx3072_SST290.00K'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K'
'CTRLv0qrhSATqdz5000_nx3072_SST297.50K'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_SST302.50K'
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K'
'CTRLv0qrhSATqdz5000_nx3072_SST310.00K'

%%T_sst_r0radconst
'CTRLv0qrhSATqdz5000_nx3072_SST275.00K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_SST285.00K_r0radconst'
%'CTRLv0qrhSATqdz5000_nx3072_SST290.00K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_SST297.50K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_SST302.50K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_SST310.00K_r0radconst'

%%T_tpp
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh225K'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'

%%u_sfc
'CTRLv0qrhSATqdz5000_nx3072_usfc.5'
'CTRLv0qrhSATqdz5000_nx3072_usfc1'
'CTRLv0qrhSATqdz5000_nx3072_usfc2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_usfc4'
'CTRLv0qrhSATqdz5000_nx3072_usfc5'
'CTRLv0qrhSATqdz5000_nx3072_usfc10'

%%Q_rad
'CTRLv0qrhSATqdz5000_nx3072_rad0.125K'
'CTRLv0qrhSATqdz5000_nx3072_rad0.25K'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_rad1.0K'
'CTRLv0qrhSATqdz5000_nx3072_rad2.0K'

%%Q_rad_r0radconst
'CTRLv0qrhSATqdz5000_nx3072_rad0.125K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_rad0.25K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_rad1.0K_r0radconst'
'CTRLv0qrhSATqdz5000_nx3072_rad2.0K_r0radconst'

%%DRAG COEFFICIENT
'CTRLv0qrhSATqdz5000_nx3072_Cddiv8'
'CTRLv0qrhSATqdz5000_nx3072_Cddiv4'
'CTRLv0qrhSATqdz5000_nx3072_Cddiv2'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Cdx2'
'CTRLv0qrhSATqdz5000_nx3072_Cdx4'
'CTRLv0qrhSATqdz5000_nx3072_Cdx8'

%%Q_rad at constant Vp
'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K'

%%Q_rad at constant Vp, lv/H constant
'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K_lv121'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K_lv81'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K_lv60'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K_lv44'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K_lv26'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K_lv19'

};
skipfornow = {

%%Other (i.e. combinations; f, l_h fixed at CTRL value)
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5' %40 m/s
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_usfc5'   %74 m/s
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_Tthresh150K_nz56'
'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_Tthresh250K_usfc1'   %71
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_nz56'
%}

%%Non-dim scaling: f and/or l_h varied from CTRL value
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh750'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv4_lh750'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv2_lh375'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fdiv4_lh375'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2_lh6000'
%'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx4_lh3000'    %too small vmax
%'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx4_lh6000'    %too small vmax
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_fx2_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh250K_usfc5_lh750'

'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh750'

'CTRLv0qrhSATqdz5000_nx3072_usfc5'
'CTRLv0qrhSATqdz5000_nx3072_usfc5_fdiv2'
'CTRLv0qrhSATqdz5000_nx3072_usfc5_fx2'
'CTRLv0qrhSATqdz5000_nx3072_usfc5_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_usfc5_lh750'

'CTRLv0qrhSATqdz5000_nx3072_fdiv2_lh375'
'CTRLv0qrhSATqdz5000_nx3072_fdiv2_lh750'
'CTRLv0qrhSATqdz5000_nx3072_fdiv4_lh375'
'CTRLv0qrhSATqdz5000_nx3072_fdiv4_lh750'
'CTRLv0qrhSATqdz5000_nx3072_fx2_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_fx2_lh6000'
'CTRLv0qrhSATqdz5000_nx3072_fx4_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_fx4_lh6000'

'CTRLv0qrhSATqdz5000_nx3072_usfc1'
'CTRLv0qrhSATqdz5000_nx3072_usfc1_fdiv2_nz56'
'CTRLv0qrhSATqdz5000_nx3072_usfc1_fx2_nz56'

'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fx2_lh3000_nz56'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fx4_lh3000_nz56'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fx2_lh6000_nz56'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fx4_lh6000_nz56'
%'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fdiv2_lh750_nz56'
'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fdiv4_lh750_nz56'
%'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_Tthresh150K_usfc3_fdiv2_lh375_nz56'

'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_fdiv2_nz56'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_fdiv2_lh750'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_fx2_nz56'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_lh3000_nz56'
'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_usfc1_lh750'
%}

%% HALF-DOMAIN %%%%%
%'CTRLv0qrhSATqdz5000_nx1536_dx4'
%'CTRLv0qrhSATqdz5000_nx1536_dx4_dz312.5'
%'CTRLv0qrhSATqdz5000_nx1536_dx4'
%'CTRLv0qrhSATqdz5000_nx3072_dx2'


}; %name of sub-directory with nc files

%{
subdirs_neweq = {
'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
'CTRLv0qrhSATqdz5000_nx3072_lh375'
'CTRLv0qrhSATqdz5000_nx3072_lh750'
'CTRLv0qrhSATqdz5000_nx3072'
'CTRLv0qrhSATqdz5000_nx3072_lh3000'
'CTRLv0qrhSATqdz5000_nx3072_lh6000'
'CTRLv0qrhSATqdz5000_nx3072_lh12000'
};
%}
subdirs_neweq = {'CTRLv0qrhSATqdz5000_nx3072'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stop_overwrite = 0;

if(save_file_newequil==1)
%check if new subdirectory exists
if(exist(sprintf('simdata_Tmean%i_%i_%i',T_mean,tf_neww-dt_final_neww,tf_neww))==7)

    %check to make sure you're not overwriting original TC_stats output
    if(exist(sprintf('simdata_Tmean%i_%i_%i/NOTANORIGINALTCSTATS.mat',T_mean,tf_neww-dt_final_neww,tf_neww))==2)   %TRUE if directory exists but does NOT contain this file
    else
        sprintf('You are trying to overwrite original TC_stats output!  Bad.')
        stop_overwrite = 1;
    end

    
else
    
    mkdir(sprintf('simdata_Tmean%i_%i_%i',T_mean,tf_neww-dt_final_neww,tf_neww))
    mkdir(sprintf('simplots_Tmean%i_%i_%i',T_mean,tf_neww-dt_final_neww,tf_neww))
    mkdir(sprintf('simsets_Tmean%i_%i_%i',T_mean,tf_neww-dt_final_neww,tf_neww))
    mkdir(sprintf('simsets_Tmean%i_%i_%i/PLOTS',T_mean,tf_neww-dt_final_neww,tf_neww))
    
    cd(sprintf('simdata_Tmean%i_%i_%i',T_mean,tf_neww-dt_final_neww,tf_neww))
    save NOTANORIGINALTCSTATS.mat T_mean tf dt_final tf_neww dt_final_neww
    cd ..
    
end
end

if(stop_overwrite == 0)

for rr=1:length(subdirs_neweq)

    clearvars -except rr subdirs_neweq subdir_pre ext_hd run_types T_mean tf tf_neww dt_equil dt_final dt_final_neww save_file_newequil

    subdir = subdirs_neweq{rr};

    if(exist(sprintf('simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdir))==2)
        sprintf('Updating equilibrium values for ax%s',subdir)

        %%Load data for given simulation to fix things within if needed
        load(sprintf('simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf-dt_final,tf,subdirs_neweq{rr}));
tf_neww
dt_final_neww
        numvars = 10;     %Vmax, rmax, r12, r0, r0Lil, Vmax_g, rmax_g, r12_g, r0_g, r0Lil_g
        for l=6:numvars     %SKIP THE FULL WIND ONES FOR NOW

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
                    sprintf('Vmax_g')
                case 7
                    var_movave = rmax_movave_g;
                    sprintf('rmax_g')
                case 8
                    var_movave = r12_movave_g;
                    sprintf('r12_g')
                case 9
                    var_movave = r0_movave_g;
                    sprintf('r0_g')
                case 10
                    var_movave = r0Lil_movave_g;
                    sprintf('r0Lil_g')
            end

            %%Extract data for (dt_final)-day period at end of simulation and calculate its mean
            tdat_end = t_day((tf_neww-dt_final_neww)/(dt/60/60/24):tf_neww/(dt/60/60/24));
            var_dat_end = var_movave((tf_neww-dt_final_neww)/(dt/60/60/24):tf_neww/(dt/60/60/24));
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
        end

    else
        

    end

    %% Equilibrium data %%%%%%%%%%%%%%%%%%%%%%%%
    %variable values
    Vmax_equil_g_sim = var_equil(6);
    rmax_equil_g_sim = var_equil(7);
    r12_equil_g_sim = var_equil(8);
    r0_equil_g_sim = var_equil(9);
    r0Lil_equil_g_sim = var_equil(10);

    %% Save data to file
    if(save_file_newequil == 1)
        dt_final_temp = dt_final;
        tf_temp = tf;
        dt_final = dt_final_neww;
        tf = tf_neww;
        save tempstuff.mat rr l
        clear rr l i
        save temp.mat
        load tempstuff.mat
        movefile('temp.mat',sprintf('simdata_Tmean%i_%i_%i/ax%s.mat',T_mean,tf_neww-dt_final_neww,tf_neww,subdir))
        delete('tempstuff.mat')
        dt_final = dt_final_temp;
        tf = tf_temp;
        clear dt_final_temp tf_temp
    end

end

end