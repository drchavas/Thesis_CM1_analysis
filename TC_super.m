%TC_super.m

%Purpose: Super script that can run TC_stats.m and all subsequent
%post-processing.

clear
clc
close all

tic

%% Which simula tions (or sets) should I run?
%name out output subdir (within simsets_Tmean#/PLOTS/[sim_set]/) where plots will be saved
%sim_sets_all = {'Lx' 'dx' 'dz' 'lh' 'lv' 'qro' 'ro' 'fcor' 'Tsst' 'Ttpp' 'usfc' 'usfc_drag' 'Qcool' 'nondim' 'nondim1.5' 'nondim2' 'mpi' 'Cd' 'QcoolVpcnst' 'QcoolVplvHcnst'}; 
%sim_sets_all = {'Tsst_drag' 'Ttpp_drag' 'usfc_drag' 'Qcool_drag' 'nondim2_drag' 'Cd_drag'}; 
%sim_sets_all = {'test'}; 
sim_sets_all = {'transient'}; 
    %IF 'single'
    sim_single = 'CTRLv0qrhSATqdz5000_nx3072_Cddiv2_drag';    %runs only this simulation
moist = 1;  %1 = moist; else = dry
    
%% Which scripts should I run?
run_TC_stats = 0;
    save_file = 0;  %for TC_stats only; note: program will not overwrite old file
run_TC_stats_dynamicequil = 0;  %overwrites old file automatically
run_TC_structure_ts = 0;    %overwrites old plot automatically
run_TC_stats_plot = 1;  %overwrites old [simset].mat file and plots automatically
    run_TC_stats_plot_dynamicequil = 1; %also run for dynamic equilibrium

%% Parameters for scripts
v_usr_fracVp = .1;  %wind speed as fraction of Vp; beyond this radius, radiative subsidence radial wind profile should apply.
T_mean = 2; %[days]; averaging time period used to calculate moving time-average radial profile from which rmax and r0 are calculated
dt_equil = 10;  %[days]; how long must be quasi-steady to define equilibrium
    %%For static equilibrium (equil_dynamic = 0):
    dt_final = 30;  %[day]; length of static equilibrium period
    tf = 100;   %[day]; end of static equilibrium period
    %%For dynamic (most stable) equilibrium (equil_dynamic = 1):
    dt_final_dynamic = 10;  %[days]; length of most stable period after day 60 over which equilibrium is calculated

%Plotting domain (TC_stats_plot.m only)
rmin_plot = 0;  %[km]
rmax_plot = 1500;    %[km]

%IMPORTANT: for calculating the outer radius using control values of the constant parameters
wrad_ctrl = .0026;   %control run RCE value
Cd_in_ctrl = .0015;
t_start_dynequil = 60;   %[days]; start here and move forward to find dynamic equilibrium

%% Basic crap
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
%subdir_pre='';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1='CHAVAS_CM1_FINAL'; 2='CHAVAS_CM1_FINAL_nodrag'
run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

%%Data time range (want data for entire simulation)
t0 = 0;

%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=100000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=1000;  %last z grid point [0,end]

%%Define subset of points for analysis
rmin_sub = 0;  %[km]; lowest value
rmax_sub = 100000;    %[km]; highest value
zmin_subsub = .5;  %[km]; lowest value
zmax_subsub = 1; %[km]; highest value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%user-defined times for radial profile
tmean0_usr = tf-dt_final;    %[day]
tmeanf_usr = tf;    %[day]

dir_home = pwd;

for jj = 1:length(sim_sets_all)

    sim_set = sim_sets_all{jj};
    run_type = run_types(jj);

    switch sim_set
        case 'test'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';
            subdirs_set = {

'CTRLv0qrhSATqdz5000_nx3072_Tthresh248K_rad4.0K_drag'

            }
            multipliers = ones(length(subdirs_set),1);
        case 'DRY'
            CTRL_val = 1
            units = '-';
            subdirs_set = {
                'CTRLv0qrhSATqdz5000_nx3072_DRYdrc'
            }
            multipliers = ones(length(subdirs_set),1);
            
        case 'lh_DRY'
            CTRL_val = 1500; %CTRL value of quantity varied across simulations
            units = 'm';
            multipliers = [-2 -1 0 1 2];
            subdirs_set = {

            'CTRLv0qrhSATqdz5000_nx3072_lh375_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_lh750_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_lh3000_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_lh6000_DRY'
            }

        case 'fcor_DRY_kerry'
            CTRL_val = 5; %CTRL value of quantity varied across simulations
            units = '**10^{-5} s^{-1}';
            %multipliers = [-2 -1 0 1 2];
            multipliers = [-2 -1 0 1 2];
            subdirs_set = {

            'CTRLv0qrhSATqdz5000_nx3072_fdiv4_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_fdiv2_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_fx2_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_fx4_DRY'
            }
        
        case 'Ttpp_DRY'
            CTRL_val = 200; %CTRL value of quantity varied across simulations
            units = 'K';
            multipliers = log2([250 225 200 175 150]/CTRL_val);
            subdirs_set = {
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh225K_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_DRY'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_DRY'
            }
        case 'Lx'
            CTRL_val = 12288; %CTRL value of quantity varied across simulations
            units = 'km';
            multipliers = [-4 -3 -2 -1 0];
            subdirs_set = {
                %%DOMAIN SIZE
                'CTRLv0qrhSATqdz5000_nx192'
                'CTRLv0qrhSATqdz5000'
                'CTRLv0qrhSATqdz5000_nx768'
                'CTRLv0qrhSATqdz5000_nx1536'
                'CTRLv0qrhSATqdz5000_nx3072'
            }
        case 'dx'
            CTRL_val = 4; %CTRL value of quantity varied across simulations
            units = 'km';
            multipliers = [-1 0 1 2 3];
            subdirs_set = {
                %%HORIZONTAL RESOLUTION
                %'CTRLv0qrhSATqdz5000_nx12288_dx1' TOO SLOW
                'CTRLv0qrhSATqdz5000_nx6144_dx2'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx1536_dx8'
                'CTRLv0qrhSATqdz5000_nx768_dx16'
                'CTRLv0qrhSATqdz5000_nx384_dx32'
            }
        case 'dz'
            CTRL_val = 625; %CTRL value of quantity varied across simulations
            units = 'm';
            multipliers = [-2 -1 0 1];
            subdirs_set = {
                %%VERTICAL RESOLUTION
                'CTRLv0qrhSATqdz5000_nx3072_dz156.25'
                'CTRLv0qrhSATqdz5000_nx3072_dz312.5'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_dz1250'

            }
        case 'lh'
            CTRL_val = 1500; %CTRL value of quantity varied across simulations
            units = 'm';
            multipliers = [-3 -2 -1 0 1 2 3];
            subdirs_set = {
                %%HORIZONTAL MIXING LENGTH
                'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
                'CTRLv0qrhSATqdz5000_nx3072_lh375'
                'CTRLv0qrhSATqdz5000_nx3072_lh750'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_lh3000'
                'CTRLv0qrhSATqdz5000_nx3072_lh6000'
                'CTRLv0qrhSATqdz5000_nx3072_lh12000'
            }
        case 'lv'
            CTRL_val = 100; %CTRL value of quantity varied across simulations
            units = 'm';
            multipliers = [-3 -2 -1 0 1 2 3];
            subdirs_set = {
                %%VERTICAL MIXING LENGTH
                'CTRLv0qrhSATqdz5000_nx3072_lv12.5'
                'CTRLv0qrhSATqdz5000_nx3072_lv25'
                'CTRLv0qrhSATqdz5000_nx3072_lv50'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_lv200'
                'CTRLv0qrhSATqdz5000_nx3072_lv400'
                'CTRLv0qrhSATqdz5000_nx3072_lv800'
            }
        case 'qro'
            CTRL_val = 200; %CTRL value of quantity varied across simulations
            units = 'km';
            multipliers = [-3 -2 -1 0 1 2 3];
            subdirs_set = {
                %%HORIZ SCALE OF MOISTURE PERTURBATION
                'CTRLv0qro25000qrhSATqdz5000_nx3072'
                'CTRLv0qro50000qrhSATqdz5000_nx3072'
                'CTRLv0qro100000qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qro400000qrhSATqdz5000_nx3072'
                'CTRLv0qro800000qrhSATqdz5000_nx3072'
                'CTRLv0qro1600000qrhSATqdz5000_nx3072'
            }
        case 'ro'
            CTRL_val = 400; %CTRL value of quantity varied across simulations
            units = 'km';
            multipliers = [-3 -2 -1 0 1 2 3];
            subdirs_set = {
                %%HORIZ SCALE OF WIND PERTURBATION
                'CTRLro50000v12.5qrh0qdz5000_nx3072'
                'CTRLro100000v12.5qrh0qdz5000_nx3072'
                'CTRLro200000v12.5qrh0qdz5000_nx3072'
                'CTRLv12.5qrh0qdz5000_nx3072'
                'CTRLro800000v12.5qrh0qdz5000_nx3072'
                'CTRLro1600000v12.5qrh0qdz5000_nx3072'
                'CTRLro3200000v12.5qrh0qdz5000_nx3072'
            }
        case 'fcor'
            CTRL_val = 5; %CTRL value of quantity varied across simulations
            units = '**10^{-5} s^{-1}';
            multipliers = [-2 -1 0 1 2 3];
            subdirs_set = {
                %%CORIOLIS
                %'CTRLv0qrhSATqdz5000_nx3072_fdiv8' %-- HITS DOMAIN WALL
                'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
                'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_fx2'
                'CTRLv0qrhSATqdz5000_nx3072_fx4'
                'CTRLv0qrhSATqdz5000_nx3072_fx8'
            }
        case 'Tsst_drag'
            CTRL_val = 300; %CTRL value of quantity varied across simulations
            units = 'K';
            %multipliers = log2([285 287.5 290 292.5 295 297.5 300 302.5 305 310]/CTRL_val);
            multipliers = log2([285 290 295 300 305 310]/CTRL_val);
            subdirs_set = {
                %%Tsst
                %'CTRLv0qrhSATqdz5000_nx3072_SST275.00K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_SST285.00K_drag'
                %'CTRLv0qrhSATqdz5000_nx3072_SST287.50K'
                'CTRLv0qrhSATqdz5000_nx3072_SST290.00K_drag'
                %'CTRLv0qrhSATqdz5000_nx3072_SST292.50K'
                'CTRLv0qrhSATqdz5000_nx3072_SST295.00K_drag'
                %'CTRLv0qrhSATqdz5000_nx3072_SST297.50K'
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                %'CTRLv0qrhSATqdz5000_nx3072_SST302.50K'
                'CTRLv0qrhSATqdz5000_nx3072_SST305.00K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_SST310.00K_drag'
            }
        case 'Ttpp_drag'
            CTRL_val = 200; %CTRL value of quantity varied across simulations
            units = 'K';
            multipliers = log2([250 225 200 175 150]/CTRL_val);
            subdirs_set = {
                %%Ttpp
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh225K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_drag'
            }
        case 'usfc_drag'
            CTRL_val = 3; %CTRL value of quantity varied across simulations
            units = 'ms^{-1}';
            %multipliers = log2([.5 1 2 3 4 5 10]/CTRL_val);
            multipliers = log2([.5 1 2 3 5 10]/CTRL_val);
            subdirs_set = {
                %%usfc
                'CTRLv0qrhSATqdz5000_nx3072_usfc.5_drag'
                'CTRLv0qrhSATqdz5000_nx3072_usfc1_drag'
                'CTRLv0qrhSATqdz5000_nx3072_usfc2_drag'
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                %'CTRLv0qrhSATqdz5000_nx3072_usfc4_drag'
                'CTRLv0qrhSATqdz5000_nx3072_usfc5_drag'
                'CTRLv0qrhSATqdz5000_nx3072_usfc10_drag'
            }
        case 'Qcool_drag'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = 'K day^{-1}';
            multipliers = [-2 -1 0 1 2];
            subdirs_set = {
                %%Qcool
                'CTRLv0qrhSATqdz5000_nx3072_rad0.125K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_rad0.25K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                'CTRLv0qrhSATqdz5000_nx3072_rad1.0K_drag'
                'CTRLv0qrhSATqdz5000_nx3072_rad2.0K_drag'
            }
        case 'Cd_drag'
            CTRL_val = 1.5; %CTRL value of quantity varied across simulations
            units = '**10^-3';
            %multipliers = [-3 -2 -1 0 1 2 3];
            multipliers = [-3 -2 -1 0 1];
            subdirs_set = {
                %%DRAG COEFFICIENT
                'CTRLv0qrhSATqdz5000_nx3072_Cddiv8_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Cddiv4_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Cddiv2_drag'
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Cdx2_drag'
                %%'CTRLv0qrhSATqdz5000_nx3072_Cdx4_drag' % -- wildly out of equilibrium
                %%'CTRLv0qrhSATqdz5000_nx3072_Cdx8_drag' % -- wildly out of equilibrium
            }
        case 'QcoolVpcnst'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = 'K day^{-1}';
            multipliers = [-1 0 1 2 3 4 5];
            subdirs_set = {
                %%Qcool at constant Vp (can't get Vp ~ 92 m/s with rad0.125)
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K'
            }
        case 'QcoolVplvHcnst'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = 'K day^{-1}';
            %multipliers = [-1 0 1 2 3 4 5];
            multipliers = [-1 0 1 2];
            subdirs_set = {
                %%Q_rad at constant Vp, lv/H constant (can't get Vp ~ 92 m/s with rad0.125)
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K_rad0.25K_lv121'
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh220K_rad1.0K_lv81'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh237K_rad2.0K_lv60'
                %'CTRLv0qrhSATqdz5000_nx3072_Tthresh247K_rad4.0K_lv44'
                %'CTRLv0qrhSATqdz5000_nx3072_Tthresh255K_rad8.0K_lv26'
                %'CTRLv0qrhSATqdz5000_nx3072_Tthresh257K_rad16.0K_lv19'
            }
        case 'nondim2_drag'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';

            subdirs_set = {
                'CTRLv0qrhSATqdz5000_nx3072_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2_lh750_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_lh750_drag'
                'CTRLv0qrhSATqdz5000_nx3072_fdiv2_lh750_drag'
                'CTRLv0qrhSATqdz5000_nx3072_fx2_lh3000_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh3000_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_drag'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_lh3000_drag'
            }
            multipliers = ones(length(subdirs_set),1);
        case 'nondim'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';

            subdirs_set = {
                'CTRLv0qrhSATqdz5000_nx3072'

%}            
                %%Ttpp
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh225K'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'
           
                %%CORIOLIS
                %%'CTRLv0qrhSATqdz5000_nx3072_fdiv8' -- HITS DOMAIN WALL
                'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
                'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_fx2'
                'CTRLv0qrhSATqdz5000_nx3072_fx4'
                'CTRLv0qrhSATqdz5000_nx3072_fx8'
                
                %%HORIZONTAL MIXING LENGTH
                'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
                'CTRLv0qrhSATqdz5000_nx3072_lh375'
                'CTRLv0qrhSATqdz5000_nx3072_lh750'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_lh3000'
                'CTRLv0qrhSATqdz5000_nx3072_lh6000'
                'CTRLv0qrhSATqdz5000_nx3072_lh12000'
            }
            multipliers = ones(length(subdirs_set),1);
        case 'nondim1.5'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';

            subdirs_set = {
            'CTRLv0qrhSATqdz5000_nx3072'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_lh3000'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2_lh750'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh3000'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_lh750'
            'CTRLv0qrhSATqdz5000_nx3072_fdiv2_lh750'
            'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2'
            'CTRLv0qrhSATqdz5000_nx3072_fx2_lh3000'
            }
            multipliers = ones(length(subdirs_set),1);
            
        case 'nondim_all'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';

            subdirs_set = {
                'CTRLv0qrhSATqdz5000_nx3072'
            
                %%1-parm
                %%Ttpp
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh225K'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh175K'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K'
           
                %%CORIOLIS
                %%'CTRLv0qrhSATqdz5000_nx3072_fdiv8' -- HITS DOMAIN WALL
                'CTRLv0qrhSATqdz5000_nx3072_fdiv4'
                'CTRLv0qrhSATqdz5000_nx3072_fdiv2'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_fx2'
                'CTRLv0qrhSATqdz5000_nx3072_fx4'
                'CTRLv0qrhSATqdz5000_nx3072_fx8'
                
                %%HORIZONTAL MIXING LENGTH
                'CTRLv0qrhSATqdz5000_nx3072_lh187.5'
                'CTRLv0qrhSATqdz5000_nx3072_lh375'
                'CTRLv0qrhSATqdz5000_nx3072_lh750'
                %'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_lh3000'
                'CTRLv0qrhSATqdz5000_nx3072_lh6000'
                'CTRLv0qrhSATqdz5000_nx3072_lh12000'

                %2-parm
                'CTRLv0qrhSATqdz5000_nx3072'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_lh3000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv2_lh750'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh3000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_lh750'
                'CTRLv0qrhSATqdz5000_nx3072_fdiv2_lh750'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2'
                'CTRLv0qrhSATqdz5000_nx3072_fx2_lh3000'
                
                %3-parm (extreme values of nondim only)
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv2_lh187.5'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv4_lh750'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fdiv4_lh750'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh187.5'

                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fdiv4'

                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_lh12000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh3000'

                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fx2_lh12000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh150K_fx4_lh6000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx2_lh12000'
                'CTRLv0qrhSATqdz5000_nx3072_Tthresh250K_fx4_lh6000'
            }
            multipliers = ones(length(subdirs_set),1);
            
        case 'transient'
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';
            subdirs_set = {
                'CTRLv0qrhSATqro400000qdz5000_nx3072_fdiv2_lh3000_drag'
                'CTRLv0qrhSATqro100000qdz5000_nx3072_fx2_lh750_drag'
            
%                'CTRLv0qrhSATqdz5000_nx3072'
%                'CTRLv0qro100000qrhSATqdz5000_nx3072_Tthresh250K_lh750_75day'
%                'CTRLv0qro100000qrhSATqdz5000_nx3072_fx2_lh750_75day'
%                'CTRLv0qro400000qrhSATqdz5000_nx3072_fdiv2_lh3000_75day'
            }
            multipliers = ones(length(subdirs_set),1);
%{
        case 'qroVpf'
            CTRL_val = .25; %CTRL value of quantity varied across simulations
            units = 'km';
            multipliers = [-3 -2 -1 0 1 2 3];
            subdirs_set = {
                %%HORIZ SCALE OF MOISTURE PERTURBATION
                'CTRLv0qrhSATqroVpfdiv128qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv64qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv32qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv16qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv8qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv4qdz5000_nx3072'
                'CTRLv0qrhSATqroVpfdiv2qdz5000_nx3072'                
                'CTRLv0qrhSATqroVpfqdz5000_nx3072'                
            }
%}
        case 'single'
            subdirs_set = {
                sim_single
            }
            CTRL_val = 1; %CTRL value of quantity varied across simulations
            units = '-';
            multipliers = [1];
        otherwise
            sprintf('Simulation set does not exist!')
    end
    
    dir_in1 = sprintf('../CM1_postproc_data/sim');
    if(mod(T_mean,1)>0)
        dir_in3 = sprintf('Tmean%3.2f_%i_%i',T_mean,tf-dt_final,tf);
        dir_in3dyn = sprintf('Tmean%3.2f_dt%i_dynamic',T_mean,dt_final_dynamic);
    else
        dir_in3 = sprintf('Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
        dir_in3dyn = sprintf('Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
    dir_in_dat = sprintf('%sdata_%s',dir_in1,dir_in3);
    dir_in_set = sprintf('%ssets_%s',dir_in1,dir_in3);
    dir_in_plt = sprintf('%splots_%s',dir_in1,dir_in3);
    dir_in_dat_dyn = sprintf('%sdata_%s',dir_in1,dir_in3dyn);
    dir_in_set_dyn = sprintf('%ssets_%s',dir_in1,dir_in3dyn);
    dir_in_plt_dyn = sprintf('%splots_%s',dir_in1,dir_in3dyn);
    assert(isdir(dir_in_dat),'Output directory doesnt exist!')
    
    for ii = 1:length(subdirs_set)
        
        subdir = subdirs_set{ii};
        
        if(run_TC_stats == 1)
            [junk] = TC_stats(subdir_pre,ext_hd,run_type,t0,tf,tmean0_usr,tmeanf_usr,v_usr_fracVp,T_mean,dt_equil,dt_final,save_file,subdir,x0,xf,y0,yf,z0,zf,rmin_sub,rmax_sub,zmin_subsub,zmax_subsub,dir_home,moist,dir_in_dat,wrad_ctrl,Cd_in_ctrl);
        end
        
        if(run_TC_stats_dynamicequil == 1)
            [junk] = TC_stats_dynamicequil(tf,T_mean,dt_final,dt_final_dynamic,subdir,dir_home,dir_in_dat,dir_in_dat_dyn,t_start_dynequil);
        end
        
        if(run_TC_structure_ts == 1)
            %%this subroutine runs for both static and dynamic equilibrium automatically
            [junk] = TC_structure_ts(run_type,T_mean,dt_final,t0,tf,dt_final_dynamic,subdir,dir_home,dir_in_dat,dir_in_dat_dyn,dir_in_plt,dir_in_plt_dyn);
        end
        

    end
    
    if(run_TC_stats_plot == 1)
        %%fixed equilibrium period
        %if(strmatch(sim_set,'nondimall')~=1)
            equil_dynamic = 0;
            [junk] = TC_stats_plot(run_type,T_mean,equil_dynamic,dt_final,tf,dt_final_dynamic,rmin_plot,rmax_plot,CTRL_val,units,multipliers,subdirs_set,sim_set,dir_home,dir_in_dat,dir_in_dat_dyn,dir_in_set,dir_in_set_dyn);
        %end
        
        %%most stable equilibrium period
        if(run_TC_stats_plot_dynamicequil == 1)
            equil_dynamic = 1;
            [junk] = TC_stats_plot(run_type,T_mean,equil_dynamic,dt_final,tf,dt_final_dynamic,rmin_plot,rmax_plot,CTRL_val,units,multipliers,subdirs_set,sim_set,dir_home,dir_in_dat,dir_in_dat_dyn,dir_in_set,dir_in_set_dyn);
        end
    end
    
%%Run separately (does not save plots)
%MPI_collapse_V
%MPI_collapse_r
%Dimensional_scaling
%Nondimensional_scaling

end

toc
