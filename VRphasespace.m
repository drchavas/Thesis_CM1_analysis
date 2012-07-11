%VRphasespace.m

%Created: 09-20-11, Dan Chavas

%Purpose: Plot storm evolution in the phase space of Vmax and R12

clear
clc
figure(1)
clf(1)     %clear figures without closing them
fclose('all')

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_types=ones(1000,1); %[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3D

t0 = 0;
tf = 100;

%mult_by = [8 4 2 1 1/2 1/4 1/8];

plot_type = 1;  %0=no plot

%%Define sub-domain of data
x0=0;   %first x grid point [0,end]
xf=1000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=1;  %first z grid point [0,end]
zf=1;  %last z grid point [0,end]
v_usr = 12; %find radius of this wind speed; 'MAX' for V_max

%IF plot_type=1: input variable and domain
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    vmin_plot = 20;  %lowest value plotted
    vmax_plot = 100; %highest value plotted
    rmin_plot = 0;  %lowest value plotted
    rmax_plot = 700;    %highest value plotted
    plot_normalized = 0;    %1=plot V/Vmax, r/rmax; 0=plot regular V(r)
    numpt_sm = 1;  %smoother uses this number of points

subdirs = {

'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_lh3000'
%'CTRLv0qrhSATqdz5000_Tthresh250K'

%{
%%lh
'CTRLv0qrhSATqdz5000_lh0'
'CTRLv0qrhSATqdz5000_lh50'
'CTRLv0qrhSATqdz5000_lh150'
'CTRLv0qrhSATqdz5000_lh375'
'CTRLv0qrhSATqdz5000_lh750'
%'CTRLv0qrhSATqdz5000'
%'CTRLv0qrhSATqdz5000_lh3000'
%'CTRLv0qrhSATqdz5000_lh6000'
%}

%{
%%SELF-SIMILARITY
%'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_x2'
'CTRLv0qrhSATqdz5000_x4'
'CTRLv0qrhSATqdz5000_x8'
'CTRLv0qrhSATqdz5000_x16'
%}


%%qro
%'CTRLv0qro0qrhSATqdz5000'
%'CTRLv0qro25000qrhSATqdz5000'
%'CTRLv0qro50000qrhSATqdz5000'
%'CTRLv0qro100000qrhSATqdz5000'
%'CTRLv0qrhSATqdz5000'
%'CTRLv0qro400000qrhSATqdz5000'
%'CTRLv0qro800000qrhSATqdz5000'
%}

%{
%%ro
'CTRLv0ro0qrh0qdz5000'
'CTRLro25000v12.5qrh0qdz5000'
'CTRLro50000v12.5qrh0qdz5000'
'CTRLro100000v12.5qrh0qdz5000'
'CTRLro200000v12.5qrh0qdz5000'
'CTRLv12.5qrh0qdz5000'
'CTRLro800000v12.5qrh0qdz5000'
%}

%{
%%qro + ro
'CTRLv12.5qrhSATqdz5000'
'CTRLro800000v12.5qro50000qrhSATqdz5000'
'CTRLro50000v12.5qro800000qrhSATqdz5000'
%}




%{
'CTRLv0qrhSATqdz5000_lh3000'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_lh750'
%}

%{
'CTRLv0qrhSATqdz5000_dx1'
'CTRLv0qrhSATqdz5000_dx2'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_dx8'
'CTRLv0qrhSATqdz5000_dx16'
%}

%{
'CTRLv0qrhSATqdz5000_SST300K_Tthresh200K'
'CTRLv0qrhSATqdz5000_SST300K_Tthresh250K'
'CTRLv0qrhSATqdz5000_SST300K_Tthresh150K_nz48_zd25'
%}


%{
'CTRLv0qrhSATqdz5000_x8_dx8_ts16_lh0_SST300_200days'
'CTRLv0qrhSATqdz5000_x8_dx8_ts16_lh0_SST300_200days_fx2'
%}

%{
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_x2'
'CTRLv0qrhSATqdz5000_x4_dx8'
'CTRLv0qrhSATqdz5000_x8_dx8'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_lh12000_200days'
%'CTRLv0qrhSATqdz5000_x64_dx8_ts16_lh12000_200days_short'
%'CTRLv0qrhSATqdz5000_div8_dx.25_nx768_short'
%}

%{
'CTRLv0qrhSATqdz5000_x16'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_200days'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_lh12000_200days'
%}

%{
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_200days'
'CTRLv0qrhSATqdz5000_x32_dx8_ts16_200days_temp'
'CTRLv0qrhSATqdz5000_x64_dx8_ts16_200days_temp'
%}

%{
%'CTRLv0qrhSATqdz5000'
%'CTRLv0qrhSATqdz5000_x16'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_fdiv2'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16'
'CTRLv0qrhSATqdz5000_x16_dx8_ts16_fx2'
%}

%{
'CTRLv0qrhSATqdz5000_RST100days'
'CTRLv0qrhSATqdz5000_RSTflux.5_1300_1400'
%'CTRLv0qrhSATqdz5000_RSTflux.5_500_600'
%'CTRLv0qrhSATqdz5000_RSTflux.5_400_500'
'CTRLv0qrhSATqdz5000_RSTflux.5_300_400'
%'CTRLv0qrhSATqdz5000_RSTflux.5_200_300'
%'CTRLv0qrhSATqdz5000_RSTflux.5_100_200'
'CTRLv0qrhSATqdz5000_RSTflux.5_50_150'
%'CTRLv0qrhSATqdz5000_RSTflux.5_25_75'
'CTRLv0qrhSATqdz5000_RSTflux.5_12.5_37.5'
%'CTRLv0qrhSATqdz5000_RSTdiv2'
%}

%'CTRLv0qrhSATqdz5000_Lx2'
%'CTRLv0qrhSATqdz5000'
%'CTRLv0qrhSATqdz5000_Ldiv2'
%'CTRLv0qrhSATqdz5000_Ldiv4'
%}


%{
%'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_div2'
'CTRLv0qrhSATqdz5000_div4'
'CTRLv0qrhSATqdz5000_div8'
'CTRLv0qrhSATqdz5000_div16'
%}

%'CTRLv0qrhSATqdz5000_div2_ptype3'
%'CTRLv0qrhSATqdz5000_div2_ptype5'
%'CTRLv0qrhSATqdz5000_div2_ptype6'

%{
'CTRLv0qrhSATqdz5000_long'
'CTRLv0qrhSATqdz5000_RSTflux.5_50_150'
'CTRLv0qrhSATqdz5000_RSTflux.5_100_200'
'CTRLv0qrhSATqdz5000_RSTflux.5_150_250'
'CTRLv0qrhSATqdz5000_RSTflux.5_200_300'
'CTRLv0qrhSATqdz5000_RSTflux.5_250_350'
'CTRLv0qrhSATqdz5000_RSTflux.5_300_400'
'CTRLv0qrhSATqdz5000_RSTflux.5_400_500'
'CTRLv0qrhSATqdz5000_RSTflux.5_500_600'
'CTRLv0qrhSATqdz5000_RSTflux.5_600_700'
%}

%{
'CTRLv0qrhSATqdz5000_RSTdiv2'
'CTRLv0qrhSATqdz5000_long'
'CTRLv0qrhSATqdz5000_div2_RSTdiv4'
'CTRLv0qrhSATqdz5000_div2'
'CTRLv0qrhSATqdz5000_div4_RSTdiv8'
'CTRLv0qrhSATqdz5000_div4'
'CTRLv0qrhSATqdz5000_div8_RSTdiv16'
'CTRLv0qrhSATqdz5000_div8'
'CTRLv0qrhSATqdz5000_div16_RSTdiv32'
'CTRLv0qrhSATqdz5000_div16_short'
%}

%{
'CTRLv0qrhSATqdz5000_lv200'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_lv50'
'CTRLv0qrhSATqdz5000_lv25'
%}

%{
'CTRLv0qrhSATqdz5000_dx2'
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_dx8'
'CTRLv0qrhSATqdz5000_dx16'
%}

%{
'CTRLv0qrhSATqdz5000_x2'
'CTRLv0qrhSATqdz5000_x4_dx8'
'CTRLv0qrhSATqdz5000_x8_dx8'
%'CTRLv0qrhSATqdz5000_x16_dx8' -- FAIL
%}


%{
'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_lh750'
'CTRLv0qrhSATqdz5000_L1536dx2'
'CTRLv0qro100000qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_lh750_L1536dx2'
'CTRLv0qro100000qrhSATqdz5000_lh750'
'CTRLv0qro100000qrhSATqdz5000_L1536dx2'
'CTRLv0qrhSATqdz5000_div2'
%}

%{
%'CTRLv0qrhSATqdz5000_div16'
%'CTRLv0qrhSATqdz5000_div8'
%'CTRLv0qrhSATqdz5000_div4'
%'CTRLv0qrhSATqdz5000_div2'
%'CTRLv0qrhSATqdz5000'
'CTRLv0qrhSATqdz5000_x2'
%'CTRLv0qrhSATqdz5000_x4'
%'CTRLv0qrhSATqdz5000_x8'
%'CTRLv0qrhSATqdz5000_x16'
%}


}; %name of sub-directory with nc files

    
numruns=length(subdirs);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    clear r_usr Vmax
%{
    if(rr==1)
        t0 = 100;    %[day], starting time for averaging
        tf = 150;   %[day], ending time for averaging
    else
        t0 = 50;
        tf = 65;
    end
%}    
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
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_t0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    elseif(run_type==1)
        numfiles=length(dir(sprintf('%s/cm1out_0*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
    end
    %xres = dx; %this is done farther down now
    dt = time; %[s]; time of cm1out_t0001.nc is defined as zero
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
    i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

    %%TIME IN DAYS
%    t_min_day = dt*i_t0/86400;
%    t_max_day = dt*i_tf/86400;

    %%initialize the output vectors
    clear v_df_qv v_units_qv rprof_tmean
    
    clear t_day
    clear r_usr
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

        %%EXTRACT DATA FOR PHASE SPACE DIAGRAM
        var = 'vinterp';

        %%Extract Vmax
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
        xres = dx;
        Vmax(ii) = max(max(data));
        
        %%Extract r_usr = r(V_usr)
        i_max = find(data==max(max(data)));  %index of max wind speed
        i_max = i_max(1);
        if(strmatch(v_usr,'MAX'))
            if(~isempty(i_max))
                r_usr(ii)=xres*i_max-.5*xres;
            else
                r_usr(ii)=NaN;
            end
        else
            if(strmatch(var,'vinterp'))
                temp1=find(data<=v_usr);
                temp2=find(temp1>i_max,1);
                i_vusr_out=temp1(temp2);
                clear temp1 temp2
                if(~isempty(i_vusr_out))
                    v_out = data(i_vusr_out);   %v at first point where v<=v_usr
                    v_in = data(i_vusr_out-1);  %v at last point where v>=v_usr
                    r_usr(ii) = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                        if(r_usr(ii)<0 || r_usr(ii)>1000000 || Vmax(ii)<20)
                            r_usr(ii)=NaN;
                        end
                else
                    v_out = NaN;
                    v_in = NaN;
                    r_usr(ii) = NaN;
                end
            end
        end                

        
    end
    
    %% PLOT
    %%PLOT THINGS

    %%Use data at reasonable temporal resolution for viewing
    ptperday=1;
    i_plot = find(mod(t_day,1/ptperday)==0);
    
    r_usr = r_usr(i_plot);
    Vmax = Vmax(i_plot);
    t_day = t_day(i_plot);
    
    r_usr_save = r_usr;
    Vmax_save = Vmax;
    t_day_save = t_day;
    
    Vmaxmax=max(Vmax);
    r_usrmax=max(r_usr);
   
    %%Remove successive points that are too close dynamically
    Vmax_last=Vmax(1);
    r_usr_last=r_usr(1);
    for i=2:length(r_usr)
        dV=abs(Vmax(i)-Vmax_last);
        dr=abs(r_usr(i)-r_usr_last);
        
        if(dV>Vmaxmax/10 | dr>r_usrmax/10)
            Vmax_last=Vmax(i);
            r_usr_last=r_usr(i);
        else    %not enough change in data to shift to new point in time evolution plot
            r_usr(i)=NaN;
            Vmax(i)=NaN;
            t_day(i)=NaN;
        end
        
    end
    
    %remove the NaNs
    Vmax=Vmax(r_usr>0);
    t_day=t_day(r_usr>0);
    r_usr=r_usr(r_usr>0);
    
    %%set up for plot
    for i=1:10
        ii_plot_temp=find(t_day<=10*i);
        ii_plot(i)=ii_plot_temp(end);   %index of last entry below day 10*i
    end
            
    colors = jet(250);
    nc= 250/10;
    line_types={'-','--',':','-.'};
    
    if(plot_normalized==1)
        figure(1)
        plot(xvals/rmaxs(rr),rprof_tmean'/Vmaxs(rr),pl_clrs{rr},'LineWidth',2)
    %    axis([xmin_sub xmax_sub 0 Vmaxs(rr)])
        %axis([0/rmaxs(rr) rmax_plot/rmaxs(rr) 0 1])
        axis([0/rmaxs(rr) 20 0 1])
        hold on

    else        
        figure(1)
        %scatter(r_usr,Vmax,4,t_day,'filled')
        hold on
        h(rr)=plot(r_usr(1:ii_plot(1)),Vmax(1:ii_plot(1)),line_types{rr},'LineWidth',1.5,'Color',colors(nc,:));
        plot(r_usr(ii_plot(1):ii_plot(2)),Vmax(ii_plot(1):ii_plot(2)),line_types{rr},'LineWidth',1.5,'Color',colors(2*nc,:))
        plot(r_usr(ii_plot(2):ii_plot(3)),Vmax(ii_plot(2):ii_plot(3)),line_types{rr},'LineWidth',1.5,'Color',colors(3*nc,:))
        plot(r_usr(ii_plot(3):ii_plot(4)),Vmax(ii_plot(3):ii_plot(4)),line_types{rr},'LineWidth',1.5,'Color',colors(4*nc,:))
        plot(r_usr(ii_plot(4):ii_plot(5)),Vmax(ii_plot(4):ii_plot(5)),line_types{rr},'LineWidth',1.5,'Color',colors(5*nc,:))
        plot(r_usr(ii_plot(5):ii_plot(6)),Vmax(ii_plot(5):ii_plot(6)),line_types{rr},'LineWidth',1.5,'Color',colors(6*nc,:))
        plot(r_usr(ii_plot(6):ii_plot(7)),Vmax(ii_plot(6):ii_plot(7)),line_types{rr},'LineWidth',1.5,'Color',colors(7*nc,:))
        plot(r_usr(ii_plot(7):ii_plot(8)),Vmax(ii_plot(7):ii_plot(8)),line_types{rr},'LineWidth',1.5,'Color',colors(8*nc,:))
        plot(r_usr(ii_plot(8):ii_plot(9)),Vmax(ii_plot(8):ii_plot(9)),line_types{rr},'LineWidth',1.5,'Color',colors(9*nc,:))
        plot(r_usr(ii_plot(9):ii_plot(10)),Vmax(ii_plot(9):ii_plot(10)),line_types{rr},'LineWidth',1.5,'Color',colors(10*nc,:))
        %scatter(r_usr_save,Vmax_save,20,t_day_save,'filled')
        scatter(r_usr,Vmax,20,t_day,'filled')
        colormap(jet)
        colorbar
        axis([rmin_plot rmax_plot vmin_plot vmax_plot])

        
    end
            
end

figure(1)
set(gca,'fontweight','bold','fontsize',11)
if(plot_normalized==1)
    input_xlabel=sprintf('r / rmax');
    xlabel(input_xlabel);
    input_ylabel=sprintf('V / Vmax',v_def);
    ylabel(input_ylabel);
    input_title1=sprintf('Normalized azimuthal velocity');
    input_title2=sprintf('z=%5.2f %s',zmin_sub,zunits);
    title({input_title1,input_title2})
else
    input_xlabel=sprintf('R%i [%s]',v_usr,'km');
    xlabel(input_xlabel);
    input_ylabel=sprintf('Vmax [%s]','ms-1');
    ylabel(input_ylabel);
    input_title1=sprintf('Time evolution of TC structure [days]');
    title(input_title1)
end
input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
legend(h,input_legend)

    
%}