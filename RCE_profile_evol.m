%RCE_profile_evol.m

%Created: 5 Apr 2012, Dan Chavas

%Purpose: to plot the mean vertical profile of qv and theta over successive 10-day
%periods to view equilibration of the atmosphere.

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc
close('all')

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='RCE/';
%subdir_pre='CTRL_icRCE/';
%subdir_pre='TRANSFER/';
ext_hd = 1; %0=local hard drive; 1='CHAVAS_CM1_FINAL'; 2='CHAVAS_CM1_FINAL_nodrag'

%SST = 300.00;  %[K]; EXTRACTED BELOW used to calculate RCE th_sfc (=SST-2K) and qv_sfc (80% RH from saturation at SST)
RH_sfc = .8;    %air-sea latent disequilibrium
dT_sfc = 2; %[K]; air-sea thermal disequilibrium

run_types=3*ones(100,1);    %[1 1 1 1 1 1 1 1 1]; %1=axisym; 3=3d
subdirs = {
'RCE_nx48_SST285.00K_Tthresh200K_usfc3_dz312.5_drag'
}; %name of sub-directory with nc files

t0a = 50;    %[day], starting time for averaging
tfa = 100;   %[day], ending time for averaging
interval_day = 5;  %[day], length of intervals

save_output_sounding = 0;   %0=no output file created; 1=yes 'input_sounding_[subdir]'
plot_type = 1;  %0=no plot; 1=plots of RCE vertical profiles of qv [g/kg] and theta [K]
    %pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
    
x0=0;   %first x grid point [0,end]
xf=10000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=10000;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=100;  %last z grid point [0,end]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numruns=length(subdirs);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)


%%Append ax or 3d to subdir names for plot
%subdirs_plot = subdirs;

for i=1:length(subdirs)
    if(run_types(i)==1) %ax
%        subdirs_plot{2*i-1}=sprintf('ax%s',subdirs{i})
%        subdirs_plot{2*i}=sprintf('ax%s_init',subdirs{i});
        subdirs_plot{i}=sprintf('ax%s',subdirs{i})
    else    %3d
%        subdirs_plot{2*i-1}=sprintf('3d%s',subdirs{i})
%        subdirs_plot{2*i}=sprintf('3d%s_init',subdirs{i});
        subdirs_plot{i}=sprintf('3d%s',subdirs{i})
    end
end
%}

%% DETERMINE IF A QUASI-STEADY STATE EXISTS, AND IF SO, WHEN
for rr=1:numruns
ts = t0a:interval_day:tfa;
for pp=1:length(ts)-1    
    t0 = ts(pp);
    tf = ts(pp+1);
    %%EXTRACT RUN_TYPE AND SUBDIR NAME
    run_type=run_types(rr);
    subdir=subdirs{rr};

    %% Extract T_sst from file name
    i_sst=strfind(subdir,'SST');
    if(isempty(i_sst))
        sst_str='300.00';
    else
        i_sstK=strfind(subdir,'K');
        i_sstK=i_sstK(find(i_sstK>i_sst,1));
        sst_str=subdir(i_sst+3:i_sstK-1);
    end
    SST = str2num(sst_str);
    
    %% OPTIONS FOR EITHER AXISYM OR 3d RUNS
    if(run_type==1)
        run_type_str='axisym';
        copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
    elseif(run_type==3)
        run_type_str='3D';
        copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
    end

    switch ext_hd
        case 0,
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/%s/CM1_output/%s',run_type_str,subdir_pre);
        case 1,
            dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL/CM1_output/%s/%s',run_type_str,subdir_pre);
        case 2,
            dir_in=sprintf('/Volumes/CHAVAS_CM1_FINAL_nodrag/CM1_output/%s/%s',run_type_str,subdir_pre);
        otherwise
            assert(2==3,'Invalid number for ext_hd!')
    end   
    
    %%DIRECTORY WITH OUTPUT DATA
    subdir_full=sprintf('%s%s',dir_in,subdir)
    
    %% EXTRACT DATA FROM ORIGINAL nc DATA

    %%EXTRACT TIMESTEP SIZE
    var = 'qvpert'; %doesn't matter, just need any variable
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
    if(run_type==3)
        numfiles=length(dir(sprintf('%s/cm1out_*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var,x0,xf,y0,yf,z0,zf);
    elseif(run_type==1)
        numfiles=length(dir(sprintf('%s/cm1out_*.nc',subdir_full)));
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var,x0,xf,y0,yf,z0,zf);
    end
    dt = time; %[s]; time of cm1out_0001.nc is defined as zero
    clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy nx_sub ny_sub xunits yunits zunits v_def v_units time t_units
    i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
    i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

    %% EXTRACT DATA FROM input_sounding
    snd_file = 'input_sounding';
    
    %%load initial sounding variables
    dir_full = strcat(dir_in,subdir);
    [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

    
    %% LOOP OVER FILES TO OBTAIN TIME SERIES DATA
    %%initialize the output vectors arbitrarily large
    clear data_hmean_qv data_hmean_th data_hmean_p
    data_hmean_qv=zeros(1000,1);
    data_hmean_th=zeros(1000,1);
    data_hmean_p=zeros(1000,1);

    for ii=1:i_tf-i_t0+1

        t_file=i_t0+ii-1;

        %%NC FILENAME (default)
        if(run_type==1)
            if(t_file<10)
                nc_file = sprintf('cm1out_000%i.nc',t_file);
            elseif(t_file<100)
                nc_file = sprintf('cm1out_00%i.nc',t_file);
            else
                nc_file = sprintf('cm1out_0%i.nc',t_file);
            end
        elseif(run_type==3)
            if(t_file<10)
                nc_file = sprintf('cm1out_000%i.nc',t_file);
            elseif(t_file<100)
                nc_file = sprintf('cm1out_00%i.nc',t_file);
            else
                nc_file = sprintf('cm1out_0%i.nc',t_file);
            end    
        end

        %%EXTRACT qv DATA: note that qvpert is in [kg/kg]
        var = 'qvpert'
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

        %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
        v_def_qv = v_def;
        v_units_qv = v_units;
        data_hmean_qv=data_hmean_qv(1:nz_sub);   %NOTE qvpert IS IN [kg/kg]!!!
        if(run_type==3)
            vert_prof_temp_qv=squeeze(mean(mean(data,1),2));   %horiz-averaged vertical profile
        elseif(run_type==1)
            xvals = xmin_sub:dx:xmax_sub;
            data=squeeze(data);
            vert_prof_temp_qv = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius
        end
        data_hmean_qv=data_hmean_qv+vert_prof_temp_qv/(i_tf-i_t0+1); %[kg/kg]
        
        %%EXTRACT theta DATA
        var = 'thpert'
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);
        zvals = zmin_sub:dz:zmax_sub;
        
        %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
        v_def_th = v_def;
        v_units_th = v_units;
        data_hmean_th=data_hmean_th(1:nz_sub);
        if(run_type==3)
            vert_prof_temp_th=squeeze(mean(mean(data,1),2));   %horiz-averaged vertical profile
        elseif(run_type==1)
            xvals = xmin_sub:dx:xmax_sub;
            data=squeeze(data);
            vert_prof_temp_th = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius
        end
        data_hmean_th=data_hmean_th+vert_prof_temp_th/(i_tf-i_t0+1); %[g/kg]
        
        %%EXTRACT pressure DATA
        var_temp = 'prspert'
%        var_temp = 'pi'
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
        [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
        
        %%CALCULATE HORIZONTALLY-AVERAGED VERTICAL PROFILE
        v_def_p = v_def;
        v_units_p = v_units;
        data_hmean_p=data_hmean_p(1:nz_sub);
        if(run_type==3)
            vert_prof_temp_p=squeeze(mean(mean(data,1),2));   %horiz-averaged vertical profile
        elseif(run_type==1)
            xvals = xmin_sub:dx:xmax_sub;
            data=squeeze(data);
            vert_prof_temp_p = ((xvals*data)/sum(xvals))';  %cylindrical integral--sum weighted by radius
        end
        data_hmean_p=data_hmean_p+vert_prof_temp_p/(i_tf-i_t0+1); %[g/kg]        
        pi = data_hmean_p;
        
        clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub xunits yunits zunits v_def v_units time t_units
        
    end
    
    vec_length = min(nz_sub,length(zz00)-z0);
          
    zz00_RCE{rr} = zz00(z0+1:z0+vec_length)'; %heights [m]
    qv00_RCE{rr} = qv00(z0+1:z0+vec_length)' + [data_hmean_qv(z0+1:z0+vec_length)];   %[kg/kg]
    th00_RCE{rr} = th00(z0+1:z0+vec_length)' + [data_hmean_th(z0+1:z0+vec_length)];   %[K] 
    p00_RCE{rr} = pp00(z0+1:z0+vec_length)' + [data_hmean_p(z0+1:z0+vec_length)];   %[Pa]
    
    %%Fix any vertical instabilities in th00_RCE (just set value to mean of levels it lies between)
    %%MOIST ONLY! For dry case, this needs to be updated since dth/dz<0 in lowest several model levels
%TEST    th00_RCE{rr}(10)=th00_RCE{rr}(9)-.001;
%TEST    th00_RCE{rr}(9)=th00_RCE{rr}(8)-.01;

        RCE_instab_remove = 0;
        instab_check = th00_RCE{rr}(2:end)-th00_RCE{rr}(1:end-1);
        max_instab = min(instab_check);
        z_instab = find(instab_check == max_instab);

        if(max_instab<0)    %there is an instability to remove

            %assert(abs(max_instab)<.1,'WARNING: THERE IS A LARGE INSTABILITY IN MOIST RCE PROFILE')

            RCE_instab_remove = 1;   %the profile will be adjusted
            
            %%Iterate upwards through the RCE profile and apply mass-weighted averaging ("mixing") over
            %%each statically-unstable layer and repeat until the entire profile is stable
            th00_RCE_adj{rr} = th00_RCE{1};  %adjusted profile begins as true RCE profile
            dth = th00_RCE_adj{rr}(2:end)-th00_RCE_adj{rr}(1:end-1);  %change in theta with height

            while(sum(dth<0)>0)   %there are unstable layers
            for j = 1:length(th00_RCE_adj{rr})-1
                dth_j = th00_RCE_adj{rr}(j+1)-th00_RCE_adj{rr}(j);
                if(dth_j<0)   %unstable
                    th_ave = (p00_RCE{rr}(j)*th00_RCE_adj{rr}(j) + p00_RCE{rr}(j+1)*th00_RCE_adj{rr}(j+1))/(p00_RCE{rr}(j)+p00_RCE{rr}(j+1));
                    th00_RCE_adj{rr}(j) = th_ave;
                    th00_RCE_adj{rr}(j+1) = th_ave;
                end
            end
            
            %%NOTE: I do NOT mix qv as well, as this is not a true
            %%instability but instead a result of numerical noise, unlike
            %%in the dry case (though in that case, qv=0 anyways)

            %update change in theta with height
            dth = th00_RCE_adj{rr}(2:end)-th00_RCE_adj{rr}(1:end-1);  %change in theta with height

            end
            
            th00_RCE{rr} = th00_RCE_adj{rr};    %only need the final adjusted profile
            
        end
%}
    
    %%Recalculate sfc theta (=SST-2K) and qv_sfc (80% RH from saturation at SST)
    %% Constants (values taken from CM1 model)
    c_CM1 = constants_CM1(); %c_CM1: [g rd cp cv p00 xlv]

    g=c_CM1(1); %[m/s2]
    Rd=c_CM1(2);  %[J/kg/K]
    Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
    Rv=c_CM1(4);   %[J/K/kg]
    p0 = c_CM1(5); %[Pa]
    Lv=c_CM1(6);   %[J/kg]

    eps=Rd/Rv;
    
    SST_C = SST-273.15;
    es_SST = 6.112*exp((17.67*SST_C)./(SST_C+243.5))*100;  %[Pa]
    qvs_SST = eps*(es_SST./(p_sfc-es_SST));  %[kg/kg]
    
    qv_sfc_RCE(rr) = RH_sfc*qvs_SST;  %80% RH
    
    th_sfc_RCE(rr) = SST-dT_sfc;  %2K disequilibrium

    %%Save input soundings for plotting
    qv_sfc_init(rr) = qv_sfc;
    th_sfc_init(rr) = th_sfc;
    zz00_init{rr} = zz00(z0+1:z0+vec_length)'; %heights [m]
    qv00_init{rr} = qv00(z0+1:z0+vec_length)';   %[kg/kg]
    th00_init{rr} = th00(z0+1:z0+vec_length)';   %[K] 


if(save_output_sounding~=0)

    %%sounding header: %1) p_sfc (mb); 2) th_sfc (K); 3) qv_sfc (g/kg)
    %%sounding columns: 1) zz00 (m); 2) th00_RCE (K); 3) qv00_RCE (g/kg); 4) u00 (m/s); 5) v00 (m/s)
    
    %%need extra layers at top above model top 
    zz00_ext = zz00_RCE{1}(end)+(zz00_RCE{1}(end)-zz00_RCE{1}(end-1))*[1 2 3];
    th00_ext = th00_RCE{1}(end)+(th00_RCE{1}(end)-th00_RCE{1}(end-1))*[1 2 3];
    qv00_ext = qv00_RCE{1}(end)*[1 1 1];
    u00_ext = u00(end)*[1 1 1];
    v00_ext = v00(end)*[1 1 1];
    
    
    sounding_all=[[zz00_RCE{1};zz00_ext'] [th00_RCE{1};th00_ext'] 1000*[qv00_RCE{1};qv00_ext'] [u00(1:vec_length)';u00_ext'] [v00(1:vec_length)';v00_ext']];
    sounding_all=double(sounding_all);
    
    pressure_all=[p00_RCE{1}];
    pressure_all=double(pressure_all);

    file_out = sprintf('input_sounding_%s',subdirs_plot{1});
%    save(file_out,'p_sfc','th_sfc','qv_sfc','sounding_all');    
    clear fid
    fid = fopen(file_out,'w');
    fprintf(fid,'   %6.2f     %6.2f     %6.2f\n',p_sfc/100,th00_RCE{1}(1),1000*qv00_RCE{1}(1));
    fprintf(fid,'   %8.5f     %8.5f     %8.5f     %8.5f    %8.5f\n',sounding_all');
    
    file_out_pres = sprintf('input_sounding_%s_pressures',subdirs_plot{1});
%    save(file_out,'p_sfc','th_sfc','qv_sfc','sounding_all');    
    clear fid
    fid = fopen(file_out_pres,'w');
    fprintf(fid,'   %8.5f     \n',pressure_all');


end


%%Make plot
if(plot_type~=0)
    
    %clf     %clear figures without closing them
    
    %% LOOP OVER RUNS TO INCUDE IN FIGURE
    for rr=1:numruns

        %%EXTRACT RUN_TYPE AND SUBDIR NAME
        run_type=run_types(rr);
        subdir=subdirs{rr};

        %% OPTIONS FOR EITHER AXISYM OR 3d RUNS
        if(run_type==1)
            run_type_str='axisym';
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/axisym/CM1_output/%s',subdir_pre);
            copyfile('nc_extract_axisym.m','nc_extract.m') %copy nc_extract_axisym.m to nc_extract
        elseif(run_type==3)
            run_type_str='3d';
            dir_in=sprintf('/Users/drchavas/Documents/Research/Thesis/CM1/v15/3d/CM1_output/%s',subdir_pre);
            copyfile('nc_extract_3d.m','nc_extract.m') %copy nc_extract_3d.m to nc_extract
        end

        %%PLOT VERTICAL PROFILES
        %PLOT DATA
        %vec_length=min(sum(~isnan(qv00_RCE(:,rr))),length(zz00));
        num_dts=length(ts)-1;
        pl_clr = (1-(pp/num_dts))*[0 1 0] + (pp/num_dts)*[1 0 0];

        figure(1)
        subplot(1,2,1)
        hold on
        %plot(1000*[qv_sfc_RCE(rr);qv00_RCE{rr}],[0;zz00_RCE{rr}/1000],pl_clrs{2*rr-1},'LineWidth',1)
        plot(1000*[qv00_RCE{rr}],[zz00_RCE{rr}/1000],'Color',pl_clr,'LineWidth',1)

            set(gca,'fontweight','bold','fontsize',11)

        figure(1)
        subplot(1,2,2)
        hold on
        %plot([th_sfc_RCE(rr);th00_RCE{rr}],[0;zz00_RCE{rr}/1000],pl_clrs{2*rr-1},'LineWidth',1)
        plot([th00_RCE{rr}],[zz00_RCE{rr}/1000],'Color',pl_clr,'LineWidth',1)
th00_all(pp,:) = th00_RCE{rr};
            set(gca,'fontweight','bold','fontsize',11)

        %%Plot initial soundings
%{
        figure(1)
        subplot(1,2,1)
        hold on
        plot(1000*[qv_sfc_init(rr);qv00_init{rr}],[0;zz00_init{rr}/1000],pl_clrs{2*rr},'LineWidth',1)

        figure(1)
        subplot(1,2,2)
        hold on
        plot([th_sfc_init(rr);th00_init{rr}],[0;zz00_init{rr}/1000],pl_clrs{2*rr},'LineWidth',1)
%}
        
    end

    hold on
    figure(1)
    subplot(1,2,1)
    set(gca,'fontweight','bold','fontsize',11)
    axis([0 1000*max(cellfun(@(x) max(x(:)),qv00_RCE)) 0 zz00(end)/1000])
    input_xlabel=sprintf('water vapor mixing ratio [g/kg]',v_def_qv);
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_legend=strrep([subdirs_plot],'_','\_');
    legend(input_legend)
    input_title=strrep(sprintf('water vapor mixing ratio [g/kg]'),'_','\_');
    title(input_title)

    figure(1)
    subplot(1,2,2)
    set(gca,'fontweight','bold','fontsize',11)
    axis([min(min(cellfun(@(x) min(x(:)),th00_RCE))) max(max(cellfun(@(x) max(x(:)),th00_RCE))) 0 zz00(end)/1000])
    input_xlabel=sprintf('potential temperature [%s]',v_units_th);
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_legend=strrep([subdirs_plot],'_','\_');
    legend(input_legend)
    input_title=strrep(sprintf('potential temperature [%s]',v_units_th),'_','\_');
    title(input_title)
 %}   
    input_title=strrep(sprintf('vert profs: horiz-mean days %i-%i',t0,tf),'_','\_');
    annotation('textbox','String',input_title,'Position',[.1 .5 .9 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')

%    title({input_title1,input_title2})

%{
    %STUFF FOR PLOT
    figure(1)
    set(gca,'fontweight','bold','fontsize',11)
    axis([0 20 0 25])
    input_xlabel=sprintf('%s [%s]',v_def_qv,v_units_qv);
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_title=strrep(sprintf('%s vertical profile: %s [%s]',subdir,v_def_qv,v_units_qv),'_','\_');
    title(input_title)
    input_legend=strrep(subdirs(1:numruns),'_','\_');
    legend(input_legend)

    figure(2)
    set(gca,'fontweight','bold','fontsize',11)
    axis([200 650 0 25])
    input_xlabel=sprintf('%s [%s]',v_def_qv,v_units_qv);
    xlabel(input_xlabel);
    input_ylabel=sprintf('Height AGL [km]');
    ylabel(input_ylabel);
    input_title=strrep(sprintf('%s vertical profile: %s [%s]',subdir,v_def_th,v_units_th),'_','\_');
    title(input_title)
    input_legend=strrep(subdirs(1:numruns),'_','\_');
    legend(input_legend)
%}

end

end

end



%}