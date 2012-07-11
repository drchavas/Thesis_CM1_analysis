%nc_plot_steady_xz.m

%Created: 10-05-11, Dan Chavas

%Purpose: to plot quasi-steady state xz

%NOTE: data at varying x (i.e. zonal cross-sec) = matrix COLUMN in MATLAB!

clear
clc
figure(1)
hold on
figure(2)
clf(1)     %clear figures without closing them
clf(2)     %clear figures without closing them


%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

run_type=1; %1=axisym; 3=3D

t0 = 70;
tf = 100;

%mult_by = [8 4 2 1 1/2 1/4 1/8];

plot_type = 1;  %0=no plot
                %1=single plot of time series of [var] for above runs
                %4=4-panel:
                    %axisym: vinterp, w, qv, uinterp
                    %3D: vinterp, mid-level w, qv, XZ-xsec qv
%IF plot_type=1: input variable and domain
%    pl_clrs={'b--' 'k' 'r--' 'b' 'g--' 'r' 'c--' 'g' 'y--' 'c' 'y--' 'y'};
    %pl_clrs={'b' 'b--' 'r' 'r--' 'g' 'g--' 'c' 'c--' 'k' 'k--' 'y' 'y--'};
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--'};
%    pl_clrs={'b--' 'k--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'b' 'k' 'r' 'g' 'c' 'y' 'm'};
    %pl_clrs={'g' 'c' 'k' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'y--'};
%    pl_clrs={'r' 'g' 'c' 'k' 'y' 'm' 'r--' 'g--' 'c--' 'k--' 'y--'};
    vars = {'angmom'}%'vinterp' 'Ri'};  %pcolor; contour on top
        %IF vars = 'vinterp' or 'vg'
        r_ts = 0; %1=plot time series of r(v_usr) 0=don't
        v_usr = 7; %find radius of this wind speed; 'MAX' for V_max
        numpt_sm = 1;  %smoother uses this number of points

    
%%Define grid of points to extract (i.e. all of them, will subset afterwards)
x0=0;   %first x grid point [0,end]
xf=10000;   %first y grid point [0,end]
y0=0;   %first y grid point [0,end]
yf=0;   %last y grid point [0,end]
z0=0;  %first z grid point [0,end]
zf=100;  %last z grid point [0,end]

%%Define subset of points to plot
rmin_plot = 0;  %[km]; lowest value plotted
rmax_plot = 4000;    %[km]; highest value plotted
zmin_plot = 0;  %[km]; lowest value plotted
zmax_plot = 20; %[km]; highest value plotted
datamin_plot = -5000;   %minimum data value plotted
datamax_plot = 5000;   %minimum data value plotted

subdir = 'CTRLv0qrhSATqdz5000_nx3072'; %name of sub-directory with nc files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numvars = length(vars);

%%Append ax or 3d to subdir names for plot
if(run_type==1) %ax
    subdir_plot=sprintf('ax%s',subdir)
else    %3d
    subdir_plot=sprintf('3d%s',subdir)
end


%% CALCULATE STEADY STATE RADIAL PROFILE
t_day_min = 0;
t_day_max = 0;

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
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_t0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
elseif(run_type==1)
    numfiles=length(dir(sprintf('%s/cm1out_0*.nc',subdir_full)));
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units fcor] = nc_extract(dir_in,subdir,'cm1out_0002.nc',var_dt,x0,xf,y0,yf,z0,zf);
end

dt = time; %[s]; time of cm1out_t0001.nc is defined as zero
clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy nx_sub ny_sub xunits yunits zunits v_def v_units time t_units
i_t0 = max(0,round(t0*24*60*60/dt)+1); %timestep corresponding to start time for averaging
i_tf = min(numfiles,round(tf*24*60*60/dt)+1); %timestep corresponding to end time for averaging

%%TIME IN DAYS
%    t_min_day = dt*i_t0/86400;
%    t_max_day = dt*i_tf/86400;

%%DIRECTORY WITH OUTPUT DATA
subdir_full=sprintf('%s%s',dir_in,subdir)

%% EXTRACT DATA FROM input_sounding
snd_file = 'input_sounding';

%%load initial sounding variables
dir_full = strcat(dir_in,subdir);
[zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);

%%initialize the output vectors
clear v_df_qv v_units_qv data_tmean
for vv=1:numvars
    data_tmean{vv}=zeros(xf-x0+1,zf-z0+1);
end


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
 
    %% Calculate various thermodynamic quantities
    %%load prspert
    var_temp = 'prspert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    pp_pert = squeeze(data);

    %%load thpert
    var_temp = 'thpert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    th_pert = squeeze(data);

    %%load qvpert
    var_temp = 'qvpert';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    qv_pert = squeeze(data);

    %%load upert
    var_temp = 'uinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    u_pert = squeeze(data);

    %%load vpert
    var_temp = 'vinterp';
    [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var_temp,x0,xf,y0,yf,z0,zf);
    v_pert = squeeze(data);

    %%calculate pressure
    p = repmat(pp00(z0+1:z0+nz_sub),nx_sub,1) + pp_pert;

    %%calculate potential temperature
    th = repmat(th00(z0+1:z0+nz_sub),nx_sub,1) + th_pert;

    %%calculate water vapor mixing ratio
    qv = repmat(qv00(z0+1:z0+nz_sub),nx_sub,1) + qv_pert;    %[kg/kg]

    %%calculate T, Tv, rho, rh
    [T rho rh the s s_sat gam_m] = thermo(p,th,qv);

    
    %%Calculate final variable data
    for vv = 1:numvars
        var = vars{vv};
        if(strcmp(var,'Ri'))    %composite variable: Richardson Number
            %%NOTE: size(s)=[384 40]; col 1 = lowest model level
            ds_dz = NaN*s;   %initialize as NaN
            ds_dz(:,2:end-1) = (s(:,3:end)-s(:,1:end-2))/(2000*dz); %[J/K/m]

            %%calculate wind speeds
            u = repmat(u00(z0+1:z0+nz_sub),nx_sub,1) + u_pert;
            v = repmat(v00(z0+1:z0+nz_sub),nx_sub,1) + v_pert;

            du_dz = (u(:,3:end)-u(:,1:end-2))/(2*1000*dz); %[s-1]
            du_dz(:,2:end+1)=du_dz;  %add extra level at bottom equal to nearest level value
            du_dz(:,end+1)=du_dz(:,end);  %add extra level at top equal to nearest level value
            dv_dz = (v(:,3:end)-v(:,1:end-2))/(2*1000*dz); %[s-1]
            dv_dz(:,2:end+1)=dv_dz;  %add extra level at bottom equal to nearest level value
            dv_dz(:,end+1)=dv_dz(:,end);  %add extra level at top equal to nearest level value

            data = gam_m.*ds_dz./(du_dz.^2 + dv_dz.^2);    %Richardson number, eqn (23) Emanuel11
            data(data>10) = NaN;
            data(data<0) = NaN;
            v_defs{vv} = 'Richardson Number';
            v_unitss{vv} = '-';

        elseif(strcmp(var,'angmom'))    %composite variable: angular momentum
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);
            
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1)';
            data = squeeze(data).*xmat*1000 + .5*fcor*(xmat*1000).^2;   %M=rV + .5fr^2
            %data = s;
            %v_defs{vv} = 'entropy';
            %v_unitss{vv} = '-';
            v_defs{vv} = 'Angular momentum';
            v_unitss{vv} = 'm^2/s^2';

        elseif(strcmp(var,'vg'))    %composite variables: gradient wind
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,'vinterp',x0,xf,y0,yf,z0,zf);
            
            xvals = xmin_sub:dx:xmax_sub;
            zvals = zmin_sub:dz:zmax_sub;
            xmat = repmat(xvals,length(zvals),1)';
            dp_dr = (p(3:end,:)-p(1:end-2,:))/(2*1000*dx); %[m/s^2]
            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value
            dp_dr(end+1,:)=dp_dr(end,:);  %add extra level at outer edge equal to nearest level value

            %calculate at mid-points between grid points
%            dp_dr = (p(2:end,:)-p(1:end-1,:))/(1000*dx); %[m/s^2]
%            dp_dr(2:end+1,:)=dp_dr;  %add extra level at center equal to nearest level value

            
            vinterp = data;
            temp2=(xmat*1000*fcor).^2+4*xmat*1000.*dp_dr;
            data =  .5*(-xmat*1000*fcor+sqrt((xmat*1000*fcor).^2+4*xmat*1000.*dp_dr));   %Vg^2 + r*f*Vg-r*(dp_dr); solve for Vg
            data = real(data);  %remove imaginary numbers, which are always very small
            %data(imag(data)~=0)=NaN;    %ignore imaginary numbers
%            data = data - squeeze(vinterp);
            v_defs{vv} = 'Gradient wind';
            v_unitss{vv} = 'm/s';
            
            
            
        else %just extract data for the input variable itself
            %%EXTRACT DATA
            clear data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units
            [data xmin_sub xmax_sub ymin_sub ymax_sub zmin_sub zmax_sub dx dy dz nx_sub ny_sub nz_sub xunits yunits zunits v_def v_units time t_units] = nc_extract(dir_in,subdir,nc_file,var,x0,xf,y0,yf,z0,zf);

            v_defs{vv} = v_def;
            v_unitss{vv} = v_units;
            
            %data = s;
        end

        %%STORE THE WIND SPEED AT THE DESIRED RADIUS
        if(r_ts==1)
            xres = dx;
            i_max = find(data==max(max(data)));  %index of max wind speed
            i_max = i_max(1);
            if(strmatch(v_usr,'MAX'))
                if(~isempty(i_max))
                    r_usr(ii)=xres*i_max-.5*xres;
                else
                    r_usr(ii)=NaN;
                end
            else
                if(strcmp(var,'vinterp') || strcmp(var,'vg'))
                    temp1=find(data<=v_usr);
                    temp2=find(temp1>i_max,1);
                    i_vusr_out=temp1(temp2);
                    clear temp1 temp2
                    if(~isempty(i_vusr_out))
                        v_out = data(i_vusr_out);   %v at first point where v<=v_usr
                        v_in = data(i_vusr_out-1);  %v at last point where v>=v_usr
                        r_usr(ii) = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                        if(r_usr(ii)<0 || r_usr(ii)>1000000 || max(max(data))<20)
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
        
        %%CALCULATE TIME-AVERAGED xz-cross-section
        data=squeeze(data);
        data_tmean{vv}=data_tmean{vv}(1:nx_sub,1:nz_sub);
        data_tmean{vv}=data_tmean{vv}+data/(i_tf-i_t0+1);
        
    end

        
end

   
%%Extract subset of data for plotting
%ALL DATA
xvals = xmin_sub:dx:xmax_sub;
zvals = zmin_sub:dz:zmax_sub;
i_xvals = max(1,x0+1):max(1,x0+1)+nx_sub;
i_zvals = max(1,z0+1):max(1,z0+1)+nz_sub;

%SUBSET DATA
indices = find(xvals>=rmin_plot & xvals<=rmax_plot);
xvals = xvals(indices);
i_xvals = i_xvals(indices);
clear indices
indices = find(zvals>=zmin_plot & zvals<=zmax_plot);
zvals = zvals(indices);
i_zvals = i_zvals(indices);
clear indices

for vv=1:numvars
    data_tmean{vv} = data_tmean{vv}(i_xvals',i_zvals');
    data_tmean{vv}=double(data_tmean{vv});  %pcolor requires double format
end

%%CALCULATE THE WIND SPEED AT THE DESIRED RADIUS FOR THE TIME-MEAN PROFILE
if(r_ts==1)
    xres = dx;
    i_max = find(data_tmean{vv}==max(max(data_tmean{vv})));  %index of max wind speed
    i_max = i_max(1);
    if(strmatch(v_usr,'MAX'))
        if(~isempty(i_max))
            r_usr_tmean=xres*i_max-.5*xres;
        else
            r_usr_tmean=NaN;
        end
    else
        if(strcmp(var,'vinterp') || strcmp(var,'vg'))
            temp1=find(data_tmean{vv}<=v_usr);
            temp2=find(temp1>i_max,1);
            i_vusr_out=temp1(temp2);
            clear temp1 temp2
            if(~isempty(i_vusr_out))
                v_out = data_tmean{vv}(i_vusr_out);   %v at first point where v<=v_usr
                v_in = data_tmean{vv}(i_vusr_out-1);  %v at last point where v>=v_usr
                r_usr_tmean = xres*(i_vusr_out-(v_out-v_usr)/(v_out-v_in)) - .5*xres;   %r=dr/2 @ i=1
                if(r_usr_tmean<0 || r_usr_tmean>100000000)
                    r_usr_tmean=NaN;
                end
            else
                v_out = NaN;
                v_in = NaN;
                r_usr_tmean = NaN;
            end
        end
    end
end

%% PLOT THINGS
%clear xmat zmat i_xvals i_zvals
xmat = repmat(xvals,length(zvals),1);
zmat = repmat(zvals,length(xvals),1)';

%%BASIC STATISTICS
for vv=1:numvars
    data_min(vv) = min(min(data_tmean{vv}));
    data_max(vv) = max(max(data_tmean{vv}));
end

%%BLOCK OUT UNDESIRED DATA VALUES
if(numvars > 1)
    data_tmean{2}(data_tmean{2}>datamax_plot)=datamax_plot;
    data_tmean{2}(data_tmean{2}<datamin_plot)=datamin_plot;
end

%1-DIM RADIAL PLOT
figure(1)

if(length(zvals)>1)   %rz color plot

    if(numvars > 1)
        pcolor(xmat,zmat,data_tmean{1}')
        shading flat
        shading interp
    else
        contourf(xmat,zmat,data_tmean{1}',20)
    end
    axis([xvals(1) xvals(end) zvals(1) zvals(end)])
    hold on
    temp=caxis;
    caxis(temp);    %are you serious?
    clear temp
    colorbar

    if(numvars>1)
        [C H] = contour(xmat,zmat,data_tmean{2}',5,'k');
        clabel(C,H)
        [C1 H1] = contour(xmat,zmat,data_tmean{2}',[1 1],'r');
        clabel(C1,H1)
    end

    %%for plot labeling
    input_ylabel=sprintf('z-distance [%s]',zunits);
    if(numvars>1)
        input_title2=sprintf('min=%5.3f max=%5.3f',data_min(2),data_max(2));
    end
    input_title3='';

else    %radial-only plot

    plot(xmat,data_tmean{1}','r--','LineWidth',2)
    %axis([rmin_plot rmax_plot datamin_plot datamax_plot])
    hold on

    %%for plot labeling
    input_ylabel=sprintf('%s [%s]',v_defs{1},v_unitss{1});
    
    input_title3=sprintf('z=%5.2f %s',zvals,zunits);

end
%[i j] = find(data_tmean{2}==min(min(data_tmean{2})));
%plot(i,j,'x')


figure(1)
plot(xvals,0*xvals,'k')
set(gca,'fontweight','bold','fontsize',11)
input_xlabel=sprintf('x-distance [%s]',xunits);
xlabel(input_xlabel);
ylabel(input_ylabel);
input_title1=sprintf('%s [%s]',v_defs{1},v_unitss{1});
input_title2=sprintf('min=%5.2f max=%5.2f',data_min(1),data_max(1));
title({input_title1,input_title2,input_title3})

%input_legend=strrep([subdirs_plot(1:numruns)],'_','\_');
%legend(input_legend)

%% FIGURE TITLE
input_title=strrep(sprintf('time-mean days %i-%i',t0,tf),'_','\_');
annotation('textbox','String',input_title,'Position',[0 .5 .2 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')
fig_title = sprintf('%s',subdir);
annotation('textbox','String',fig_title,'Position',[0 -.45 .5 .5],'Linestyle','none','Interpreter','none','HorizontalAlignment','center','FontSize',10,'Fontweight','bold')

%%TIME SERIES OF r_usr
if(r_ts==1)
    i_plot = ~isnan(r_usr);
    figure(2)
%    t_day_min = min(t_day);
%    t_day_max = max(t_day_max,max(t_day));

    plot(t_day(i_plot)',smooth(r_usr(i_plot),numpt_sm)','LineWidth',2)
%    axis([t_day_min t_day_max rmin_plot rmax_plot])
    hold on
    smooth_str=sprintf('%i-pt smooth',numpt_sm);
    r_usr_means = mean(r_usr);

    set(gca,'fontweight','bold','fontsize',11)
    input_xlabel=sprintf('time [day]');
    xlabel(input_xlabel);
    if(strcmp(v_usr,'MAX'))
        input_ylabel=sprintf('radius of %s wind [km]',v_usr);
        input_title1=sprintf('Time series: Radius of %s wind [km] (%s)',v_usr,smooth_str);
    else
        input_ylabel=sprintf('radius of%5.1f m/s wind [km]',v_usr);
        input_title1=sprintf('Time series: Radius of%5.1f m/s wind [km] (%s)',v_usr,smooth_str);
    end
    ylabel(input_ylabel);
    input_title2=sprintf('z=%5.2f %s; from %5.1f to %5.1f',zmin_sub,zunits,min(r_usr),max(r_usr));

    title({input_title1,input_title2})

end

v_usr
r_usr_tmean

%}