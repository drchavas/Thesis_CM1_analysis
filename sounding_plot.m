%sounding_plot.m

%Plot a sounding from an input datafile

%Created: 12 Mar 2012, Dan Chavas

clc
clear
figure(1)
clf(1)

%% USER INPUT %%%%%%%%%%%%%%%%%%
dz = .625;  %desired vertical resolution
z_top = 25000;  %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
Rd=287;  %[J/kg/K]
Rv=461.5;   %[J/K/kg]
Cpd=1005.7; %[J/kg/K]; spec heat of dry air
epsilon=Rd/Rv;
g=9.81; %[m/s2]

snd_files = {
%'input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3'
}

pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'm.' 'y.'};
pl_clrs2={'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'm-.' 'y-.'};

numruns=length(snd_files);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)


%%Append ax or 3d to subdir names for plot
for i=1:numruns

    snd_file = snd_files{i};
    
    nz_sub = (z_top/1000)/dz;


    %%load initial sounding variables
    dir_full = '/Users/drchavas/Documents/Research/Thesis/CM1/v15/analysis/input_soundings/';
    if(i==5)
        dir_full = '/Volumes/CHAVAS_CM1_FINAL/CM1_output/axisym/CTRL_icRCE/'
    end
    
    [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);
    zz00(end)


    %% Plot sounding
    figure(1)
    subplot(1,2,1)
    hold on
    plot(1000*[qv00],[zz00/1000],pl_clrs{i},'LineWidth',1)

        set(gca,'fontweight','bold','fontsize',11)
        axis([0 20 0 z_top/1000])
        input_xlabel=sprintf('water vapor mixing ratio [g/kg]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('Height AGL [km]');
        ylabel(input_ylabel);
        input_title=strrep(sprintf('water vapor mixing ratio [g/kg]'),'_','\_');
        title(input_title)

    figure(1)
    subplot(1,2,2)
    hold on
    plot([th00],[zz00/1000],pl_clrs{i},'LineWidth',1)
    plot([T00],[zz00/1000],pl_clrs2{i},'LineWidth',1)

        set(gca,'fontweight','bold','fontsize',11)
        axis([50 600 0 z_top/1000])
        input_xlabel=sprintf('T and theta [K]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('Height AGL [km]');
        ylabel(input_ylabel);
        input_title=strrep(sprintf('T and theta [K]'),'_','\_');
        title(input_title)

end

        input_legend=strrep(snd_files(1:numruns),'_','\_');
        legend(input_legend)
        
[zz00'/1000 T00']
nz_sub
Ttpp = str2num(snd_file(strfind(snd_file,'Tthresh')+7:strfind(snd_file,'Tthresh')+9))
if(isempty(Ttpp))
    Ttpp = 200;
end
ztpp = zz00(find(T00<Ttpp,1)-1)/1000+dz*(T00(find(T00<Ttpp,1)-1)-Ttpp)/(T00(find(T00<Ttpp,1)-1)-T00(find(T00<Ttpp,1)))
