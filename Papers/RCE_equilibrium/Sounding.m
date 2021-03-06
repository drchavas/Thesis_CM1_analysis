%Sounding.m

%Created: 28 Mar 2012, Dan Chavas

%This file is the same as sounding_plot.m, except that produces a
%single-figure sounding plot for the Control case

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)


%% USER INPUT %%%%%%%%%%%%%%%%%%
dz = .625;  %desired vertical resolution
z_top = 20000;  %[m]

snd_files = {    
'input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag'
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
Rd=287;  %[J/kg/K]
Rv=461.5;   %[J/K/kg]
Cpd=1005.7; %[J/kg/K]; spec heat of dry air
epsilon=Rd/Rv;
g=9.81; %[m/s2]

pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'm.' 'y.'};
pl_clrs2={'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'm-.' 'y-.'};

numruns=length(snd_files);  %total number of runs you want to plot (taken from below)
%numruns=6;  %total number of runs you want to plot (taken from below)



set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',2)
h=figure(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 20 20];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.1    0.1    0.8    0.8]);


%%Append ax or 3d to subdir names for plot
for i=1:numruns

    snd_file = snd_files{i};
    
    nz_sub = (z_top/1000)/dz;


    %%load initial sounding variables
    dir_full = '/Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/input_soundings/';
    [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);
    zz00(end)

    
    
    grid on

    hl1 = line(1000*[qv00],[zz00/1000],'Color','b');
    ax1 = gca;
    set(ax1,'XColor','b','YColor','k')
    xlabel('water vapor mixing ratio [g/kg]')
    ylabel('altitude [km]')

    ax2 = axes('Position',get(ax1,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','r','YColor','k');
    hl2 = line([th00],[zz00/1000],'Color','r','Parent',ax2);
    hl3 = line([T00],[zz00/1000],'Color','r','LineStyle','--','Parent',ax2);
    %title('Control RCE sounding')
    xlabel('Temperature (dash) / Potential temperature (solid) [K]')

    %%Remove the top altitude label (overlaps with x-axis)
    temp = get(ax1,'YTick')
    set(ax1,'YTick',temp(1:end-1))

    temp = get(ax2,'YTick')
    set(ax2,'YTick',temp(1:end-1))

end

%        input_legend=strrep(snd_files(1:numruns),'_','\_');
%        legend(input_legend)
        
cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig0.5_sounding.pdf

