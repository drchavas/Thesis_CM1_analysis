%Fig0half_Sounding_ctrl.m

%Created: 24 Sep 2012, Dan Chavas


clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis


%% USER INPUT %%%%%%%%%%%%%%%%%%
dz = .625;  %desired vertical resolution
z_top = 20000;  %[m]

snd_files = {    
'input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3_drag'
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants (values taken from CM1 model)
c_CM1 = constants_CM1(); %c_CM1: [g rd cp cv p00 xlv]

g=c_CM1(1); %[m/s2]
Rd=c_CM1(2);  %[J/kg/K]
Cpd=c_CM1(3); %[J/kg/K]; spec heat of dry air
Rv=c_CM1(4);   %[J/K/kg]
p0 = c_CM1(5); %[Pa]
Lv=c_CM1(6);   %[J/kg]

eps=Rd/Rv;

pl_clrs={'b' 'r' 'g' 'c' 'k' 'y' 'm' 'm.' 'y.'};
pl_clrs2={'b--' 'r--' 'g--' 'c--' 'k--' 'y--' 'm--' 'm-.' 'y-.'};

numruns=length(snd_files);  %total number of runs you want to plot (taken from below)

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
    xlabel('water vapor mixing ratio [g/kg]','FontSize',18)
    ylabel('altitude [km]','FontSize',18)

    ax2 = axes('Position',get(ax1,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','r','YColor','k');
    hl2 = line([th00],[zz00/1000],'Color','r','Parent',ax2);
    hl3 = line([T00],[zz00/1000],'Color','r','LineStyle','--','Parent',ax2);
    %title('Control RCE sounding')
    xlabel('temperature (dash) / potential temperature (solid) [K]','FontSize',18)

    %%Remove the top altitude label (overlaps with x-axis)
    temp = get(ax1,'YTick');
    set(ax1,'YTick',temp(1:end-1))

    temp = get(ax2,'YTick');
    set(ax2,'YTick',temp(1:end-1))
    
    %%Define x-axis labels so that grid lines overlap
    num_grid = 9;
    th_min = 170;
    th_max = 440;
    qv_min = 0;
    qv_max = 18;
    
    ticks_th = th_min:(th_max-th_min)/num_grid:th_max;
    ticks_qv = qv_min:(qv_max-qv_min)/num_grid:qv_max;
    set(ax1,'XTick',ticks_qv)
    set(ax1,'XLim',[qv_min qv_max])
    set(ax2,'XTick',ticks_th)
    set(ax2,'XLim',[th_min th_max])


    
end
        
cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig0half_Sounding_ctrl.pdf

