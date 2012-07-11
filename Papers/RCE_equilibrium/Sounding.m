%Sounding.m

%Created: 28 Mar 2012, Dan Chavas

%This file is the same as sounding_plot.m, except that produces a
%single-figure sounding plot for the Control case

clear
clc
figure(1)
clf(1)

cd ../..

%set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)


%% USER INPUT %%%%%%%%%%%%%%%%%%
dz = .625;  %desired vertical resolution
z_top = 20000;  %[m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
Rd=287;  %[J/kg/K]
Rv=461.5;   %[J/K/kg]
Cpd=1005.7; %[J/kg/K]; spec heat of dry air
epsilon=Rd/Rv;
g=9.81; %[m/s2]

snd_files = {    
'input_sounding_3dControl_RCE'
%'input_sounding_3dRCE_nx48_SST300.00K_Tthresh200K_usfc3'
%'sounding_test_dan'
%'input_sounding_3dRCE_nx48_SST315.00K_Tthresh200K_usfc3'

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
    dir_full = '/Users/drchavas/Documents/Research/Thesis/CM1/v15/analysis';
    [zz00 pp00 th00 qv00 u00 v00 T00 Tv00 thv00 rho00 qvs00 rh00 pi00 p_sfc th_sfc qv_sfc] = snd_extract(dir_full,snd_file,dz,nz_sub);
    zz00(end)

    
    
    h=figure(1)
    set(h,'Position',[360 278 375 350])
    grid on
%{
    %subplot(3,1,1)
    axes(ax1)
%    axis([-2 2 -2 2])
    ylabel('Potential Temperature [K]')
    xlabel({sprintf('log_2(V_p/V^*_p)')})
    input_title1=sprintf('$V_m$');
    text1=text(-1.7,1.7,input_title1,'FontSize',17);
    set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
    %set(text1,'EdgeColor','k')
    set(ax1,'YTick',[-2 -1 0 1 2],'XTick',[-2 -1 0 1 2])
    grid on
    box on

    h=legend({'T_{sst}' 'T_{tpp}' 'Q_{cool}' 'u_s'},'Orientation','horizontal','Position',[0.15    0.15    0.70    0.02],'EdgeColor','white')
    grid on
%}

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
title('Control RCE sounding')
xlabel('potential temperature [K]')

    %{
    %% Plot sounding
    figure(1)
    subplot(1,2,1)
    hold on
    plot(1000*[qv00],[zz00/1000],pl_clrs{i},'LineWidth',2)

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
    plot([th00],[zz00/1000],pl_clrs{i},'LineWidth',2)
    plot([T00],[zz00/1000],pl_clrs2{i},'LineWidth',2)

        set(gca,'fontweight','bold','fontsize',11)
        axis([50 600 0 z_top/1000])
        input_xlabel=sprintf('T and theta [K]');
        xlabel(input_xlabel);
        input_ylabel=sprintf('Height AGL [km]');
        ylabel(input_ylabel);
        input_title=strrep(sprintf('T and theta [K]'),'_','\_');
        title(input_title)
    %}

end

%        input_legend=strrep(snd_files(1:numruns),'_','\_');
%        legend(input_legend)
        
[zz00'/1000 T00']
