%Fig6_Dimensional_scaling.m

%Created: 5 Sep 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%clf(1)

%'dx' 'dz' 'domain' 'lh' 'lv' 'fcor' 'qro' 'ro' 'Qrad' 'Tsst' 'Ttpp' 'usfc' 'mpi'
sim_sets = {'lh_drag' 'lv_drag' 'fcor_drag' 'qro_drag' 'ro_drag' 'Ttpp_drag'}

T_mean = 2; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat_max=0;
dat_min=0;

%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out2 = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

if(wrad_const == 1)
    wrad_str = 'ctrl';
else
    wrad_str = 'rce';
end


for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out2,sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    
    if(strcmp(sim_set,'ro_drag'))
        i_ctrl = find(strcmp(subdirs_set,'CTRLv12.5qrh0qdz5000_nx3072_drag')==1,1);
    else
        i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    end
    if(isempty(i_ctrl))
        i_ctrl = find(multipliers==0,1);
    end
    mpi_ctrl = mpi_all(i_ctrl);

    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        [junk i_sort] = sort(mpi_all);
        clear junk
        multipliers = log2(mpi_all(i_sort)./mpi_ctrl);
        i_ctrl = find(multipliers==0,1);
    else
        i_sort = 1:length(mpi_all);
    end
    
    xvals_pl = multipliers;    %values defined by user at top except for 'mpi'
    
    figure(1)
%    subplot(3,1,1)
    if(m==1)
    %ax1=axes('position',[0.1    0.6    0.35    0.35]);
    ax1=axes('position',[0.2    0.75   0.70    0.23]);
    end
    axes(ax1)
    data_temp = Vmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));

    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on

        %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        %options = fitoptions('Method','Smooth','SmoothingParam',0.3)
        %f = fit(xvals_pl_all', data_pl_all', 'smooth',options)
        f = fit(xvals_pl', data_pl', 'poly1')
        % plot(f, xvals_pl_all, data_pl_all)
        %%Linear model: f(x) = p1*x + p2
         xmin_pl = min(xvals_pl);
         xmax_pl = max(xvals_pl);
         xdiff_pl = xmax_pl-xmin_pl;
         xfit = xmin_pl-xdiff_pl/3:xdiff_pl/30:xmax_pl+xdiff_pl/3
         yfit = f.p1.*xfit + f.p2;
         plot(xfit,yfit,'m-')
    end
    
%    subplot(3,1,2)
    if(m==1)
    %ax2=axes('position',[0.55    0.6    0.35    0.35]);
    ax2=axes('position',[0.2    0.45   0.70    0.23]);
    end
    axes(ax2)
    data_temp = rmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on

    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        %options = fitoptions('Method','Smooth','SmoothingParam',0.3)
        %f = fit(xvals_pl_all', data_pl_all', 'smooth',options)
        f = fit(xvals_pl', data_pl', 'poly1')
        % plot(f, xvals_pl_all, data_pl_all)
        %%Linear model: f(x) = p1*x + p2
         xmin_pl = min(xvals_pl);
         xmax_pl = max(xvals_pl);
         xdiff_pl = xmax_pl-xmin_pl;
         xfit = xmin_pl-xdiff_pl/3:xdiff_pl/30:xmax_pl+xdiff_pl/3
         yfit = f.p1.*xfit + f.p2;
         plot(xfit,yfit,'m-')
    end

%    subplot(3,1,3)
    if(m==1)
        %ax3=axes('position',[0.1    0.2    0.35    0.35]);
        ax3=axes('position',[0.2    0.15   0.70    0.23]);
    end
    axes(ax3)
    data_temp = r0Lil_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        hpl(m)=plot(xvals_pl,data_pl,'x','Color',pl_clrs{m});
    else
        hpl(m)=plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m});
    end
    hold on
   
    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        %options = fitoptions('Method','Smooth','SmoothingParam',0.3)
        %f = fit(xvals_pl_all', data_pl_all', 'smooth',options)
        f = fit(xvals_pl', data_pl', 'poly1')
        % plot(f, xvals_pl_all, data_pl_all)
        %%Linear model: f(x) = p1*x + p2
         xmin_pl = min(xvals_pl);
         xmax_pl = max(xvals_pl);
         xdiff_pl = xmax_pl-xmin_pl;
         xfit = xmin_pl-xdiff_pl/3:xdiff_pl/30:xmax_pl+xdiff_pl/3
         yfit = f.p1.*xfit + f.p2;
         plot(xfit,yfit,'m-')
    end
%{    
%    subplot(3,1,3)
    if(m==1)
        ax4=axes('position',[0.55    0.2    0.35    0.35]);
    end
    axes(ax4)
    data_temp = r0Lil_Lilctrl_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on
    
    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp_drag'))    %sort mpi files in ascending order
        %options = fitoptions('Method','Smooth','SmoothingParam',0.3)
        %f = fit(xvals_pl_all', data_pl_all', 'smooth',options)
        f = fit(xvals_pl', data_pl', 'poly1')
        % plot(f, xvals_pl_all, data_pl_all)
        %%Linear model: f(x) = p1*x + p2
         xmin_pl = min(xvals_pl);
         xmax_pl = max(xvals_pl);
         xdiff_pl = xmax_pl-xmin_pl;
         xfit = xmin_pl-xdiff_pl/3:xdiff_pl/30:xmax_pl+xdiff_pl/3
         yfit = f.p1.*xfit + f.p2;
         plot(xfit,yfit,'m-')
    end
%}   
    
end

h=figure(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 15 30];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

%subplot(3,1,1)
axes(ax1)
%set(ax1,'FontWeight','normal','FontSize',16)
axis([-3 3 -3 3])
ylab = ylabel('$log_2(V_m/V^*_m)$')
set(ylab,'Interpreter','Latex','FontSize',18);
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$V_m$');
text1=text(-2.8,2.8,input_title1,'FontSize',18);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','EdgeColor','k');
%set(text1,'EdgeColor','k')
set(ax1,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on

%subplot(3,1,2)
axes(ax2)
%set(ax2,'FontWeight','normal','FontSize',16)
axis([-3 3 -3 3])
ylab = ylabel('$log_2(r_m/r^*_m)$')
set(ylab,'Interpreter','Latex','FontSize',18);
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_m$');
%title(input_title1)
text2=text(-2.8,2.8,input_title1,'FontSize',18);
set(text2,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','EdgeColor','k');
set(ax2,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on

%subplot(3,1,3)
axes(ax3)
%set(ax3,'FontWeight','normal','FontSize',16)
axis([-3 3 -3 3])
ylab = ylabel('$log_2(r_0/r^*_0)$')
set(ylab,'Interpreter','Latex','FontSize',18);
xlab = xlabel({sprintf('$log_2(X/X^*)$')})
set(xlab,'Interpreter','Latex','FontSize',18);
input_title1=sprintf('$r_0$');
%title(input_title1)
text3=text(-2.8,2.8,input_title1,'FontSize',18);
set(text3,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex','EdgeColor','k');
set(ax3,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on

%{
%subplot(3,1,3)
axes(ax4)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_{0ctrl}$');
%title(input_title1)
text3=text(-2.7,2.7,input_title1,'FontSize',17);
set(text3,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
set(ax3,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on
%}

leg_array={'$l_h$' '$l_v$' '$f$' '$r_{0q}$' '$r_{0u}$' '$V_p$'};
LegendHandle=legend(leg_array,'Orientation','horizontal','Position',[.55 .05 .5 .05],'EdgeColor','white','FontSize',18)
set(LegendHandle,'Interpreter','Latex','Box','off');
%%note: matlab has a bug and the 'Position' input above will only change
%%the location, not the size of the box (or at least it won't make the box
%%smaller) DRC 6 Sep 2012

%% Break up the legend into multiple rows
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/157048
%{
% Block 1
% Axes handle 1 (this is the visible axes)
ah1 = gca;
% Legend at axes 1
legend(ah1,hpl(1:3),'l_h','l_v','f',1,'Orientation','horizontal','Position',[175    30    .2    0.02],'EdgeColor','black','FontSize',14)

% Block 2
% Axes handle 2 (unvisible, only for place the second legend)
ah2=axes('position',get(gca,'position'), 'visible','off');
% Legend at axes 2
legend(ah2,hpl(4:6),'r_{0q}','r_{0u}','V_p',2,'Orientation','horizontal','Position',[175    10    .2    0.02],'EdgeColor','black','FontSize',14)
%}

%% Change the length of the lines in the legend

%http://www.mathworks.com/matlabcentral/newsreader/view_thread/155063
%find the children
Children = get(LegendHandle, 'Children') ;
%Fot line case, the number of childre always n*3, n is the
%number of line legend. In another word (every 3 objects create one lenged). 

%recogonize the children(every 3 objects). Using get(Children(Index)) ;
%The first child is the marker.
%The second one is the line
%The last one is the text
compress = 1.2; %1 = no compression; >1 = horiz compression of legend
for i=1:length(leg_array)
    
    %%rescale line lengths
    i_chil = length(Children)-3*(i-1);
    XData = get(Children(i_chil-1), 'XData') ;
    if(i==1)    %determine length of original line
        XScale = XData(2) - XData(1) ;
        XScale_new = .3*XScale;   %length of new line
        dXScale = XScale - XScale_new;
    end
    XData(1) = XData(1) - compress*dXScale*(i-1) ;
    XData(2) = XData(1) + XScale_new;
    set(Children(i_chil-1), 'XData', XData) ;

    %%rescale marker position
    XData = get(Children(i_chil-2), 'XData') ;
    XData = XData - compress*dXScale*(i-1) - dXScale/2 ; % midpoint of shifted line
    set(Children(i_chil-2), 'XData', XData) ;

    %%rescale text location
    XData = get(Children(i_chil), 'Position') ;
    XData(1) = XData(1) - compress*dXScale*(i-1) - 2*dXScale/2 ; % shift text with marker
    set(Children(i_chil), 'Position', XData) ;

    %%shift texts with subscripts downwards
%{
    if(i_chil~=12)  %skip 'f'
        postemp = get(Children(i_chil),'Position');
        postemp(2) = .7*postemp(2);
        set(Children(i_chil),'Position',postemp);
    else
        postemp = get(Children(i_chil),'Position');
        postemp(2) = 1.05*postemp(2);
        set(Children(i_chil),'Position',postemp);        
    end    
%}    
end

%%add a line to Vp legend item
set(Children(2),'LineStyle','-')

%}

grid on

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig6_Dimensional_scaling.pdf
