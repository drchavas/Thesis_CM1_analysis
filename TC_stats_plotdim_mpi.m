%TC_stats_plotdim_mpi.m

%Created: 14 Feb 2012, Dan Chavas

%This file is the same as TC_stats_plotdim.m, except that it plots with color-coding the 4 input variables that modulate MPI
%and displays them all togther based upon their respective mpi values.

clear
clc
figure(1)
clf(1)

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%clf(1)

%%variables of interest (sim_set name): 'dx' 'dz' 'domain' 'lh' 'lv' 'H' 'Qrad' 'Vpot' 'cor' 'qro' 'ro' 'rodrmax'
sim_sets = {'Tsst' 'Ttpp' 'Qrad' 'usfc' 'Cd'}
%sim_sets = {'test'}

dat_max=0;
dat_min=0;

for m=1:length(sim_sets)
    
    sim_set = sim_sets{m};  %string
    load(sprintf('simsets/%s.mat',sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    pl_shapes={'x' '*' 's' 'd' '+'};
    
    i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    mpi_ctrl = mpi_all(i_ctrl);

    [junk i_sort] = sort(mpi_all);
    clear junk
    multipliers = log2(mpi_all(i_sort)/mpi_ctrl);
    
    xvals_pl = multipliers;    %values defined by user at top
    
    figure(1)
%    subplot(3,1,1)
    if(m==1)
    ax1=axes('position',[0.15    0.71    0.70    0.27]);
    end
    axes(ax1)
    data_temp = Vmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    scatter(xvals_pl,data_pl,'Marker',pl_shapes{m},'MarkerEdgeColor',pl_clrs{m})
    hold on

%    subplot(3,1,2)
    if(m==1)
    ax2=axes('position',[0.15    0.41    0.70    0.27]);
    end
    axes(ax2)
    data_temp = rmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    scatter(xvals_pl,data_pl,'Marker',pl_shapes{m},'MarkerEdgeColor',pl_clrs{m})
    hold on

%    subplot(3,1,3)
    if(m==1)
    ax3=axes('position',[0.15    0.11    0.70    0.27]);
    end
    axes(ax3)
    data_temp = r0Lil_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    scatter(xvals_pl,data_pl,'Marker',pl_shapes{m},'MarkerEdgeColor',pl_clrs{m})
    hold on
   
end

h=figure(1)
set(h,'Position',[360 278 375 700])

%subplot(3,1,1)
axes(ax1)
axis([-2 2 -2 2])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$V_m$');
text1=text(-1.7,1.7,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
%set(text1,'EdgeColor','k')
set(ax1,'YTick',[-2 -1 0 1 2],'XTick',[-2 -1 0 1 2])
grid on

%subplot(3,1,2)
axes(ax2)
axis([-2 2 -2 2])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_m$');
%title(input_title1)
text2=text(-1.7,1.7,input_title1,'FontSize',17);
set(text2,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
set(ax1,'YTick',[-2 -1 0 1 2],'XTick',[-2 -1 0 1 2])
grid on

%subplot(3,1,3)
axes(ax3)
axis([-2 2 -2 2])
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_0$');
%title(input_title1)
text3=text(-1.7,1.7,input_title1,'FontSize',17);
set(text3,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
set(ax1,'YTick',[-2 -1 0 1 2],'XTick',[-2 -1 0 1 2])
grid on
box on

h=legend({'T_{sst}' 'T_{tpp}' 'Q_{rad}' 'u_s' 'C_d'},'Orientation','horizontal','Position',[0.15    0.035    0.70    0.02],'EdgeColor','white')
grid on
