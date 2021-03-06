%Dimensional_scaling.m

%Created: 11 Jan 2012, Dan Chavas

%This file is the same as TC_stats_plotset.m, except that it only produces
%the plot for all the relevant dimensional scalings (Figure 1).

clear
clc
figure(1)
clf(1)

cd ../..

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

%clf(1)

%'dx' 'dz' 'domain' 'lh' 'lv' 'fcor' 'qro' 'ro' 'Qrad' 'Tsst' 'Ttpp' 'usfc' 'mpi'
%sim_sets = {'lh' 'lv' 'fcor' 'qro' 'ro' 'mpi' 'Cd'}
sim_sets = {'lh' 'lv' 'fcor' 'qro' 'ro' 'Ttpp'}

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
    
    i_ctrl = find(strcmp(subdirs_set,'CTRLv0qrhSATqdz5000_nx3072')==1,1);
    if(isempty(i_ctrl))
        i_ctrl = find(multipliers==0,1);
    end
    mpi_ctrl = mpi_all(i_ctrl);

    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
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
    ax1=axes('position',[0.15    0.73   0.70    0.23]);
    end
    axes(ax1)
    data_temp = Vmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));

    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on

        %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
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
    ax2=axes('position',[0.15    0.43   0.70    0.23]);
    end
    axes(ax2)
    data_temp = rmax_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on

    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
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
        ax3=axes('position',[0.15    0.13   0.70    0.23]);
    end
    axes(ax3)
    data_temp = r0Lil_equil_g(i_sort);
    data_pl = log2(data_temp./data_temp(i_ctrl));
    dat_max = max(dat_max,max(data_pl));
    dat_min = min(dat_min,min(data_pl));
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on
   
    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
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
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
        plot(xvals_pl,data_pl,'x','Color',pl_clrs{m})
    else
        plot(xvals_pl,data_pl,'x-','Color',pl_clrs{m})
    end
    hold on
    
    %% Add best fit line for mpi (should look identical to MPI_collapse plots)
    if(strcmp(sim_set,'Ttpp'))    %sort mpi files in ascending order
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
%set(h,'Position',[100 278 900 700])
set(h,'Position',[0 278 300 700])

%subplot(3,1,1)
axes(ax1)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$V_m$');
text1=text(-2.7,2.7,input_title1,'FontSize',17);
set(text1,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
%set(text1,'EdgeColor','k')
set(ax1,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on

%subplot(3,1,2)
axes(ax2)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)')
%xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_m$');
%title(input_title1)
text2=text(-2.7,2.7,input_title1,'FontSize',17);
set(text2,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
set(ax2,'YTick',[-3 -2 -1 0 1 2 3],'XTick',[-3 -2 -1 0 1 2 3])
grid on
box on

%subplot(3,1,3)
axes(ax3)
axis([-3 3 -3 3])
ylabel('log_2(Y/Y*)')
xlabel({sprintf('log_2(X/X*)')})
input_title1=sprintf('$r_0$');
%title(input_title1)
text3=text(-2.7,2.7,input_title1,'FontSize',17);3
set(text3,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','Latex');
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

h=legend({'l_h' 'l_v' 'f' 'r_{0q}' 'r_{0u}' 'V_p'},'Orientation','horizontal','Position',[0.15    0.05    0.70    0.02],'EdgeColor','white')
grid on

cd Papers/RCE_equilibrium/
