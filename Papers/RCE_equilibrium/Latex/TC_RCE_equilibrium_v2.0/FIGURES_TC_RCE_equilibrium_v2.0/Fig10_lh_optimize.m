%Fig10_lh_optimize.m

%Purpose: to estimate the "correct" value of l_h that gives Vm/Vp = 1/sqrt(2)
%from ER11 Eq. 41

%Created: 24 Sep 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%
%%TESTING ONLY!!! %%%%%
lh_all = 1500*[1/8 1/4 1/2 1 2 4 8];
mpi = 92;
%%%%%%%%%%%%%%%%%%%%

subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'lh'};
T_mean = 2; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',2)
h=figure(1)

%resize the figure to be pdf print-ready
set(h,'Units','centimeters');
hpos = [0 0 20 20];
set(h,'Position',hpos);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',hpos);
set(h,'PaperSize',hpos(3:4));

ax1=axes('position',[0.1    0.1    0.85    0.85]);

%%Determine output subdirectory pathname for given sim_set
if(equil_dynamic == 1)
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
    else
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
    end
else
    if(wrad_const == 1)
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
    else
        subdir_out = sprintf('../CM1_postproc_data/simsets_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
    end
end

if(wrad_const == 1)
    wrad_str = 'ctrl';
else
    wrad_str = 'rce';
end

sim_set=sim_sets{1};
load(sprintf('%s/%s.mat',subdir_out,sim_set))

for m=1:length(sim_sets)
%{
    sim_set = sim_sets{m};  %string
    load(sprintf('%s/%s.mat',subdir_out,sim_set));
    pl_clrs={'b' 'r' 'g' 'c' 'k' 'm' 'y' 'b--' 'r--' 'g--' 'c--' 'k--' 'm--' 'y--'};
    numruns = length(subdirs_set);
    
    clear fcor_all lh_all mpi_all
    for ss=1:numruns

        subdir = subdirs_set{ss};
        %%Determine output subdirectory pathname for given sim_set
        if(equil_dynamic == 1)
            if(wrad_const == 1)
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic_wradconst',T_mean,dt_final_dynamic);
            else
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_dt%i_dynamic',T_mean,dt_final_dynamic);
            end
        else
            if(wrad_const == 1)
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i_wradconst',T_mean,tf-dt_final,tf);
            else
                subdir_out2 = sprintf('../CM1_postproc_data/simdata_Tmean%i_%i_%i',T_mean,tf-dt_final,tf);
            end
        end

        if(wrad_const == 1)
            wrad_str = 'ctrl';
        else
            wrad_str = 'rce';
        end
        %%Load data for given simulation
        load(sprintf('%s/ax%s.mat',subdir_out2,subdir));

        %% Other
        fcor_all(ss) = fcor;
        mpi_all(ss) = mpi;
        lh_all(ss) = lh;

    end
%}    

xvals_pl(1:numruns,m) = lh_all;        %values defined by user at top
x_mpi = 0:10:12000;
y_mpi=mpi/sqrt(2)*ones(length(x_mpi),1);

axes(ax1)
clr1 = [0 0 1];

%% Best fit curve to data
%http://www.mathworks.com/matlabcentral/fileexchange/29545-power-law-exponential-and-logarithmic-fit
% logfit(X,Y), will search through all the possible axis scalings and 
%              finish with the one that incurs the least error (with error 
%              measured as least squares on the linear-linear data.)

%A power law relationship 
x_fit = min(xvals_pl):10:max(xvals_pl);

[slope, intercept] = logfit(xvals_pl,Vmax_equil_g,'loglog'); 
           y_fit = (10^intercept)*x_fit.^(slope);
           
%An exponential relationship 
%[slope, intercept] = logfit(xvals_pl,Vmax_equil_g,'logy'); 
%           y_fit = (10^intercept)*(10^slope).^x_fit;

%Find optimal l_h at crossing point
lh_optimal = x_fit(find(y_fit<y_mpi(1),1))

plot(xvals_pl,Vmax_equil_g,'x','MarkerEdgeColor',clr1,'MarkerSize',20)
hold on
plot(x_mpi,y_mpi,'r--')
plot(x_fit,y_fit,'b')
plot(lh_optimal,y_mpi(1),'g.','MarkerSize',20)
y_optim = min(get(ax1,'YLim')):(max(get(ax1,'YLim'))-min(get(ax1,'YLim')))/50:max(get(ax1,'YLim'));
plot(lh_optimal*ones(length(y_optim)),y_optim,'--','Color',[.8 .8 .8],'LineWidth',1)
xlabel('$l_h$','Interpreter','Latex','FontSize',18)
ylabel('$V_m/V_p$','Interpreter','Latex','FontSize',18)
legend({'$V_m$','$\frac{V_p}{\sqrt{2}}$'},'Interpreter','Latex','FontSize',18,'FontWeight','bold')
%title('Testing ER11 Eq 41 scaling for V_m/V_p : l_h-adjusted')
%axis([-pl_edge pl_edge -pl_edge pl_edge])

end
%}

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig10_lh_optimize.pdf