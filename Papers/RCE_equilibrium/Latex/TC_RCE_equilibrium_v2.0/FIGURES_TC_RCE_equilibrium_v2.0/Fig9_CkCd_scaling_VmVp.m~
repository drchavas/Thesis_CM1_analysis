%Fig9_CkCd_scaling_VmVp.m

%Compare the scaling of V_m/V_p with C_k/C_d (ER11 Eqn 41)

%Created: 26 Mar 2012, Dan Chavas
%Updated: 19 Sep 2012, Dan Chavas

clear
clc
figure(1)
clf(1)

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis

%% USER INPUT %%%%%%%%%%%%%%%%%%
subdir_pre='CTRL_icRCE/';    %general subdir that includes multiple runs within
ext_hd = 1; %0=local hard drive; 1=external hard drive

sim_sets = {'Cd'};
T_mean = 500; %[day]
equil_dynamic = 1;  %1 = use dynamic equilibrium
    %%IF 0:
    dt_final = 50;
    tf = 150;
    %%IF 1:
    dt_final_dynamic = 30;  %[days]; new length of period over which equilibrium is calculated
wrad_const = 0; %1 = use CTRL value for wrad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)
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
load(sprintf('%s/%s',subdir_out,sim_set))

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

%clearvars -except Vmax_equil_g mpi_all

xvals_pl(1:numruns,m) = multipliers;        %values defined by user at top

%%NEED CORRECT MPI VALUES!
mpi_all = [261.1 200.5 133.1 mpi_all(4) 63.6 44.7 31.6];
Ck = .0015;
Cd_all = .0015*[0.125 .25 .3501 1 2 4 8];
CkCd = Ck./Cd_all;

%make Vmax_equil_g/Vp match the predicted value by tuning lh -- but don't
%need to rerun the model, just use the scaling law
%lh_orig = 1500; %[m]
%lh_new = lh_orig*(mpi_all(4)/(sqrt(2)*Vmax_equil_g(4)))^1.15
Vm_adj = mpi_all(4)/(sqrt(2)*Vmax_equil_g(4))   %equals fractional adjustment due to scaling of Vm with lh

pl_edge = max([abs(floor(min(xvals_pl(1:numruns,m)))) abs(ceil(max(xvals_pl(1:numruns,m)))) 3]);

axes(ax1)
x = log2((CkCd)/(1));

CkCd_smooth = min(CkCd):(max(CkCd)-min(CkCd))/100:max(CkCd);
x_smooth = log2((CkCd_smooth)/(1));
ER11_eq41 = ((CkCd_smooth)./2).^((CkCd_smooth)./(2*(2-(CkCd_smooth))));

clr1 = [0 0 1];

plot(x_smooth,ER11_eq41,'k--',x,Vm_adj.*Vmax_equil_g./mpi_all,'x','MarkerEdgeColor',clr1,'MarkerSize',20,'LineWidth',2)
xlabel('$\\log_2(C_k/C_d)$','Interpreter','Latex','FontSize',18)
ylabel(sprintf('$V_m/V_p$'),'Interpreter','Latex','FontSize',18)
legend('ER11','CM1','FontSize',18,'FontWeight','bold')
%title('Testing ER11 Eq 41 scaling for V_m/V_p : l_h-adjusted')
%axis([-pl_edge pl_edge -pl_edge pl_edge])

end
%}

cd /Users/drchavas/Documents/Research/Thesis/CM1/v15/Thesis_CM1_analysis/Papers/RCE_equilibrium/Latex/TC_RCE_equilibrium_v2.0/FIGURES_TC_RCE_equilibrium_v2.0

print -dpdf -r300 Fig9_CkCd_scaling_VmVp.pdf