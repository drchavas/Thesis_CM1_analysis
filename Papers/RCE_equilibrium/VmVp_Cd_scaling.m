%VmVp_Cd_scaling.m

%Compare the scaling of V_m/V_p with C_k/C_d (ER11 Eqn 41)

%Created: 26 Mar 2012, Dan Chavas

clear
clc

cd ../..

load simsets/Cd.mat
clearvars -except Vmax_equil_g mpi_all

%make Vmax_equil_g/Vp match the predicted value by tuning lh -- but don't
%need to rerun the model, just use the scaling law

%lh_orig = 1500; %[m]
%lh_new = lh_orig*(mpi_all(4)/(sqrt(2)*Vmax_equil_g(4)))^1.15
Vm_adj = mpi_all(4)/(sqrt(2)*Vmax_equil_g(4))   %equals fractional adjustment due to scaling of Vm with lh


Ck = 1.5;
Cd = 1.5*[0.125 .25 .501 1 2 4 8];
x = log2((Ck./Cd)/(1));
ER11_eq41 = ((Ck./Cd)./2).^((Ck./Cd)./(2*(2-(Ck./Cd))));


set(0,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultlinelinewidth',1)

h=figure(1)
clf(1)

set(h,'Position',[160 578 575 400])

ax1=axes('position',[0.15    0.15    0.75    0.75]);
axes(ax1)


plot(x,ER11_eq41,'bx-',x,Vm_adj.*Vmax_equil_g./mpi_all,'gx-')
xlabel('log_2[(Ck/Cd)/1]')
ylabel('V_m/V_p')
legend('ER11','model')
title('Testing ER11 Eq 41 scaling for V_m/V_p : l_h-adjusted')


%}