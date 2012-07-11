%plot_things.m

%clf(3)

clc
clear

load r_v
%%vars: v_usr_v, r_usr_means_v, rmaxs
load r_satdef
%%vars: v_usr_satdef, r_usr_means_satdef


scale = [-4 -3 -2 -1 0 1 2 3 4];
rmax=rmaxs;
rmax_CTRL = rmax(find(scale==0));
rmax=log2(rmax/rmax_CTRL);
%fig_title = 'Sensitivity of rmax to self-sim';

r12 = r_usr_means_v;    %from output of nc_plot_steady.m
r12_CTRL = r12(find(scale==0));
r12=log2(r12/r12_CTRL);

rsd = r_usr_means_satdef;    %from output of nc_plot_steady.m
rsd_CTRL = rsd(find(scale==0));
rsd=log2(rsd/rsd_CTRL);


%r12 = [18 59.5 82 197.7 150.8 160.2 328.9]
%wts = [2^3 2^2 2^1 2^0 2^(-1) 2^(-2) 2^(-3)];

figure(3)
set(gca,'fontweight','bold','fontsize',12,'fontname','Arial')
hold off
%subplot(2,1,1)
plot(scale,r12,'b.','MarkerSize',20)
hold on
plot(scale,rmax,'r.','MarkerSize',20)
plot(scale,rsd,'g.','MarkerSize',20)

set(gca,'fontweight','bold','fontsize',12,'fontname','Arial')

plot(scale,scale,'k','LineWidth',2)


xlabel('log_2(R* / R*_{CTRL})')
ylabel('log_2(r / r_{CTRL})')

%title({'Scaling of steady-state radii with halving/doubling of horiz scales (R*)'})

legend('r_{12}','r_{max}','r_{satdef}','Location','SouthEast')
