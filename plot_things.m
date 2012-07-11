%plot_things.m

%clf(3)

%%vars: rmaxs, r_usr_means, Vmaxs

scale = [-3 -2 -1 0 1 2 3 4];
rmax=rmaxs;
rmax_CTRL = rmax(find(scale==0));
rmax=log(rmax/rmax_CTRL);
%fig_title = 'Sensitivity of rmax to self-sim';

r12 = r_usr_means;    %from output of nc_plot_steady.m
r12_CTRL = r12(find(scale==0));
r12=log(r12/r12_CTRL);


%r12 = [18 59.5 82 197.7 150.8 160.2 328.9]
%wts = [2^3 2^2 2^1 2^0 2^(-1) 2^(-2) 2^(-3)];

figure(3)
set(gca,'fontweight','bold','fontsize',12,'fontname','Arial')
hold off
%subplot(2,1,1)
plot(scale,scale)
hold on
plot(scale,r12,'b.')

plot(scale,rmax,'r.')


xlabel('log(r_{in} / r_{in, CTRL})')
ylabel('log(r_{out} / r_{out, CTRL})')

title({'Scaling of steady-state radii with halving/doubling'})

