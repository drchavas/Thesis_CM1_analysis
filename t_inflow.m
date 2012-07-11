%t_inflow.m

%Purpose: to calculate the timescale for inflow from some radius to Rmax

%Created: 9 Aug 2011, Dan Chavas

%NOTE: run this program after having run nc_plot_steady for the variable 'uinterp'
%at the desired model level (i.e. the lowest one) for the appropriate averaging
%time (e.g. the final 30 days of simulation)

clc

%%variables:
%rprof_tmean = u data
%dx = resolution
%xvals = radii; dx/2:dx:end-dx/2


%%Estimate Rmax as region of peak convergence (-1/r)*(d(ru)/dr)
xvals_half = xvals(1:end-1)+2;  %dx:dx:end-dx
conv=(-1./xvals_half).*(diff(xvals.*rprof_tmean')); 
i_convmax=find(conv==max(conv));
Rmax = xvals_half(i_convmax)    %[km]

%%Find radius where inflow begins (i.e. u<0 permanently)
inflow=rprof_tmean(i_convmax+1:end);
i_outer=i_convmax+find(inflow>0,1)-1;
Rout = xvals_half(i_outer)    %[km]
if(isempty(i_outer))
    Rwall = .9349*(xvals(end))-mod(.9349*xvals(end),dx)+dx/2
    i_outer = find(xvals==Rwall)
end

%%Integrate dr*(1/v) inwards (dr=constant) from outer pt to Rmax
T_inflow = -1000*dx*sum(1./rprof_tmean(i_convmax+1:i_outer));  %[s]
T_inflow = T_inflow/60/60/24    %[day]