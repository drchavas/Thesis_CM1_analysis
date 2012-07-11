%VmaxRmax.m

%Purpose to interpolate to the true Vmax and rmax in the radial profile

%Created: 12 May 2011, Dan Chavas

function [Vmax rmax] = VmaxRmax(xvals,rprof_in)
%%load test data
%load vmaxrmaxdat
%%vars: xvals rprof_tmean
data = rprof_in;
clear rprof_in

%%find Vmax, rmax from data
i_max = find(data==max(data));  %index of max wind speed
%i_max2 = find(data==max(data(data<max(data))))  %index of second highest wind speed

if(i_max==1)    %occurs in super-low-res runs sometimes
    Vmax = data(i_max);
    rmax = xvals(i_max);
else
    %%figure out if rmax in data is to inside or outside of true rmax
    v3temp = data(i_max);   %Vmax from data
    v2temp = data(i_max-1);
    if(i_max>2)
        v1temp = data(i_max-2);
    else
        v1temp = v2temp-(v3temp-v2temp);
    end

    %if rmax_data is inside of true rmax, then slope r<rmax_data should be approx constant
    if((v3temp-v2temp)/(v2temp-v1temp)<.98) %rmax_data outside eyewall
        %%save relevant variables
        v3 = data(i_max);   %Vmax from data
        r3 = xvals(i_max);  %Rmax from data
        v2 = data(i_max-1);
        r2 = xvals(i_max-1);
        v4 = data(i_max+1);
        r4 = xvals(i_max+1);
        if(i_max==2)
            v1 = data(1)-(data(2)-data(1));
            r1 = xvals(1)-(xvals(2)-xvals(1));
        else
            v1 = data(i_max-2);
            r1 = xvals(i_max-2);
        end    
    else    %rmax_data inside eyewall
        %%save relevant variables
        v3 = data(i_max+1);   
        r3 = xvals(i_max+1);
        v2 = data(i_max);   %Vmax from data
        r2 = xvals(i_max);  %Rmax from data
        v4 = data(i_max+2);
        r4 = xvals(i_max+2);
        v1 = data(i_max-1);
        r1 = xvals(i_max-1);
    end

    %%calculate linearly interpolated velocities, v2p and v1p, at r2 and r1 respectively
    %%note: assumes constant radial grid resolution!
    v2p = v3+(v3-v4);
    v3p = v2+(v2-v1);

    %%calculate coefficients a,b,c,d of y=ax+b and y=cx+d lines for each linear interpolation
    %Inside eyewall: line from (r2,v2) --> (r3,v3p) y=x+b
    a=(v2-v3p)/(r2-r3);
    b=v2-r2*(v2-v3p)/(r2-r3);

    %Outside eyewall: line from (r3,v3) --> (r2,v2p) y=cx+d
    c=(v3-v2p)/(r3-r2);
    d=v3-r3*(v3-v2p)/(r3-r2);

    %%Calculate interpolated Vmax and rmax from above coefficients
    rmax = (d-b)/(a-c);
    Vmax = c*rmax+d;

    %%Check: if answer is smaller than max value in actual data, then go with actual data
    if(Vmax < data(i_max))
        Vmax = data(i_max);
        rmax = xvals(i_max);
    end
end

%{
%%Plot things to check
plot(xvals,data','r','LineWidth',2)
axis([0 400 0 2*max(rprof_tmean)])
set(gca,'fontweight','bold','fontsize',11)
xunits = 'km';
input_xlabel=sprintf('x-distance [%s]',xunits);
xlabel(input_xlabel);
v_def = 'wind speed';
v_units = 'm/s';
input_ylabel=sprintf('%s [%s]',v_def,v_units);
ylabel(input_ylabel);
input_title1=sprintf('%s [%s]',v_def,v_units);
title(input_title1)

hold on
plot(rmax,Vmax,'*')
plot(r3,v3p,'o')
plot(r2,v2p,'s')
%}

end
