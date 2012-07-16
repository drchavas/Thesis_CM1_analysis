%r0ER11_calc.m

%Created: 16 Jul 2012, Dan Chavas

%Purpose: calculate the theoretical r0 from Emanuel and Rotunno (2011),
%r0ER11, by fitting the ER11 model to the CM1 r(Vg=.1*Vp).


%%Use the ER11 estimate for Vm based on Vp, Ck/Cd; because here we are
%%interested in the far outer region of the storm where the horizontal
%%mixing length does not affect the equilibrium storm much, we don't want
%%to use the actual CM1 Vm to estimate r0ER11 since Vm is highly sensitive
%%to lh!

%%TESTING
%{
clear
clc
mpi = 91;    %[ms-1]
fcor = 5e-5;    %[s-1]
Cd_in = .0015;
v_usr_fracVp = .1;  %wind speed as fraction of Vp; beyond this radius, radiative subsidence radial wind profile should apply.
r_rad = 1000;    %r(Vg=.1*Vp)
%}

function [r0ER11] = r0ER11_calc(r_rad,v_usr_fracVp,mpi,fcor,Cd_in);

%%%%%%%%%%%%%%%%

CkCd_ctrl = .0015;
CkCd = CkCd_ctrl/Cd_in;
rr = 1:10000;    %[km]
rmax_guess_0 = 50;  %[km]
V_usr = v_usr_fracVp*mpi;   %[ms-1]

%%ER11 Eq. 41
Vmax_ER11 = mpi*(.5*CkCd)^(.5*(CkCd/(2-CkCd)));

rmax_guess = rmax_guess_0;
error = 10000;
niter = 0;

while(abs(error)>1 && niter < 10000)

    %%ER11 Eq. 36
    V_ER11 = Vmax_ER11*(rmax_guess./rr).*((2*(rr./rmax_guess).^2)./(2-CkCd+CkCd*(rr./rmax_guess).^2)).^(1/(2-CkCd)) - .5*fcor*1000*rr;
    i_r_ER11 = find(V_ER11>V_usr);
    r_ER11 = rr(i_r_ER11(end));
    error = r_ER11 - r_rad;
    
    if(error>0)
        rmax_guess = rmax_guess - .1;
    else
        rmax_guess = rmax_guess + .1;
    end
    niter = niter + 1;
end

if(niter > 9999)
    r0ER11 = NaN;
else
    i_r0ER11 = find(V_ER11>0);
    r0ER11 = rr(i_r0ER11(end));
end
%niter
%r0ER11
%{
figure(10)
clf(10)
plot(rr(V_ER11>0),V_ER11(V_ER11>0))
hold on
plot(r_rad,V_usr,'rx')
%}

end