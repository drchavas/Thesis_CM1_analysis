%r0_calc.m

% Purpose: to calculate the simple asymptotic and full radial windspeed
% profile solutions for a given input (r,V) pair at a given latitude.  The
% full solution of eqn. 6 (Dean et al 2009) is solved numerically and
% converges from above; the simple asymptotic solution is given as eqn. 8
% (Dean et al 2009).

% Ref 1: Dean et. al, 2009. "On the size distribution of Atlantic tropical
% cyclones." Geo Res Let 36, L14803
% Ref 2: Emanuel, 2004 (a.k.a "the Lilly paper"). "Tropical cyclone
% energetics and structure".

% Created: 15 Feb 2010, Dan Chavas

% Input: r_user [km], V_user [ms-1], TC_lat [deg]
% Output: 1) Full solution: r_full [km], V_full [ms-1]
%         2) Asymptotic solution: r_Lil [km], V_Lil [ms-1]
%         3) Error flag (if res too large): res_error [0 or 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
%clc

%function [r_full,V_full,r_0_full,res_error] = r0_calc(r_user,V_user,fcor,Cd_in,wrad);
function [r_0_full,res_error,i_error] = r0_calc(r_user,V_user,fcor,Cd_in,wrad);

%% User input
%r_user=200;  %[km]
%V_user=10;      %[ms-1]
%TC_lat=20; %[deg]

dr=100;        %[m]; radial increments
r_0_i=80000000;    %[m] INITIAL GUESS (needs to be much larger than the actual outer radius)
iter_thresh=1000;
%%%%%%%%%%%%%%%

r_user=r_user*1000; %[km] --> [m]

%% Constants
w_rad=wrad; %.016 -- old value; %[ms-1]; equilibrium subsidence velocity
omega=7.292e-5; %[rad s-1]
f=fcor; %abs(2*omega*sind(TC_lat)); %[s-1]; Coriolis parameter
%C_D=1.5e-3;   %drag coefficient
C_D=Cd_in;

%% Numerical solution to full Lilly solution
%initialization
r_0=r_0_i;    %start with initial guess
r_0_new=r_0_i;  %start with initial guess
res=1e10;   %initial residual (simply needs to be some large number)
n_iter=0;   %initialize number of iterations to zero
r_0_save=r_0;

while(abs(res)>100*dr && n_iter<iter_thresh)  %repeat loop until solution converges sufficiently
    n_iter=n_iter+1;
    
    r_0_old=r_0;    %save old one
    r_0=r_0_new;    %set new one

    r2=r_user:dr:r_0;    %[m]; r>=r_user
    V=0*r2;  %initialize V
    V(1)=V_user;    %[ms-1]

    %calculate radial profile for given r_0
    for i=1:length(r2)-1
        a=2*C_D*r2(i)^2/(w_rad*(r_0^2-r2(i)^2));
        b=-1;
        c1=-r2(i)*f;

        dVdr=(a*V(i)^2+b*V(i)+c1)/r2(i);
        V(i+1)=V(i)+dVdr*dr;
    end

    indices=find(V<=0); %want index where V crosses zero
    if(numel(indices)>0)
        r_0_prof=r2(indices(1));
        r_0_new=(r_0+r_0_prof)/2;
        res=r_0-r_0_prof;
        r_0_save = r_0; %keep the last value for which an r_0 was found
    else    %jumped too far, no values with V<0
        r_0_new=(r_0+r_0_save)/2;   %must move downwards from above
        %keep residual same as before
    end
    
end
r_0_full=r_0;   %[m]

%{
%extend final profile inwards from r_user (simply repeats whats in while loop)
r_0_old=r_0;    %save old one
r_0=r_0_new;    %set new one

r1=r_user/2:dr:r_user;  %[m]; r<=r_user
V_inner=0*r1;
V_inner(end)=V_user;    %[ms-1]

%calculate radial profile for given r_0 (start at r_user moving INWARDS)
for i=1:length(r1)-1
    a=2*C_D*r1(end+1-i)^2/(w_rad*(r_0^2-r1(end+1-i)^2));
    b=-1;
    c1=-r1(end+1-i)*f;

    dVdr=(a*V_inner(end+1-i)^2+b*V_inner(end+1-i)+c1)/r1(end+1-i);
    V_inner(end-i)=V_inner(end+1-i)-dVdr*dr;
end

%FINAL radial profile
r_full=[r1(1:end-1) r2];
V_full=[V_inner(1:end-1) V];
%}
%double check everything
% r=r_full;
% V=V_full;
% dVdr2=(V(2:end)-V(1:end-1))/dr;
% dVdr2=[dVdr2(1) dVdr2];
% term1=2*C_D.*r.^2.*V.^2./(w_rad.*(r_0_full^2-r.^2));
% term2=-f.*r;
% term3=-V;
% term4=r.*dVdr2;
% residual=2*C_D.*r.^2.*V.^2./(w_rad.*(r_0_full^2-r.^2))-f.*r-V-r.*dVdr2;
res_error=0;
if(res > .1*r_0)    %give error message of the residual is approaching the same order as the solution
    sprintf('WARNING: Large Residual on Lilly numerical solution')
    res_error=1;    %saved in final output All_data.mat

i_error = 1;
else
i_error = 0;
end

%convert [m] --> [km]
%r_full=r_full/1000;

%% Simple Lilly profile (assume d(rV)/dr small)
%{
%%1) r_0, calculated from above radius of V_user (eqn 8, Dean et al 2009)
r_0_Lil=(r_user)^2+2*C_D*(r_user)*V_user^2/(f*w_rad);
r_0_Lil=sqrt(r_0_Lil); %[m]

dr=(r_0_Lil-r_user)/200;
r_Lil=r_user/2:dr:r_0_Lil;    %[m]

%calculate radial profile
V_Lil=sqrt((f*w_rad/(2*C_D))*(r_0_Lil^2-r_Lil.^2)./r_Lil);  %[ms-1]
r_Lil=[r_Lil r_0_Lil];  %ensure that V=0 actually shows up on plot
V_Lil=[V_Lil 0];

%convert [m] --> [km]
r_Lil=r_Lil/1000;
%}

%% Plot
%{
figure(1)
%subplot(2,2,4)
hold off
h=plot(r_Lil/1000,V_Lil,'g');
set(h,'linewidth',2)
set(gca,'fontweight','bold')
axis([r_user/2000 max(r_0_Lil/1000,r_0_full/1000) 0 30])
hold on
plot(r/1000,V,'r')
legend('simple','full')
input_title1=sprintf('Lilly outer wind profile comparison');
input_title2=sprintf('V_{user}=%i ms^{-1}; r(V_{user})=%2.2f km; TC_{lat}=%2.2f deg',V_user,r_user/1000,TC_lat);
title({input_title1,'',input_title2})
xlabel('radial distance [km]')
ylabel('azimuthal wind speed [m s^{-1}]')
scatter(r_user/1000,V_user,'b*')
%}

end
