%vmax_RE.m

%Purpose: calculate the true vmax and rmax given an input vmax and rmax
%from RE87

%Created: 22 Mar 2011, Dan Chavas

clc
clear
%clf

%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vmaxs=[2.5 5 7.5 10 12.5 15];    %[ms-1]; desired true vmax; v(rmax)=vmax
r0s=[100 200 400 600 800 1000 1200 1600 2000]*1000;    %[m]
fcor=5e-5;  %[s-1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:length(vmaxs)
    for j=1:length(r0s)

        vmax=vmaxs(i);
        r0=r0s(j);
        rmax=r0/5; %assume rmax scales with r0
        rm_REs=5*1000:.01*1000:r0/2;
        
        %% Solve for RE87 input vm_RE and rm_RE that will give desired true vmax and rmax

        %%Step 1: vm_RE as a function of rm_REs
        %from RE87: v(rmax)=vmax
        dd1=2*rm_REs./(rmax+rm_REs);
        dd2=2*rm_REs./(r0+rm_REs);
        vm_REs=sqrt(((vmax+.5*fcor*rmax)^2-.25*fcor^2*rmax^2)./((rmax./rm_REs).^2.*(dd1.^3-dd2.^3)));

        %%%%Step 2: dv/dr = 0 at r=rmax
        rref=rmax;
        %0 = - fcor/2 - ((24*rm_REs*rref^2*vm_REs^2)/(rm_REs + rref)^4 - (fcor^2*rref)/2 + (2*rref*vm_REs^2*((8*rm_REs^3)/(r0 + rm_REs)^3 - (8*rm_REs^3)/(rm_REs + rref)^3))/rm_REs^2)/(2*((fcor^2*rref^2)/4 - (rref^2*vm_REs^2*((8*rm_REs^3)/(r0 + rm_REs)^3 - (8*rm_REs^3)/(rm_REs + rref)^3))/rm_REs^2)^(1/2))
        dvdr = - fcor/2 - ((24.*rm_REs.*rref.^2.*vm_REs.^2)./(rm_REs + rref).^4 - (fcor^2.*rref)./2 + (2.*rref.*vm_REs.^2.*((8*rm_REs.^3)./(r0 + rm_REs).^3 - (8*rm_REs.^3)./(rm_REs + rref).^3))./rm_REs.^2)./(2*((fcor^2.*rref.^2)/4 - (rref.^2.*vm_REs.^2.*((8*rm_REs.^3)./(r0 + rm_REs).^3 - (8*rm_REs.^3)./(rm_REs + rref).^3))./rm_REs.^2).^(1/2));

        %plot(rm_REs,dvdr)
        index=find(dvdr>0);
        rm_RE=rm_REs(index(1));
        vm_RE=vm_REs(index(1));
        %}

        %% Plot the RE87 profile for givein input vmax and rmax
        %vm_RE=15.6655;    %[ms-1]
        %rm_RE=39.1*1000; %[m]
        %fcor=5e-5;  %[s-1]
        %r0=400*1000;    %[m]
        
        rref=0:.1*1000:r0;

        %from RE87
        clear dd1 dd2 vr
        dd1=2*rm_RE./(rref+rm_RE);
        dd2=2*rm_RE./(r0+rm_RE);
        vr=sqrt(vm_RE^2.*(rref./rm_RE).^2.*(dd1.^3-dd2.^3)+.25*fcor^2.*rref.^2)-.5*fcor.*rref;

        %Check that it worked
%        plot(rref/1000,vr)
        vmax_real(i,j)=max(vr);
        rmax_real(i,j)=rref(vr==max(vr));
        
        %save RE87 inputs
        rm_REall(i,j)=rm_RE;
        vm_REall(i,j)=vm_RE;
        
    end
end

rmaxs=r0s/5;   %[m]; desired true rmax

vmaxs=vmaxs';

for i=1:length(vmaxs)
    for j=1:length(rmaxs)
        
        data_all((i-1)*length(rmaxs)+j,:)=[vmaxs(i) rmaxs(j)/1000 vm_REall(i,j) rm_REall(i,j)/1000];
        
    end 
end

data_all
%}
