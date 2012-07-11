%TC_selfsim.m


%%inputs
r0 = 100;    %*rmax
Ck=1;
Cd=1;


%%initialize radius array
r = 1:.1:r0;

%%calculate fraction of angular momentum lost from r0 to rmax (from RE11)
chi = (.5*Ck/Cd)^(1/(2-Ck/Cd));

V = (1-(1/chi)*(r/r0).^2)./r;

figure(4)
plot(r,V)
axis([0 19 0 1])