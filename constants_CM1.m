%constants_CM1.m



function [c_CM1] = constants_CM1();

%from CM1 'include/constants.incl':
g      = 9.81;
to     = 273.15;
rd     = 287.04;
rv     = 461.5;
cp     = 1005.7;
cpinv  = 1.0/cp;
cv     = cp-rd;
cpv    = 1870.0;
cvv    = cpv-rv;
p00    = 1.0e5;
rp00   = 1.0/p00;
rcp    = 1.0/cp;
pi     = 3.14159265358979323;
cpdcv  = cp/cv;
rovcp  = rd/cp;
rddcp  = rd/cp;
rddcv  = rd/cv;
rddrv  = rd/rv;
cvdrd  = cv/rd;
cpdrd  = cp/rd;
eps    = rd/rv;
reps   = rv/rd;
repsm1 = rv/rd-1.0;
cpt    = (cpv/cp)-1.0;
cvt    = (cvv/cv)-1.0;
pnum   = (cp*rv)-(cpv*rd);
xlv    = 2501000.0;
lathv  = xlv;
xls    = 2834000.0;
lvdcp  = xlv/cp;
condc  = xlv*xlv/(rv*cp);
cpl    = 4190.0;
cpi    = 2106.0;
lv1    = xlv+(cpl-cpv)*to;
lv2    = cpl-cpv;
ls1    = xls+(cpi-cpv)*to;
ls2    = cpi-cpv;
govtwo = 0.5*g;
piddeg = pi/180.0;
degdpi = 180.0/pi;
clwsat = 1.e-6;
rhow   = 1.0e3;

%c_CM1: [g rd cp cv p00 xlv]
c_CM1(1) = g;
c_CM1(2) = rd;
c_CM1(3) = cp;
c_CM1(4) = cv;
c_CM1(5) = p00;
c_CM1(6) = xlv;


end