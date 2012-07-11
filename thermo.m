%thermo.m

%Purpose: calculate some relevant thermodynamic quantities

%Created: 20 May 2011, Dan Chavas


function [T rho rh the s s_sat gam_m thv] = thermo(p,th,qv,ql)


%constants
Rd = 287;   %[J/K/kg]
Cpd = 1004; %[J/K/kg]
Cpv = 1952; %[J/K/kg]
Rv = 461;   %[J/K/kg]
eps = Rd/Rv;
p0 = 100000;    %[Pa]
Lv = 2.5e6;    %[J/kg]
g = 9.8067; %[m/s^2]

%%calculate actual temperature from potential temperature
T = th.*(p/p0).^(Rd/Cpd);

%%calculate virtual temperature from actual temperature and wv mix ratio
Tv = T.*(1+.608*qv);

%%calculate density from pressure and virtual temperature
rho = p./(Rd*Tv);

%%calculate saturation vapor pressures
T_C = T-273.15;
es = 611.2*exp((17.67*T_C)./(T_C+243.5));  %[Pa]

%%calculate saturation mixing ratios
qvs = eps*(es./(p-es));  %[kg/kg]

%%calculate relative humidity
rh = qv./qvs;

%%calculate equivalent potential temperature
Cp = Cpd + qv*Cpv;
the=th.*exp((Lv.*qv)./(Cp.*T));
%the2(i,j)=thet(i,j).*exp(LO.*q(i,j)./1000/(CP*(A*log(q(i,j)./1000)+B))); (from Brian Tang's ASPECH model)

%%calculate virtual potential temperature (BR09b, p.3)
thv=th.*(1+(Rv/Rd)*qv)./(1+qv+ql);

%%calculate entropy (eqn (1) Emanuel11)
s = Cpd*log(T) - Rd*log(p) + Lv*qv./T;

%%calculate saturation entropy (eqn (1) Emanuel11)
s_sat = Cpd*log(T) - Rd*log(p) + Lv*qvs./T;

%%calculate moist adiabatic lapse rate (see eqn 4.7.3 Emanuel Convection
%%book): ignoring liquid water content contribution
gam_d = (g/Cpd)*(1+qv)./(1+qv*(Cpv/Cpd));
gam_m = gam_d.*(1+(Lv*qv)./(Rd.*T))./(1+Lv^2*qv.*(1+qv/eps)./(Rv*T.^2.*(Cpd+qv*Cpv)));

end