%MPI_nodrag.m

%Calculate the theoretical MPI assuming no contribution of the gustiness
%wind to surface drag.  This is calculated from the normal MPI.

%Created: 27 Jul 2012, Dan Chavas

%Theory: Vb^2 = Vb_drag^2*(1+usfc/Vb)  -- from basic MPI theory, assuming
%usfc is included in the calculation of enthalpy flux but not surface drag

%clear
%clc
%%TESTING
%{
Vp_nominal = 60;
usfc = 3;
%}
function [Vp_usfcnodrag] = MPI_nodrag(Vp_nominal,usfc)

%polynomial coefficients: Vp^3 - Vp_nominal^2*V_p - u*Vp_nominal^2 = 0
coefs = [1 0 -Vp_nominal^2 -usfc*Vp_nominal^2];
theroots = roots(coefs);
Vp_usfcnodrag = theroots(theroots>Vp_nominal);

end