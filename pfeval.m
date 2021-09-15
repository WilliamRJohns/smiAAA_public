function [pf] = pfeval(ZZ, pj, rj, polypart)
%Evaluates the partial fraction + polynomail representation of a rational
%function over ZZ where pj = the poles : rj = residues : polypart = the polynomial

%Cauchy Matrix
CC = (1./bsxfun(@minus, ZZ, pj)).';

%Evaluate our function
pf = CC*rj;
polyfren = polyval(polypart,ZZ);
pf = pf.'+polyfren;
end