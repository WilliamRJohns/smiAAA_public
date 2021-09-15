% Thesis example with mismatch num and den degree where przd deflates 3
% extraneous eigenvalues.
%Requires Chebfun toolbox to be installed, http://www.chebfun.org/
format long
vZ =1i* linspace(-10.,10.,100);
vZ2 = linspace(-10.,10.,2000);
% Define denominator and find its zeros (poles of f)
fden = @(z) (-3 +2.*z  + z.^2) .* (1+z.^2);
syms x;
dd = sym2poly(fden(x));
fpoles = roots(dd);
datpoles = polyval(polyder(dd), fpoles);
% Modify numerator to have an exact order (pow+2,4) 
    fnum = @(z) (1.23 + z) .* (1 + z).*(2+z).*(5+z).*(8+z).^3;
    f  = @(z)  fnum(z)./fden(z);
    % Compute zeros and residues of fin
    c = sym2poly(fnum(x));
    fzeros = roots(c);
    fres = fnum(fpoles)./datpoles;
    % Run aaa 
    [r0,pol0,res0,zer0,z0,f0,w0,errvec0] = aaa(f,vZ, 'tol', 10^-11,'mmax', 30, 'lawson',0);
    %przd routine
    [newpoles]=przd(z0,w0);
    % Results
    disp('PRZ poles')
    disp(pol0)
    disp('PRZD poles')
    disp(newpoles)
    
 
