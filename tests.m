%%Try International Space Station example from 
%%http://slicot.org/20-site/126-benchmark-examples-for-model-reduction
%%http://guettel.com/rktoolbox/examples/html/example_iss.html
%%
%%Generates data and plots for Internation Space Station Example in thesis;
%%Requires 'rktoolbox' to be installed. http://guettel.com/rktoolbox/
%%Requires 'matrix fitting toolbox' to be installed. https://www.sintef.no/projectweb/vectorfitting/downloads/matrix-fitting-toolbox/
close all;clear all;
load('iss.mat');

ell = 3*3;
N   = 2*length(w);
F   = zeros(N/2, ell);

for j = 1:N/2
  % We now evaluate the responses at 1i*w(j) from the state
  % space representation. The responses at -1i*w(j) is the
  % conjugated one.
  resp = full(C*((A-1i*w(j)*speye(length(A)))\B));
  F(j, :) = resp(:).';
end
f=F.';
s=1i*w.';
freq=s/(2*1i);

%poles comparison with FastAAA
%polecompare(f,s);

%FastAAA demonstration
%FastAAAcompare(f,s,1e-4,false,0);

% %Stable approximation of multiple functions with common poles can also be
% %constructed similar to vector fitting (this does not include the symetry
% %of vector fitting)
% [saaafl,wj,saaaf,zj,wj1,fj,err] = smiaaa(f,s,1e-6,false,10,1);
% wden = wj(1:length(zj)); wnum = wj(length(zj)+1:end);
% [spoles_aaa,sres_aaa,pfaaaf,~,bestpoly]=properrational(zj.',wnum,wden,fj.',saaaf,s);
% disp('smiAAA 9 functions Runtime')
% mm = @() smiaaa(f,s,1e-6,false,0,1); % handle to function
% timeit(mm)
% 
% figure();
% loglog(freq,abs(f-saaaf)+1e-13)
% title('smiAAA tol=1e-6 9 iss functions')
% 
% figure()
% semilogx(freq,max(abs(f-saaafl),[],1)+10^-13,'b',freq,max(abs(f-saaaf),[],1),'r','Linewidth',1.5)
% title('smiAAA-L Errors tol=10^{-6}')
% xlabel('Frequency Hz')
% ylabel('Max Abs(Error)')
% legend('smiAAA-L','smiAAA')

%symetric smiaaa on the symetric problem
[symaaal,pwj,symaaa,pzj,~,pfj] = symmetricsmiaaah2(f,s,1e-4,false,25,1);
nn=length(pwj)/2;
[ppoles_aaa,~,ppfaaaf,~,~]=properrational(pzj.',pwj(nn+1:end),pwj(1:nn),pfj.',f,s);
[rmse_laaa,~,H2_aaal]=comp_error(f,symaaal);
[rmse_aaa,~,H2_aaa]=comp_error(f,symaaa);

[nsymaaal,npwj,nsymaaa,npzj,~,npfj] = symmetricsmiaaah2o(f,s,1e-4,false,25,1);
nn=length(npwj)/2;
[nppoles_aaa,~,nppfaaaf,~,~]=properrational(npzj.',npwj(nn+1:end),npwj(1:nn),npfj.',f,s);
[rmse_nlaaa,~,H2_naaal]=comp_error(f,nsymaaal);
[rmse_naaa,~,H2_naaa]=comp_error(f,nsymaaa);

