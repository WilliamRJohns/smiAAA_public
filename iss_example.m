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
disp('smiAAA symetric runtime');
mm = @() symmetricsmiaaah2(f,s,1e-4,false,1,1); % handle to function
%timeit(mm)

disp('smiAAAl symetric runtime');
mm = @() symmetricsmiaaah2(f,s,1e-4,false,25,1); % handle to function
%timeit(mm)

%Symetric Problem With RKfit real option
Amat = util_build_real_matrix(1i*w);
for j = 1:ell
  Fmat{j} = util_build_real_matrix(F(:, j));
end

m   = length(ppoles_aaa);
tol = 1e-4;
b   = zeros(N, 1);
b(1:2:end) = 1;

init = util_ieee_poles(-2, 3, 1e-2, m);
%init=ppoles_aaa;

[xi, ratfun, misfit, out] = rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     10, ...
				    'tol',       tol, ...
				    'real',      1, ...
				    'k',         0, ...
				    'safe',      1, ...
				    'stable',    0, ...
				    'reduction', 1));

                            
disp('RKfit stable/ieee poles runtime')
mm = @() rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     10, ...
				    'tol',       tol, ...
				    'real',      1, ...
				    'k',         0, ...
				    'safe',      1, ...
				    'stable',    0, ...
				    'reduction', 1)); % handle to function
%timeit(mm)

figure(); hold on; rkmat=[];
for i=1:ell
rk=ratfun{i};
%rk_res=residue(rk);
%figure();scatter(real(rk_res),imag(rk_res));
loglog(freq,abs(rk(s)-f(i,:)),'r','Linewidth',1.5);
rkmat=[rkmat;rk(s)];
%semilogy(freq,abs(rk(s)-F(:,i)),'rx');
end
[rmse_rk,rel_rk,H2_rk]=comp_error(f,rkmat);

%Covert rkfit to partial fractions and try Lawson
clear poles; %rkfit poles function uses this name
rkpoles=poles(ratfun{1});
rkres=[];
for ii=1:6
    bob=ratfun{ii};
    rkres=[rkres residue(bob)];
end

init = inf*ones(1,m);
[xi, inf_ratfun, misfit, out] = rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     20, ...
				    'tol',       tol, ...
				    'real',      1, ...
				    'k',         0, ...
				    'safe',      1, ...
				    'stable',    0, ...
				    'reduction', 1));
                
disp('RKfit stable/inf poles runtime')
mm = @() rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     10, ...
				    'tol',       tol, ...
				    'real',      1, ...
				    'k',         0, ...
				    'safe',      1, ...
				    'stable',    0, ...
				    'reduction', 1)); % handle to function
%timeit(mm)
                
figure(); hold on; inf_rkmat=[];
for i=1:ell
rk=inf_ratfun{i};
%rk_res=residue(rk);
%figure();scatter(real(rk_res),imag(rk_res));
loglog(freq,abs(rk(s)-f(i,:)),'r','Linewidth',1.5);
inf_rkmat=[inf_rkmat;rk(s)];
%semilogy(freq,abs(rk(s)-F(:,i)),'rx');
end
[rmse_rk_inf,rel_rk_inf,H2_rk_inf]=comp_error(f,inf_rkmat);



figure()
loglog(freq,max(abs(f-symaaal)+10^-13,[],1),'b','LineWidth',2.0);hold on;
loglog(freq,max(abs(f-symaaa)+10^-13,[],1),'g','LineWidth',2.0);
loglog(freq,max(abs(f-rkmat)+10^-13,[],1),'r','LineWidth',2.0);
title('Log Errors symetric smiAAA')
xlabel('Log error')


%Parameters for VF
%Complex starting poles :
N=length(ppoles_aaa); %order of approximation
Niter1=5; %Fitting column sum: n.o. iterations
Niter2=3; %Fitting column: n.o. iterations
Ns=length(w);
Nc=ell;


w=s/i;
bet=linspace(w(1),w(Ns),N/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
end


%weight=ones(1,Ns);
%weight=1./abs(f);
weight=1./sqrt(abs(f));

%Fitting options
opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=3;      %Fitting includes D and E (%William Changed)
opts.spy1=0; 
opts.spy2=1; 
opts.logx=0; 
opts.logy=1; 
opts.errplot=1;
opts.phaseplot=1;

opts.skip_pole=0; 
opts.skip_res=1;
opts.cmplx_ss=1;  %=1 --> Will generate state space model with diagonal A
opts.legend=1;
opts.weightparam=1; %william test
opts.weight=[];

%make thing for vfdriver
count=1;
for ik=1:3
    for ki=1:3
        bigY(ik,ki,:)=F(:,count);
        count=count+1;
    end
end
  

%Forming (weighted) column sum:
g=0;
for n=1:Nc
  %g=g+f(n,:); %unweighted sum     
  g=g+f(n,:)/norm(f(n,:));
  %g=g+f(n,:)/sqrt(norm(f(n,:)));     
end
weight_g=1./abs(g);
weight=ones(1,Ns);
%Fully Symetric Problem for direct comparison with VF
[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts);
disp('VFdriver runtime');
mm = @() VFdriver(bigY,s,poles,opts); % handle to function
%timeit(mm)
vfmat=[];
for ii=1:3
    for jj=1:3
        bob=bigYfit(ii,jj,:);
        vfmat=[vfmat ; bob(:).'];
    end
end
[rmse_vf,rel_vf,H2_vf]=comp_error(f,vfmat);

figure()
semilogy(freq,max(abs(f-symaaa)+10^-13,[],1),'g','LineWidth',2.0);hold on;
semilogy(freq,max(abs(f-symaaal),[],1),'b','LineWidth',2.0);
semilogy(freq,max(abs(f-vfmat),[],1),'c','LineWidth',2.0);
semilogy(freq,max(abs(f-rkmat),[],1),'r','LineWidth',2.0);
semilogy(freq,max(abs(f-inf_rkmat),[],1),'k','LineWidth',2.0);

title('Log Errors symetric problem comparison')
ylabel('Log error')
xlabel('Frequency')
legend('smiAAA','smiAAA-L','vf','RKfit(ieee poles)','RKfit(inf poles)')

%Linf errors
Linf_aaa=max(abs(symaaa-f),[],'all');
Linf_aaal=max(abs(symaaal-f),[],'all');
Linf_rk=max(abs(rkmat-f),[],'all');
Linf_vf=max(abs(vfmat-f),[],'all');
Linf_rk_inf=max(abs(inf_rkmat-f),[],'all');

figure()
loglog(freq,max(abs(f-symaaa)+10^-13,[],1),'g','LineWidth',2.0);hold on;
loglog(freq,max(abs(f-symaaal),[],1),'b','LineWidth',2.0);
loglog(freq,max(abs(f-vfmat),[],1),'c','LineWidth',2.0);
loglog(freq,max(abs(f-rkmat),[],1),'r','LineWidth',2.0);
loglog(freq,max(abs(f-inf_rkmat),[],1),'k','LineWidth',2.0);

title('Log Errors ISS example')
ylabel('Log error')
xlabel('Frequency')
legend('smiAAA','smiAAA-L','vf','RKfit(ieee poles)','RKfit(inf poles)')

figure()
loglog(freq,abs(f-symaaa)+10^-13,'b','LineWidth',2.0);hold on;
loglog(freq,abs(f-inf_rkmat)+10^-13,'r','LineWidth',2.0);
figure()
loglog(freq,abs(f-inf_rkmat)+10^-13,'r','LineWidth',2.0);

