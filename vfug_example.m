%Generates the plot and data for Vector Fitting Users Guide Example
% Adapted from VF example as below
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
clear all


%N=55; %order of approximation
%Niter1=5; %Fitting column sum: n.o. iterations
%Niter2=3; %Fitting column: n.o. iterations

disp('Reading data from file ...') %--> s(1,Ns), bigY(Nc,Nc,Ns)
fid1=fopen('fdne.txt','r');
Nc=fscanf(fid1,'%f',1);
Ns=fscanf(fid1,'%f',1);
bigY=zeros(Nc,Nc,Ns); s=zeros(1,Ns);
for k=1:Ns
  s(k)=fscanf(fid1,'%e',1);
  for row=1:Nc
    for col=1:Nc
      dum1=fscanf(fid1,'%e',1);
      dum2=fscanf(fid1,'%e',1);   
      bigY(row,col,k)=dum1+j*dum2;
    end
  end
end
s=i*s;
fclose(fid1);

tell=0;
for col=1:Nc
  for row=col:Nc
    tell=tell+1;
    f(tell,:)=squeeze(bigY(row,col,:)).'; %stacking elements into a single vector
  end
end
freq=s/(2*1i*pi);

%fit the data with aaa
%reoganize the data
k=1; F=[];
for ii=1:6
    for j=1:6
        F(k,:)=bigY(ii,j,:);
        k=k+1;
    end
end
%remove dupliate values(Why do these exists?)
[s,ia,~]=unique(s);
F=F(:,ia);

%Convergence VS poles comparison
polecompare(f,s);

%FastAAA demonstration
FastAAAcompare(f,s,1e-6,false,0);

%symetric smiaaa on the symetric problem
[symaaal,pwj,symaaa,pzj,~,pfj] = symmetricsmiaaah2(f,s,1e-4,false,25,1);
nn=length(pwj)/2;
[ppoles_aaa,~,ppfaaaf,~,~]=properrational(pzj.',pwj(nn+1:end),pwj(1:nn),pfj.',f,s);
[rmse_laaa,~,H2_laaa]=comp_error(f,symaaal);
[rmse_aaa,~,H2_aaa]=comp_error(f,symaaa);
disp('smiAAA symetric runtime');
mm = @() symmetricsmiaaah2(f,s,1e-4,false,1,1); % handle to function
timeit(mm)



disp('smiAAAl symetric runtime');
mm = @() symmetricsmiaaah2(f,s,1e-4,false,25,1); % handle to function
timeit(mm)

%Symetric Problem With RKfit real option
N=2*length(s);
Amat = util_build_real_matrix(s);
for j = 1:21
  Fmat{j} = util_build_real_matrix(f(j,:).');
end
tic
m   = length(ppoles_aaa);
tol = 1e-9;
b   = zeros(N, 1);
b(1:2:end) = 1;

init = util_ieee_poles(-2, 3, 1e-2, m);
[xi, ratfun, misfit, out] = rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     20, ...
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
timeit(mm)

figure(); hold on; rkmat=[];
for i=1:21
rk=ratfun{i};
%rk_res=residue(rk);
%figure();scatter(real(rk_res),imag(rk_res));
loglog(freq,abs(rk(s)-f(i,:)),'r','Linewidth',1.5);
rkmat=[rkmat;rk(s)];
%semilogy(freq,abs(rk(s)-F(:,i)),'rx');
end
[rmse_rk,rel_rk,H2_rk]=comp_error(f,rkmat);
Linf_rk_ieee=max(abs(rkmat-f),[],'all');

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
                
figure(); hold on; inf_rkmat=[];

disp('RKfit stable/inf poles runtime')
mm = @() rkfit(Fmat, Amat, b, init, ...
	         struct('maxit',     10, ...
				    'tol',       tol, ...
				    'real',      1, ...
				    'k',         0, ...
				    'safe',      1, ...
				    'stable',    0, ...
				    'reduction', 1)); % handle to function
timeit(mm)
inf_rkmat=[];
for i=1:21
rk=inf_ratfun{i};
%rk_res=residue(rk);
%figure();scatter(real(rk_res),imag(rk_res));
loglog(freq,abs(rk(s)-f(i,:)),'r','Linewidth',1.5);
inf_rkmat=[inf_rkmat;rk(s)];
%semilogy(freq,abs(rk(s)-F(:,i)),'rx');
end
[rmse_rk_inf,rel_rk_inf,H2_rk_inf]=comp_error(f,inf_rkmat);



%Parameters for VF
%Complex starting poles :
N=length(ppoles_aaa); %order of approximation
Niter1=5; %Fitting column sum: n.o. iterations
Niter2=3; %Fitting column: n.o. iterations

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
%opts.weightparam=1; %william test
%opts.weight=[];
  

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
timeit(mm)
vfmat=[];
for ii=1:6
    for jj=1:6
        bob=bigYfit(ii,jj,:);
        vfmat=[vfmat ; bob(:).'];
    end
end
[rmse_vf,rel_vf,H2_vf]=comp_error(F,vfmat);

figure()
semilogy(freq,max(abs(f-symaaa)+10^-13,[],1),'g','LineWidth',2.0);hold on;
semilogy(freq,max(abs(f-symaaal),[],1),'b','LineWidth',2.0);
semilogy(freq,max(abs(F-vfmat),[],1),'c','LineWidth',2.0);
semilogy(freq,max(abs(f-rkmat),[],1),'r','LineWidth',2.0);
semilogy(freq,max(abs(f-inf_rkmat),[],1),'k','LineWidth',2.0);

title('Log Errors VF User''s Guide example')
ylabel('Log error')
xlabel('Frequency')
legend('smiAAA','smiAAA-L','vf','RKfit(ieee poles)','RKfit(inf poles)')

%Linf errors
Linf_aaa=max(abs(symaaa-f),[],'all');
Linf_aaal=max(abs(symaaal-f),[],'all');
Linf_rk_inf=max(abs(inf_rkmat-f),[],'all');
Linf_vf=max(abs(vfmat-F),[],'all');


