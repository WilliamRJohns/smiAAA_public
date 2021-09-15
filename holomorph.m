%Define some function
clear all
close all
format long
alf=1e-3;
%f1 = @(x) 14./((x-1-.1*alf)*(x-1-.01*alf)*(x-1-.00001*alf))+13*x+3 %1 pole appears in the domain with 1e-6 aaa tol, no eigs in the domain

%alf=1e-6;
%A1=[-1 0 0;0 -1 0;0 0 -1];
%f1 = @(x) 14./((x-1-.1*alf)*(x-1-.001*alf)*(x-1-.000001*alf))+13*x+3 %2 pole appears in the domain with 1e-6 aaa tol, 1 eig in the domain

A1=[-1 0 0;0 -1 0;0 0 -1];
f1 = @(x) 14./((x-1-.1*alf).*(x-1-.01*alf))+1100*x+3 %2 pole appears in the domain with 1e-6 aaa tol, 1 eig in the domain
Flam = @(x) f1(x)*A1;
%Zint = rand(1,100).*exp(rand(1,100)*2*pi*1i);
%Zbnd = 1*exp(1i*linspace(0,2*pi,2*100));
%dom=union(Zint,Zbnd);
load dom.mat;
F=f1(dom);

[laaaf,wj,aaaf,zj,~,fj]=miaaa(F(:).',dom(:).',1e-6,false,40);
nn=length(wj)/2;
[ppoles_aaa,~,ppfaaaf,~,~]=properrational(zj.',wj(nn+1:end),wj(1:nn),fj.',F(:).',dom(:).');

figure()
%add a circle
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
h = plot(xunit, yunit);
hold on

int_poles=ppoles_aaa(real(ppoles_aaa).^2+imag(ppoles_aaa).^2 < 1);
ext_poles=ppoles_aaa(real(ppoles_aaa).^2+imag(ppoles_aaa).^2 > 1);
scatter(real(int_poles),imag(int_poles),'rx');
scatter(real(ext_poles),imag(ext_poles),'gx');


%add the poles of f1
scatter([1+.1*alf,1+0.1*alf,1+.01*alf],[0 0 0],'bx')

%Now lets slap a matrix on there

A1=A1./norm(A1,'fro');
coefnorms={A1};
%Compute the Linearization
[Am, Bm] = liAAALinearize(zj.', fj.', wj.', coefnorms, 0);
[Vw,Dw]=eig(Am,Bm);
evs=diag(Dw);
evecs=Vw(1:3,:);

idx=find(sqrt(real(evs).^2+imag(evs).^2)<1);
evs=evs(idx);
evecs=evecs(:,idx);

%figure();
scatter(real(evs),imag(evs),'k*');
legend('','Poles in D','Poles !in D','Poles F1','Eigs')

[claaaf,cwj,caaaf,czj,~,cfj]=circlemiaaa(F(:).',dom(:).',1e-6,false,40);
nn=length(cwj)/2;
[cppoles_aaa,~,cppfaaaf,~,~]=properrational(czj.',cwj(nn+1:end),cwj(1:nn),cfj.',F(:).',dom(:).');

figure();
h = plot(xunit, yunit);hold on;
scatter(real(cppoles_aaa),imag(cppoles_aaa),'bo');
scatter(real(ppoles_aaa),imag(ppoles_aaa),'rx');
legend('Unit Circle','sv-aaa poles','smiAAA poles');
title('Pole Comparison sv-aaa vs smiAAA ')

%Compute the Linearization with circlemiaaa
[cAm, cBm] = liAAALinearize(czj.', cfj.', cwj.', coefnorms, 0);
[cVw,cDw]=eig(cAm,cBm);
cevs=diag(cDw);
cevecs=cVw(1:3,:);

cidx=find(sqrt(real(cevs).^2+imag(cevs).^2)<1);
cevs=cevs(cidx);
cevecs=cevecs(:,cidx);

figure();
h = plot(xunit, yunit);hold on;
scatter(real(cevs),imag(cevs),'bo');
scatter(real(evs),imag(evs),'rx');
legend('Unit Circle Boundry','sv-aaa eigs','ciclemiaa eigs')


normF = 0; 
for J =1:length(dom)
    normF = max(normF, norm(Flam(dom(J))));
end
Nevs = length(evs);
resids = zeros(Nevs, 1);
for kk = 1:Nevs
     nF = max(norm(Flam(evs(kk))),normF);
     resids(kk) =norm(Flam(evs(kk))*evecs(:,kk))/(nF*norm(evecs(:,kk))); %backward error for eigenpair (evs(kk), evecs(kk))
end
beidx=find(resids<=1e-11);
evs=evs(beidx);
evecs=evecs(:,beidx);

%Moebius Transform stuff
m = @(z) -2*z/(z+1i);
mdom=m(dom);



function [LE, LF] = liAAALinearize(z, f, lw, coeffs, sparseFlag)
%William's test linearization for lawson optimized solutions
n = size(coeffs{1},1);
[m,d] = size(f);
% If it is large and not sparse, we make it sparse nonetheless to avoid
% swapping
if (m+1)*n > 1000
    sparseFlag = 1;
end

if sparseFlag
    F = -speye(m) + spdiags(ones(m-1,1),-1, m,m);
    E = -spdiags(z,0,m,m);
    if m > 1
        % In this way we cover the case when m = 1, where the second
        % addendum is 0. The dense case does not need this distinction.
        E = E + spdiags(z(1:m-1),-1,m,m);
    end
    In = speye(n);
    b = spalloc(m,1,1); b(1) = 1;
    LE = sparse(n*(m+1), n*(m+1));
    LF = LE;
else
    F = -eye(m) + diag(ones(1,m-1),-1);
    E = -diag(z) + diag(z(1:m-1), -1);
    In = eye(n);
    b = zeros(m,1); b(1) = 1;
    LE = zeros(n*(m+1), n*(m+1));
    LF = LE;
end
F(1,1) = 0;
E(1,:) = lw(1:length(lw)/2);


a = f.*lw(length(lw)/2+1:end); % the columns of a are the vectors a_i: f_i(z_j)*w_j
% Now build Etilde and Ftilde
Et = [-b E];
Ft = [zeros(m,1) F];

aux = kron(a(:,1).', coeffs{1});
for j = 2:d
    aux = aux + kron(a(:,j).', coeffs{j});
end
LE(1:n,n+1:end) = aux;
LE(n+1:end,:) = kron(Et, In);
LF(n+1:end,:) = kron(Ft, In);

end

